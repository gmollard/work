#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
gen_find_statement.py - READ-ONLY.  Given an action block's object id and a
"statement id", find out what that statement id corresponds to in the Gen
Toolset OLE API (OpenAPI / MetaModelObject).

  python gen_find_statement.py --block-id 12345 --statement-id 67890

It answers three questions:
  1. Is the statement id an object id?   (GetMetaModelObject + type code)
  2. If so, is that object reachable from the action block, and through
     which chain of associations?        (breadth-first walk)
  3. What are the block's direct associations (codes, counts, targets), so
     you can compare positions/ids if the statement id is NOT an object id.

Name lookups through pywin32 attribute access do not work on this Toolset,
so every call goes through GetIDsOfNames / Invoke by member id (class Late).
"""
import argparse
import sys
from collections import deque

import pythoncom
import win32com.client as wc

ANY = pythoncom.DISPATCH_METHOD | pythoncom.DISPATCH_PROPERTYGET
RC_INVALID_FOR_OBJECT = (16, 21)   # ATC_INVALID, ASSOC_NOT_FOUND
RC_NO_SUCH_LINK = 19               # code valid, link absent


class Late(object):
    """Late-bound COM wrapper that calls members by dispid."""

    def __init__(self, d):
        self.o = d._oleobj_ if hasattr(d, "_oleobj_") else d
        self.ids = {}

    def id(self, name):
        if name not in self.ids:
            r = self.o.GetIDsOfNames(name)
            self.ids[name] = r[0] if isinstance(r, (tuple, list)) else r
        return self.ids[name]

    def call(self, name, *args):
        r = self.o.Invoke(self.id(name), 0, ANY, 1, *args)
        return Late(r) if type(r).__name__ == "PyIDispatch" else r

    def rc(self):
        return self.call("LastReturnCode")


def items(coll):
    """Iterate a MetaModelObjects collection through its default member."""
    n = coll.call("Count")
    if not n:
        return

    def get(i):
        r = coll.o.Invoke(0, 0, ANY, 1, i)     # dispid 0 = default member
        return Late(r) if type(r).__name__ == "PyIDispatch" else r

    try:                                       # ASSUMPTION: 1-based like Models
        get(n)
        base = 1
    except Exception:
        base = 0
    for i in range(base, base + n):
        yield get(i)


def connect():
    ts = wc.Dispatch("COOLgen.Toolset")
    model = Late(ts.Models(1))
    return model.call("OpenAPI")


def get_obj(api, oid):
    o = api.call("GetMetaModelObject", oid)
    if api.rc() != 0 or not isinstance(o, Late):
        return None
    return o


class Namer(object):
    def __init__(self, api, mnemonic):
        self.api = api
        try:
            self.code = api.call("FetchSchemaPrpTypeByMnem", mnemonic)
            if api.rc() != 0:
                self.code = None
        except Exception:
            self.code = None

    def __call__(self, o):
        if self.code is None:
            return ""
        try:
            s = o.call("FetchCharArrayPrp", self.code)
            return s if o.rc() == 0 and isinstance(s, str) else ""
        except Exception:
            return ""


class Walker(object):
    def __init__(self, max_code):
        self.max_code = max_code
        self.valid = {}                       # (type, fetch) -> codes

    def edges(self, o):
        otype = o.call("ObjectTypeCode")
        for fetch in ("FetchCardManyAsc", "FetchCardOneAsc"):
            key = (otype, fetch)
            known = self.valid.get(key)
            probing = known is None
            codes = range(1, self.max_code + 1) if probing else sorted(known)
            found = set()
            for code in codes:
                try:
                    res = o.call(fetch, code)
                except Exception:
                    continue
                r = o.rc()
                if r in RC_INVALID_FOR_OBJECT:
                    continue
                if probing and r in (0, RC_NO_SUCH_LINK):
                    found.add(code)
                if r != 0 or res is None:
                    continue
                if fetch == "FetchCardManyAsc":
                    for t in items(res):
                        yield code, t
                else:
                    yield code, res
            if probing:
                self.valid[key] = found


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--block-id", type=int, required=True)
    ap.add_argument("--statement-id", type=int, required=True)
    ap.add_argument("--max-depth", type=int, default=3)
    ap.add_argument("--max-nodes", type=int, default=400)
    ap.add_argument("--max-code", type=int, default=1500)
    ap.add_argument("--name-mnemonic", default="NAME")   # ASSUMPTION
    args = ap.parse_args()

    api = connect()
    name = Namer(api, args.name_mnemonic)

    def desc(o):
        return "#%s type=%s '%s'" % (o.call("ID"), o.call("ObjectTypeCode"), name(o))

    block = get_obj(api, args.block_id)
    if block is None:
        sys.exit("Block id %s not found (rc=%s)." % (args.block_id, api.rc()))
    print("BLOCK    ", desc(block))

    # 1. Is the statement id an object id?
    stmt = get_obj(api, args.statement_id)
    if stmt is None:
        print("STATEMENT id %s is NOT an object id in this model (rc=%s)."
              % (args.statement_id, api.rc()))
        print("-> it is probably a sequence/line number or an id from another "
              "source. See the direct associations below.")
    else:
        print("STATEMENT", desc(stmt), "(exists as an object)")

    walker = Walker(args.max_code)

    # 2. Breadth-first walk from the block
    parent = {args.block_id: None}
    queue = deque([(block, 0)])
    direct = {}                      # assoc code -> list of (id, type, name)
    scanned = 0
    while queue and scanned < args.max_nodes:
        o, depth = queue.popleft()
        scanned += 1
        if depth >= args.max_depth:
            continue
        oid = o.call("ID")
        for code, t in walker.edges(o):
            tid = t.call("ID")
            if depth == 0:
                direct.setdefault(code, []).append(
                    (tid, t.call("ObjectTypeCode"), name(t)))
            if tid not in parent:
                parent[tid] = (oid, code)
                queue.append((t, depth + 1))

    if stmt is not None:
        if args.statement_id in parent:
            chain, cur = [], args.statement_id
            while parent[cur] is not None:
                pid, code = parent[cur]
                chain.append((pid, code, cur))
                cur = pid
            print("\nREACHABLE from the block:")
            for pid, code, cid in reversed(chain):
                print("  #%s --assoc %s--> #%s" % (pid, code, cid))
        else:
            print("\nNot reached within depth %d / %d nodes (scanned %d). "
                  "Raise --max-depth or --max-nodes, or the id belongs to "
                  "another block." % (args.max_depth, args.max_nodes, scanned))

    # 3. Direct associations of the block
    print("\nDIRECT associations of the block (assoc code: count -> first targets)")
    for code in sorted(direct):
        lst = direct[code]
        head = ", ".join("#%s(t%s %s)" % (i, t, n) for i, t, n in lst[:5])
        print("  %-5s %3d -> %s%s" % (code, len(lst), head,
                                      " ..." if len(lst) > 5 else ""))


if __name__ == "__main__":
    main()
