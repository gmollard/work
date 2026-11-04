#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
gen_inspect.py - READ-ONLY inspector for the CA Gen Toolset OLE API.
Prints the content of meta-model objects: every property that answers and
every association with its targets.  All calls go through
GetIDsOfNames/Invoke by member id (this Toolset does not support
pywin32's normal attribute access).

Commands
  types                      count objects per object-type code (overview)
  object --id N              full content of one object
  tree   --id N [--depth D]  content of N and everything reachable from it
  scan   --type T [--limit K]   content of the first K objects of type T
                                (T = type code, or a mnemonic such as ACBLKBSD)

Options (before the command)
  --mnemonics FILE   text file, one mnemonic per line (copy them from
                     odrpta.chm).  Each is resolved to its numeric code so
                     the output can show  57(ACBLKBSD)  instead of  57.
  --json OUT.json    also write everything printed to a JSON file.
  --max-prop N       highest property type code to probe   (default 800)
  --max-assoc N      highest association type code to probe (default 1500)

How discovery works: there is no call that lists an object's properties
or associations, so codes 1..N are probed once per object type and the
codes that answer are cached for the other objects of that type.
Properties are read with FetchCharArrayPrp / FetchSingleCharPrp /
FetchSingleNumericPrp; the label (text/char/num) shows which one answered.
A code answering under several labels needs a look at odrpta.chm.
"""
import argparse
import json
import sys
from collections import Counter, deque

import pythoncom
import win32com.client as wc

ANY = pythoncom.DISPATCH_METHOD | pythoncom.DISPATCH_PROPERTYGET
PROP_FETCHERS = (("FetchCharArrayPrp", "text"),
                 ("FetchSingleCharPrp", "char"),
                 ("FetchSingleNumericPrp", "num"))
RC_OK, RC_PROP_NOT_FOUND, RC_NO_SUCH_LINK = 0, 14, 19
RC_ASSOC_INVALID = (16, 21)


class Late(object):
    def __init__(self, d):
        self.o = d._oleobj_ if hasattr(d, "_oleobj_") else d
        self.ids = {}

    def id(self, name):
        if name not in self.ids:
            r = self.o.GetIDsOfNames(name)
            self.ids[name] = r[0] if isinstance(r, (tuple, list)) else r
        return self.ids[name]

    def has(self, name):
        try:
            self.id(name)
            return True
        except Exception:
            return False

    def call(self, name, *args):
        r = self.o.Invoke(self.id(name), 0, ANY, 1, *args)
        return Late(r) if type(r).__name__ == "PyIDispatch" else r

    def rc(self):
        return self.call("LastReturnCode")


def items(coll):
    """Iterate a MetaModelObjects collection (Count + default member)."""
    def get(i):
        r = coll.o.Invoke(0, 0, ANY, 1, i)
        return Late(r) if type(r).__name__ == "PyIDispatch" else r
    try:
        n = coll.call("Count")
    except Exception:
        n = None
    if n is None:                      # no Count: iterate until it fails
        i = 1
        while True:
            try:
                yield get(i)
            except Exception:
                return
            i += 1
    if not n:
        return
    try:                               # ASSUMPTION: 1-based like Models
        get(n)
        base = 1
    except Exception:
        base = 0
    for i in range(base, base + n):
        yield get(i)


def connect():
    ts = wc.Dispatch("COOLgen.Toolset")
    return Late(ts.Models(1)).call("OpenAPI")


def first_method(api, names):
    for n in names:
        if api.has(n):
            return n
    return None


# --------------------------------------------------------------------------
class Inspector(object):
    def __init__(self, api, a):
        self.api = api
        self.max_prop, self.max_assoc = a.max_prop, a.max_assoc
        self.maxlen = a.max_len
        self.prop_valid = {}      # (otype, fetcher) -> codes
        self.asc_valid = {}       # (otype, fetcher) -> codes
        self.pm, self.om, self.am = {}, {}, {}      # code -> mnemonic
        self.name_code = None
        self.records = {}         # id -> dict for JSON
        self._load_mnemonics(a.mnemonics)
        self._find_name_code()

    # ---- mnemonics ------------------------------------------------------
    def _resolve(self, method, m):
        if not self.api.has(method):
            return None
        try:
            code = self.api.call(method, m)
            return code if self.api.rc() == 0 else None
        except Exception:
            return None

    def _load_mnemonics(self, path):
        if not path:
            return
        with open(path) as f:
            for line in f:
                m = line.strip()
                if not m or m.startswith("#"):
                    continue
                for method, target in (("FetchSchemaPrpTypeByMnem", self.pm),
                                       ("FetchSchemaObjTypeByMnem", self.om),
                                       ("FetchSchemaAscTypeByMnem", self.am)):
                    code = self._resolve(method, m)
                    if code is not None:
                        target[code] = m

    def _find_name_code(self):
        code = self._resolve("FetchSchemaPrpTypeByMnem", "NAME")  # ASSUMPTION
        self.name_code = code

    @staticmethod
    def lab(m, code):
        return "%s(%s)" % (code, m[code]) if code in m else str(code)

    # ---- reading --------------------------------------------------------
    def name(self, o):
        if self.name_code is None:
            return ""
        try:
            s = o.call("FetchCharArrayPrp", self.name_code)
            return s if o.rc() == 0 and isinstance(s, str) else ""
        except Exception:
            return ""

    def props(self, o):
        otype = o.call("ObjectTypeCode")
        out = []
        for fn, kind in PROP_FETCHERS:
            key = (otype, fn)
            known = self.prop_valid.get(key)
            probing = known is None
            codes = range(1, self.max_prop + 1) if probing else sorted(known)
            ok, nf = set(), set()
            for code in codes:
                try:
                    v = o.call(fn, code)
                except Exception:
                    continue
                r = o.rc()
                if probing:
                    if r == RC_OK:
                        ok.add(code)
                    elif r == RC_PROP_NOT_FOUND:
                        nf.add(code)
                if r == RC_OK:
                    out.append((code, kind, v))
            if probing:
                # if "not found" answers for most codes, it does not mean "valid"
                self.prop_valid[key] = ok if len(nf) > self.max_prop // 2 else ok | nf
        return sorted(out, key=lambda x: (x[0], x[1]))

    def assocs(self, o):
        """List of (code, kind, [target Late])"""
        otype = o.call("ObjectTypeCode")
        out = []
        for fn, kind in (("FetchCardManyAsc", "many"), ("FetchCardOneAsc", "one")):
            key = (otype, fn)
            known = self.asc_valid.get(key)
            probing = known is None
            codes = range(1, self.max_assoc + 1) if probing else sorted(known)
            ok, nl = set(), set()
            for code in codes:
                try:
                    res = o.call(fn, code)
                except Exception:
                    continue
                r = o.rc()
                if r in RC_ASSOC_INVALID:
                    continue
                if probing:
                    if r == RC_OK:
                        ok.add(code)
                    elif r == RC_NO_SUCH_LINK:
                        nl.add(code)
                if r != RC_OK or res is None:
                    continue
                targets = list(items(res)) if fn == "FetchCardManyAsc" else [res]
                if targets:
                    out.append((code, kind, targets))
            if probing:
                self.asc_valid[key] = ok if len(nl) > self.max_assoc // 2 else ok | nl
        return out

    def short(self, t):
        return "#%s type=%s '%s'" % (t.call("ID"), self.lab(self.om, t.call("ObjectTypeCode")),
                                     self.name(t))

    def clip(self, v):
        s = repr(v)
        return s if len(s) <= self.maxlen else s[:self.maxlen] + "...'"

    # ---- printing -------------------------------------------------------
    def show(self, o, pad="", with_assocs=True):
        oid = o.call("ID")
        print("%s%s" % (pad, self.short(o)))
        rec = {"id": oid, "type": o.call("ObjectTypeCode"), "name": self.name(o),
               "props": [], "assocs": []}
        for code, kind, v in self.props(o):
            print("%s  prop %-14s [%s] = %s" % (pad, self.lab(self.pm, code), kind, self.clip(v)))
            rec["props"].append({"code": code, "kind": kind, "value": v})
        children = []
        if with_assocs:
            for code, kind, targets in self.assocs(o):
                for t in targets:
                    print("%s  assoc %-13s (%s) -> %s" % (pad, self.lab(self.am, code), kind, self.short(t)))
                    children.append((code, t))
                rec["assocs"].append({"code": code, "kind": kind,
                                      "targets": [t.call("ID") for t in targets]})
        self.records[oid] = rec
        return children


def object_cmd(insp, args):
    o = insp.api.call("GetMetaModelObject", args.id)
    if insp.api.rc() != 0:
        sys.exit("Object %s not found (rc=%s)" % (args.id, insp.api.rc()))
    insp.show(o)


def tree_cmd(insp, args):
    root = insp.api.call("GetMetaModelObject", args.id)
    if insp.api.rc() != 0:
        sys.exit("Object %s not found (rc=%s)" % (args.id, insp.api.rc()))
    seen = {args.id}
    queue = deque([(root, 0)])
    shown = 0
    while queue and shown < args.max_nodes:
        o, depth = queue.popleft()
        children = insp.show(o, pad="    " * depth)
        shown += 1
        if depth < args.depth:
            for _, t in children:
                tid = t.call("ID")
                if tid not in seen:
                    seen.add(tid)
                    queue.append((t, depth + 1))
        print()
    if queue:
        print("-- stopped at --max-nodes %d (%d objects still queued)"
              % (args.max_nodes, len(queue)))


def resolve_type(insp, text):
    if text.isdigit():
        return int(text)
    code = insp._resolve("FetchSchemaObjTypeByMnem", text)
    if code is None:
        sys.exit("Cannot resolve object type %r" % text)
    return code


def scan_cmd(insp, args):
    tcode = resolve_type(insp, args.type)
    if not insp.api.has("FetchAllModelTypeObjs"):
        sys.exit("FetchAllModelTypeObjs is not available on OpenAPI.")
    coll = insp.api.call("FetchAllModelTypeObjs", tcode)
    n = 0
    for o in items(coll):
        insp.show(o)
        print()
        n += 1
        if n >= args.limit:
            break
    print("-- %d objects of type %s shown" % (n, tcode))


def types_cmd(insp, args):
    name = first_method(insp.api, ("FetchAllModelObjects", "FetchAllModelObjs"))
    if name is None:
        sys.exit("No FetchAllModel... method found on OpenAPI.")
    counts, example = Counter(), {}
    for o in items(insp.api.call(name)):
        t = o.call("ObjectTypeCode")
        counts[t] += 1
        example.setdefault(t, (o.call("ID"), insp.name(o)))
    print("%-18s %8s  %s" % ("TYPE", "COUNT", "EXAMPLE"))
    for t, c in sorted(counts.items()):
        eid, en = example[t]
        print("%-18s %8d  #%s '%s'" % (insp.lab(insp.om, t), c, eid, en))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mnemonics")
    ap.add_argument("--json")
    ap.add_argument("--max-prop", type=int, default=800)
    ap.add_argument("--max-assoc", type=int, default=1500)
    ap.add_argument("--max-len", type=int, default=120, help="clip long values")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("types")
    o = sub.add_parser("object")
    o.add_argument("--id", type=int, required=True)
    t = sub.add_parser("tree")
    t.add_argument("--id", type=int, required=True)
    t.add_argument("--depth", type=int, default=2)
    t.add_argument("--max-nodes", type=int, default=60)
    s = sub.add_parser("scan")
    s.add_argument("--type", required=True)
    s.add_argument("--limit", type=int, default=5)
    args = ap.parse_args()

    insp = Inspector(connect(), args)
    {"types": types_cmd, "object": object_cmd, "tree": tree_cmd,
     "scan": scan_cmd}[args.cmd](insp, args)
    if args.json:
        with open(args.json, "w") as f:
            json.dump(list(insp.records.values()), f, indent=1, default=str)
        print("JSON written to", args.json)


if __name__ == "__main__":
    main()
