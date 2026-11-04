#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
gen_api_test.py - step-by-step test harness for the CA Gen Toolset OLE
automation API (OpenAPI / MetaModelObjects / MetaModelObject objects).

Run it on the Windows machine where the Gen Toolset is running with a
TEST model open.  Requires:  pip install pywin32

Commands (everything is read-only except `swap --apply`):

  connect                     Attach to the Toolset, count the objects.
  dump [--filter TEXT]        List objects: type code, name, object id.
  assocs --id N               List every association leaving object N.
  whereused --id N            Find every object holding an association
                              that points to object N (e.g. the USE
                              statement that points to the old block).
  swap --holder H --code C --old O --new N [--apply]
                              Repoint association C of object H from O
                              to N.  Dry run unless --apply is given.

Sources for the calls used (Gen 8.5 Encyclopedia API Reference Guide,
appendix "CA Gen Toolset Automation"):
  CreateObject("COOLgen.Toolset") -> .Models.Item(1) -> .OpenAPI
  OpenAPI: CountModelObjs, FetchAllModelObjs, GetMetaModelObject,
           FetchSchemaPrpTypeByMnem, Rollback
  MetaModelObject: ID, ObjectTypeCode, LastReturnCode, FetchCharArrayPrp,
           FetchCardManyAsc, FetchCardOneAsc, AddAssoc, DelAssoc
Anything marked ASSUMPTION below is not confirmed by the documentation.
"""
import argparse
import sys

try:
    import win32com.client as wc
    import pywintypes
except ImportError:
    sys.exit("pywin32 is required:  pip install pywin32")

PROGID = "COOLgen.Toolset"

# Return codes (documented: "Return Codes" page of the guide)
RC_NAMES = {
    0: "SUCCESSFUL", 11: "MODEL_ALREADY_LOCKED", 13: "OBJECT_NOT_FOUND",
    14: "PROP_NOT_FOUND", 15: "PTC_INVALID", 16: "ATC_INVALID",
    17: "OTC_INVALID", 18: "MNEMONIC_INVALID", 19: "DEST_OBJ_NOT_FOUND",
    21: "ASSOC_NOT_FOUND", 32: "MODEL_LOCKED", 38: "MODEL_NOT_LOCKED",
    40: "OBJECT_PROTECTED", 42: "FROM_OBJECT_NOT_FOUND",
    43: "TO_OBJECT_NOT_FOUND", 44: "INVALID_CARDINALITY",
    45: "ASSOC_PROTECTED", 53: "VALUE_INVALID", 54: "NOT_AUTHORIZED",
    55: "ASSOC_ALREADY_EXISTS", 56: "COMMIT_UNSUCCESSFUL",
}
SUCCESS = 0
INVALID_FOR_OBJECT = (16, 21)   # ATC_INVALID, ASSOC_NOT_FOUND
NO_SUCH_LINK = 19               # DEST_OBJ_NOT_FOUND (code valid, link absent)


def rc_name(r):
    return "%s (%s)" % (r, RC_NAMES.get(r, "see Return Codes page"))


# --------------------------------------------------------------------------
# COM helpers.  Depending on the type library, pywin32 exposes a no-argument
# member either as a property (value returned) or as a method (callable
# returned); these helpers accept both.
# --------------------------------------------------------------------------
def scalar(v):
    if v is None or isinstance(v, (int, float, str, bool)):
        return v
    return v()


def com(v):
    return v if hasattr(v, "_oleobj_") else v()


def rc(o):
    return scalar(o.LastReturnCode)


def items(coll):
    """Iterate a MetaModelObjects collection (count + Item)."""
    n = scalar(coll.Count)
    if not n:
        return
    # ASSUMPTION: base is not documented for MetaModelObjects.Item.
    # Models.Item(1) is documented as 1-based, so try that first.
    try:
        coll.Item(n)
        base = 1
    except Exception:
        base = 0
    for i in range(base, base + n):
        yield coll.Item(i)


# --------------------------------------------------------------------------
# Connection
# --------------------------------------------------------------------------
def connect():
    try:
        toolset = wc.GetActiveObject(PROGID)
        how = "running instance (GetActiveObject)"
    except Exception:
        toolset = wc.Dispatch(PROGID)
        how = "Dispatch (CreateObject equivalent)"
    models = com(toolset.Models)
    try:
        n_models = scalar(models.Count)
    except Exception:
        n_models = "?"
    print("Attached via %s; models open: %s" % (how, n_models))
    if n_models == 0:
        raise RuntimeError("No model is open in the Toolset. Open the test "
                           "model first.")
    model = com(models.Item(1))          # documented: the item will be 1
    api = com(model.OpenAPI)             # documented: OpenAPI property
    return api


def name_prop_code(api, mnemonic):
    # ASSUMPTION: the property mnemonic of an object's name is "NAME".
    code = scalar(api.FetchSchemaPrpTypeByMnem(mnemonic))
    if rc(api) != SUCCESS:
        print("WARNING: mnemonic %r not accepted (%s). Names will be blank; "
              "look the mnemonic up in odrpta.chm and pass --name-mnemonic."
              % (mnemonic, rc_name(rc(api))))
        return None
    return code


def obj_name(o, ncode):
    if ncode is None:
        return ""
    try:
        s = o.FetchCharArrayPrp(ncode)
    except Exception:
        return ""
    return s if rc(o) == SUCCESS and isinstance(s, str) else ""


def get_obj(api, oid):
    o = com(api.GetMetaModelObject(oid))
    if rc(api) != SUCCESS:
        raise RuntimeError("GetMetaModelObject(%s) -> %s" % (oid, rc_name(rc(api))))
    return o


def describe(o, ncode):
    return "#%s type=%s '%s'" % (scalar(o.ID), scalar(o.ObjectTypeCode),
                                 obj_name(o, ncode))


# --------------------------------------------------------------------------
# Association walker.  The API has no "list all associations" call, so we
# probe association type codes 1..max_code with FetchCardManyAsc /
# FetchCardOneAsc.  Codes that are invalid for an object type are cached so
# the expensive probing happens once per type.
# --------------------------------------------------------------------------
class Walker(object):
    def __init__(self, max_code):
        self.max_code = max_code
        self.valid = {}      # (object type, fetch function) -> set(codes)

    def edges(self, o):
        """Yield (assoc_code, fetch_function, target_object)."""
        otype = scalar(o.ObjectTypeCode)
        for fetch in ("FetchCardManyAsc", "FetchCardOneAsc"):
            key = (otype, fetch)
            known = self.valid.get(key)
            probing = known is None
            codes = range(1, self.max_code + 1) if probing else sorted(known)
            found = set()
            for code in codes:
                try:
                    res = getattr(o, fetch)(code)
                except Exception:
                    continue
                r = rc(o)
                if r in INVALID_FOR_OBJECT:
                    continue
                if probing and r in (SUCCESS, NO_SUCH_LINK):
                    found.add(code)
                if r != SUCCESS or res is None:
                    continue
                if fetch == "FetchCardManyAsc":
                    for t in items(res):
                        yield code, fetch, t
                else:
                    yield code, fetch, res
            if probing:
                self.valid[key] = found


# --------------------------------------------------------------------------
# Commands
# --------------------------------------------------------------------------
def cmd_connect(api, args):
    n = scalar(api.CountModelObjs)
    print("CountModelObjs = %s, rc = %s" % (n, rc_name(rc(api))))
    print("OK - the Toolset automation interface is reachable.")


def cmd_dump(api, args):
    ncode = name_prop_code(api, args.name_mnemonic)
    coll = com(api.FetchAllModelObjs)
    rows = []
    for o in items(coll):
        rows.append((scalar(o.ObjectTypeCode), obj_name(o, ncode), scalar(o.ID)))
    rows.sort(key=lambda r: (r[0], r[1].lower(), r[2]))
    flt = args.filter.lower() if args.filter else None
    shown = 0
    print("%-8s %-10s %s" % ("TYPE", "ID", "NAME"))
    for otype, name, oid in rows:
        if flt and flt not in name.lower():
            continue
        print("%-8s %-10s %s" % (otype, oid, name))
        shown += 1
    print("-- %d shown, %d objects in model" % (shown, len(rows)))


def cmd_assocs(api, args):
    ncode = name_prop_code(api, args.name_mnemonic)
    o = get_obj(api, args.id)
    print("Object", describe(o, ncode))
    seen = set()
    for code, fetch, t in Walker(args.max_code).edges(o):
        key = (code, scalar(t.ID))
        if key in seen:
            continue
        seen.add(key)
        print("  assoc %-5s via %-17s -> %s" % (code, fetch, describe(t, ncode)))
    if not seen:
        print("  (no associations found; try a larger --max-code)")


def cmd_whereused(api, args):
    ncode = name_prop_code(api, args.name_mnemonic)
    target = get_obj(api, args.id)
    print("Looking for objects with an association to", describe(target, ncode))
    walker = Walker(args.max_code)
    objs = list(items(com(api.FetchAllModelObjs)))
    hits = set()
    for i, o in enumerate(objs, 1):
        if i % 25 == 0:
            print("  ... %d / %d objects scanned" % (i, len(objs)))
        for code, fetch, t in walker.edges(o):
            if scalar(t.ID) == args.id:
                hits.add((scalar(o.ID), code, fetch))
    if not hits:
        print("No holder found. Either nothing references this object, or "
              "the link is not exposed as an association (STOP and read the "
              "meta-model docs in odrpta.chm).")
        return
    for hid, code, fetch in sorted(hits):
        h = get_obj(api, hid)
        print("HOLDER %s  assoc code %s (found via %s)"
              % (describe(h, ncode), code, fetch))


def cmd_swap(api, args):
    ncode = name_prop_code(api, args.name_mnemonic)
    holder = get_obj(api, args.holder)
    old = get_obj(api, args.old)
    new = get_obj(api, args.new)
    print("holder:", describe(holder, ncode))
    print("old   :", describe(old, ncode))
    print("new   :", describe(new, ncode))

    if scalar(old.ObjectTypeCode) != scalar(new.ObjectTypeCode) and not args.force:
        sys.exit("ABORT: old and new objects have different type codes "
                 "(use --force to override).")

    # Pre-check: does the holder really link to old with this code?
    linked = False
    for code, fetch, t in Walker(args.max_code).edges(holder):
        if code == args.code and scalar(t.ID) == args.old:
            linked = True
            break
    print("holder --%s--> old exists: %s" % (args.code, linked))
    if not linked:
        sys.exit("ABORT: that association does not exist (check ids/code).")

    if not args.apply:
        print("\nDRY RUN - nothing changed. Re-run with --apply on a COPY "
              "of the test model.")
        return

    print("\nAPPLYING: DelAssoc, then AddAssoc ...")
    # nTriggerFlag = False -> delete ONLY the association, never related objects
    holder.DelAssoc(args.code, args.old, False)
    r1 = rc(holder)
    print("DelAssoc ->", rc_name(r1))
    if r1 != SUCCESS:
        sys.exit("ABORT: delete failed, model unchanged.")

    holder.AddAssoc(args.code, args.new)
    r2 = rc(holder)
    print("AddAssoc ->", rc_name(r2))
    if r2 != SUCCESS:
        print("AddAssoc failed after the delete succeeded: rolling back.")
        try:
            api.Rollback()          # documented title; syntax shows UpRollback
        except Exception:
            api.UpRollback()
        sys.exit("ROLLED BACK. Per the docs a rollback does not fully undo "
                 "changes: do NOT save this model; close it without saving "
                 "and reopen your backup copy.")

    print("\nDone. Verify:")
    print("  1. python gen_api_test.py assocs --id %s" % args.holder)
    print("  2. Close and reopen the caller's action diagram in the Toolset")
    print("     (documented: a plug-in's changes only show after re-entry).")
    print("  3. Check the USE statement and its view mappings, run the "
          "Toolset verification, then File > Save.")


def build_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--name-mnemonic", default="NAME",
                   help="property mnemonic holding an object's name (ASSUMPTION: NAME)")
    p.add_argument("--max-code", type=int, default=1500,
                   help="highest association type code to probe (default 1500)")
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("connect")
    d = sub.add_parser("dump")
    d.add_argument("--filter", help="only names containing this text")
    a = sub.add_parser("assocs")
    a.add_argument("--id", type=int, required=True)
    w = sub.add_parser("whereused")
    w.add_argument("--id", type=int, required=True)
    s = sub.add_parser("swap")
    s.add_argument("--holder", type=int, required=True, help="id of the USE statement object")
    s.add_argument("--code", type=int, required=True, help="association type code")
    s.add_argument("--old", type=int, required=True)
    s.add_argument("--new", type=int, required=True)
    s.add_argument("--apply", action="store_true")
    s.add_argument("--force", action="store_true")
    return p


COMMANDS = {"connect": cmd_connect, "dump": cmd_dump, "assocs": cmd_assocs,
            "whereused": cmd_whereused, "swap": cmd_swap}


def main():
    args = build_parser().parse_args()
    try:
        api = connect()
        COMMANDS[args.cmd](api, args)
    except pywintypes.com_error as e:
        print("COM error:", e)
        print("Hints: 'Class not registered' -> try 32-bit Python (the Toolset "
              "is registered under Wow6432Node); no model -> open it in the "
              "Toolset first; still failing -> register the script as a "
              "Toolset plug-in (see the guide).")
        sys.exit(1)
    except RuntimeError as e:
        print("ERROR:", e)
        sys.exit(1)


if __name__ == "__main__":
    main()
