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
