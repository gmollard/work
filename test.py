import win32com.client as wc
ts = wc.GetActiveObject("COOLgen.Toolset")
models = ts.Models
print("Count:", models.Count)

for label, fn in {
    "models(1)": lambda: models(1),
    "models(0)": lambda: models(0),
    "iterate":   lambda: list(models)[0],
}.items():
    try:
        print(label, "->", fn())
    except Exception as e:
        print(label, "failed:", e)

Set ts = CreateObject("COOLgen.Toolset")
WScript.Echo "Count: " & ts.Models.Count
Set m = ts.Models.Item(1)
If m Is Nothing Then
  WScript.Echo "Item(1) is Nothing"
Else
  WScript.Echo "Model OK: " & TypeName(m)
  WScript.Echo "Objects: " & m.OpenAPI.CountModelObjs
End If

import comtypes.client as cc
ts = cc.CreateObject("COOLgen.Toolset", dynamic=True)
print(ts.Models.Count, ts.Models.Item(1))


import pythoncom
import win32com.client as wc

ts = wc.GetActiveObject("COOLgen.Toolset")
model = ts.Models(1)
ole = model._oleobj_

# 1. Normal late binding
try:
    api = model.OpenAPI
    print("model.OpenAPI ->", api)
except Exception as e:
    print("model.OpenAPI failed:", e)

# 2. Can this object resolve names at all?
for name in ("OpenAPI", "Name", "Item"):
    try:
        print(name, "->", ole.GetIDsOfNames(name))
    except Exception as e:
        print(name, "-> failed:", e)

# 3. Probe member ids directly, property reads only
for dispid in range(0, 60):
    try:
        r = ole.Invoke(dispid, 0, pythoncom.DISPATCH_PROPERTYGET, 1)
        print("dispid", dispid, "->", type(r).__name__, r)
    except pythoncom.com_error:
        pass