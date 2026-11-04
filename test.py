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