import pythoncom
import win32com.client as wc

ts = wc.GetActiveObject("COOLgen.Toolset")
models = ts.Models
ole = models._oleobj_

print("typeinfo count:", ole.GetTypeInfoCount())
print("Count:", models.Count)

dispid = ole.GetIDsOfNames("Item")[0]
modes = {
    "METHOD":  pythoncom.DISPATCH_METHOD,
    "PROPGET": pythoncom.DISPATCH_PROPERTYGET,
    "BOTH":    pythoncom.DISPATCH_METHOD | pythoncom.DISPATCH_PROPERTYGET,
}
for label, flags in modes.items():
    for idx in (1, 0):
        try:
            r = ole.Invoke(dispid, 0, flags, 1, idx)
            print(f"{label:8} Item({idx}) ->", repr(r))
        except Exception as e:
            print(f"{label:8} Item({idx}) failed:", e)