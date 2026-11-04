import win32com.client as wc
ts = wc.GetActiveObject("COOLgen.Toolset")   # or Dispatch, same result for you
models = ts.Models
print("Count:", models.Count)

attempts = {
    "Item(1)":      lambda: models.Item(1),
    "Item(0)":      lambda: models.Item(0),
    "models(1)":    lambda: models(1),
    "iterate":      lambda: next(iter(models)),
}
for label, fn in attempts.items():
    try:
        print(f"{label:10} ->", fn())
    except Exception as e:
        print(f"{label:10} failed:", e)

# What the Models object actually exposes (invkind 1 = method, 2 = property get)
ti = models._oleobj_.GetTypeInfo()
for i in range(ti.GetTypeAttr().cFuncs):
    fd = ti.GetFuncDesc(i)
    print(ti.GetNames(fd.memid), "invkind", fd.invkind)