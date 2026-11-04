import psutil, win32com.client as wc

def toolset_pids():
    return {p.pid: p.info["exe"] for p in psutil.process_iter(["name", "exe"])
            if (p.info["name"] or "").lower() == "toolset.exe"}

before = toolset_pids()
print("Before:", before)

ts = wc.Dispatch("COOLgen.Toolset")

after = toolset_pids()
print("After :", after)
print("NEW processes:", set(after) - set(before))
print("Models.Count:", ts.Models.Count)