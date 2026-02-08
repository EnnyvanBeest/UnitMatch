# tools/bump_version.py
from pathlib import Path
import tomllib, re, sys

pyproj = Path("UnitMatchPy/pyproject.toml")
data = tomllib.loads(pyproj.read_text(encoding="utf-8"))
ver = data["project"]["version"]
maj, minr, patch = map(int, ver.split("."))
if len(sys.argv) > 1 and sys.argv[1] in {"major","minor","patch"}:
    kind = sys.argv[1]
else:
    kind = "patch"

if kind == "major":
    maj, minr, patch = maj+1, 0, 0
elif kind == "minor":
    minr, patch = minr+1, 0
else:
    patch += 1

newver = f"{maj}.{minr}.{patch}"
print(f"{ver} -> {newver}")

txt = pyproj.read_text(encoding="utf-8")
txt = re.sub(r'(?m)^(version\s*=\s*")[0-9]+\.[0-9]+\.[0-9]+(")', rf'\g<1>{newver}\2', txt)
pyproj.write_text(txt, encoding="utf-8")
