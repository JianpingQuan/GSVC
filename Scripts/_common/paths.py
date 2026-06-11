import os
from pathlib import Path

def gsvc_root() -> Path:
    env = os.environ.get("GSVC_ROOT") or os.environ.get("PHAGE_PROJECT_ROOT")
    if env and Path(env).is_dir():
        return Path(env).resolve()
    here = Path(__file__).resolve().parents[2]
    return here

PHAGE_ROOT = gsvc_root()
OUT_DIR = PHAGE_ROOT / "Scripts" / "output"
OUT_DIR.mkdir(parents=True, exist_ok=True)
os.chdir(PHAGE_ROOT)
