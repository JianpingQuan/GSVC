"""Read/write text with automatic Chinese Windows encoding detection."""
from __future__ import annotations

from pathlib import Path


def read_text_auto(path: Path) -> str:
    """Decode script text from UTF-8, UTF-8-BOM, or GBK/GB18030 (common on Windows R)."""
    raw = path.read_bytes()
    if raw.startswith(b"\xef\xbb\xbf"):
        return raw.decode("utf-8-sig")
    for enc in ("utf-8", "gb18030", "gbk", "cp936"):
        try:
            text = raw.decode(enc)
        except UnicodeDecodeError:
            continue
        # Reject UTF-8 that only decodes with replacement characters
        if enc == "utf-8" and text.count("\ufffd") > 2:
            continue
        return text
    return raw.decode("utf-8", errors="replace")


def write_text_utf8(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")
