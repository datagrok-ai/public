"""Quick query client for the long-running KG server.

On the first call it self-installs the query engine (kuzu) into `.kg/.venv`,
restores `kg.kuzu` from the committed `kg.kuzu.xz`, and starts `kg_server.py`
in the background. Every call after that responds in tens of milliseconds.

Usage:
    py qq.py                                                 # show server info
    py qq.py "MATCH (p:Package) RETURN p.name LIMIT 10"      # one-off query
    py qq.py --json "<cypher>"                               # raw JSON output
    py qq.py --stop                                          # shut server down
    py qq.py --restart                                       # stop + start

The server binds to 127.0.0.1 only. PID is recorded in `.kg_server.pid`.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parent   # .kg/scripts/
ROOT = SCRIPTS_DIR.parent                        # .kg/
SERVER = SCRIPTS_DIR / "kg_server.py"
PID_FILE = ROOT / ".kg_server.pid"
LOG_FILE = ROOT / "kg_server.log"
PORT = int(os.environ.get("KG_PORT", "7475"))
URL = f"http://127.0.0.1:{PORT}"

VENV_PY_WIN = ROOT / ".venv" / "Scripts" / "python.exe"
VENV_PY_NIX = ROOT / ".venv" / "bin" / "python"

DB_PATH = ROOT / "kg.kuzu"
UNPACK = SCRIPTS_DIR / "unpack.py"
VENV_DIR = ROOT / ".venv"
REQS = SCRIPTS_DIR / "requirements.txt"


def _python() -> str:
    if VENV_PY_WIN.exists(): return str(VENV_PY_WIN)
    if VENV_PY_NIX.exists(): return str(VENV_PY_NIX)
    return sys.executable


def _venv_python() -> Path:
    return VENV_PY_WIN if sys.platform == "win32" else VENV_PY_NIX


def _ensure_venv() -> None:
    """Ensure .venv exists and has the pinned kuzu. pip + requirements.txt are
    the source of truth: `pip install -r requirements.txt` installs the exact
    pinned version (upgrading if a previous checkout's venv drifted after a pull)
    and is a fast no-op when already satisfied — so no manual version checking is
    needed. qq.py is stdlib-only, so it runs under system python and bootstraps
    the venv that kg_server.py (which imports kuzu) uses."""
    vpy = _venv_python()
    if not vpy.exists():
        print("[qq] first run: setting up the query engine (one-time, ~30s)…", file=sys.stderr)
        if subprocess.run([sys.executable, "-m", "venv", str(VENV_DIR)]).returncode != 0:
            sys.exit("[qq] could not create .venv. Install manually: "
                     "python -m venv .venv && .venv/bin/pip install -r .kg/scripts/requirements.txt")
    pip = [str(vpy), "-m", "pip", "install", "-q", "--disable-pip-version-check"]
    pip += ["-r", str(REQS)] if REQS.exists() else ["kuzu"]
    if subprocess.run(pip).returncode != 0:
        sys.exit(f"[qq] failed to install kuzu. Run: {vpy} -m pip install -r {REQS}")


def is_up(timeout_s: float = 0.3) -> bool:
    try:
        with urllib.request.urlopen(f"{URL}/ping", timeout=timeout_s) as r:
            return r.status == 200
    except (urllib.error.URLError, ConnectionRefusedError, OSError):
        return False


def start_background() -> None:
    """Spawn kg_server.py detached so it survives the parent exiting."""
    args = [_python(), str(SERVER)]
    log = LOG_FILE.open("ab", buffering=0)
    if sys.platform == "win32":
        # CREATE_NO_WINDOW | DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP
        flags = 0x00000008 | 0x00000200 | 0x08000000
        subprocess.Popen(
            args, stdout=log, stderr=log,
            cwd=str(ROOT), creationflags=flags, close_fds=True,
        )
    else:
        subprocess.Popen(
            args, stdout=log, stderr=log,
            cwd=str(ROOT), start_new_session=True, close_fds=True,
        )


def _ensure_db() -> None:
    """If the materialized DB is missing (fresh clone), restore it from the
    committed kg.kuzu.xz artifact before the server opens it."""
    if DB_PATH.exists(): return
    if not UNPACK.exists():
        sys.exit(f"[qq] {DB_PATH.name} missing and no scripts/unpack.py — the .kg bundle is incomplete")
    print(f"[qq] {DB_PATH.name} missing — restoring from kg.kuzu.xz …", file=sys.stderr)
    if subprocess.run([_python(), str(UNPACK)]).returncode != 0:
        sys.exit("[qq] could not restore kg.kuzu — kg.kuzu.xz is missing from the bundle.")


def ensure_up(max_wait_s: float = 20.0) -> None:
    if is_up(): return
    _ensure_venv()
    _ensure_db()
    print(f"[qq] starting kg_server in background (logs: {LOG_FILE})…", file=sys.stderr)
    start_background()
    deadline = time.time() + max_wait_s
    while time.time() < deadline:
        time.sleep(0.25)
        if is_up():
            print(f"[qq] server ready in {round(max_wait_s - (deadline - time.time()), 1)}s",
                  file=sys.stderr)
            return
    sys.exit(f"[qq] server did not respond within {max_wait_s}s. See {LOG_FILE}")


def stop() -> None:
    if not is_up():
        print("[qq] server not running", file=sys.stderr); return
    try:
        urllib.request.urlopen(f"{URL}/stop", timeout=2).read()
    except Exception:
        pass
    # Wait a beat, fall back to killing by PID if necessary
    for _ in range(20):
        time.sleep(0.1)
        if not is_up(): break
    if PID_FILE.exists():
        try:
            pid = int(PID_FILE.read_text())
            if sys.platform == "win32":
                subprocess.run(["taskkill", "/F", "/PID", str(pid)],
                               capture_output=True)
            else:
                os.kill(pid, 9)
        except Exception:
            pass
        try: PID_FILE.unlink()
        except OSError: pass
    print("[qq] server stopped", file=sys.stderr)


def query(cypher: str, limit: int = 1000) -> dict:
    body = json.dumps({"cypher": cypher, "limit": limit}).encode("utf-8")
    req = urllib.request.Request(
        f"{URL}/query", data=body,
        headers={"Content-Type": "application/json"},
    )
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            return json.loads(r.read())
    except urllib.error.HTTPError as e:
        return json.loads(e.read())


def info() -> dict:
    with urllib.request.urlopen(f"{URL}/info", timeout=5) as r:
        return json.loads(r.read())


def render_table(cols: list[str], rows: list[list]) -> str:
    if not rows: return "(no rows)"
    def fmt(v):
        if v is None: return ""
        if isinstance(v, (list, dict)):
            return json.dumps(v, ensure_ascii=False)
        return str(v)
    header = cols[:]
    string_rows = [[fmt(v) for v in r] for r in rows]
    widths = [max(len(header[i]), max((len(r[i]) for r in string_rows), default=0))
              for i in range(len(header))]
    sep = "  "
    out = [sep.join(h.ljust(w) for h, w in zip(header, widths))]
    out.append(sep.join("-" * w for w in widths))
    for r in string_rows:
        out.append(sep.join(c.ljust(w) for c, w in zip(r, widths)))
    return "\n".join(out)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("cypher", nargs="?")
    ap.add_argument("--limit", type=int, default=1000)
    ap.add_argument("--json", action="store_true",
                    help="Emit raw JSON instead of a tabular view")
    ap.add_argument("--stop", action="store_true")
    ap.add_argument("--restart", action="store_true")
    ap.add_argument("--no-start", action="store_true",
                    help="Fail if server isn't already running (don't auto-spawn)")
    args = ap.parse_args()

    if args.stop:
        stop(); return 0
    if args.restart:
        stop(); time.sleep(0.5); ensure_up(); print(json.dumps(info())); return 0

    if args.no_start:
        if not is_up():
            sys.exit("[qq] server not running and --no-start passed")
    else:
        ensure_up()

    if not args.cypher:
        print(json.dumps(info(), indent=2))
        return 0

    res = query(args.cypher, limit=args.limit)
    if args.json:
        print(json.dumps(res, indent=2, ensure_ascii=False))
        return 1 if "error" in res else 0
    if "error" in res:
        print(f"ERROR ({res.get('type', '')}): {res['error']}", file=sys.stderr)
        return 1
    print(render_table(res["columns"], res["rows"]))
    foot = f"\n{res['row_count']} row(s)  ·  {res['elapsed_ms']}ms"
    if res.get("truncated"):
        foot += "  ·  truncated"
    print(foot, file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
