"""Long-running HTTP query server for the Kuzu knowledge graph.

Avoids the ~1.5s Python+Kuzu cold-start overhead per query. Keeps the DB
loaded in memory and serves Cypher queries over HTTP at 127.0.0.1:7475.

Endpoints:
    GET  /ping           -> {"ok": true}
    GET  /info           -> {"nodes": N, "edges": M, "uptime_s": …}
    GET  /stop           -> shuts the server down
    POST /query          -> {"cypher": "...", "limit": 1000}
                           returns {"columns": [...], "rows": [[...], ...], "row_count": N}

Start manually with `py kg_server.py`, or via `qq.py` which auto-starts it
in the background if it isn't already running.

Bind address is loopback only — never exposed to the network.
"""

from __future__ import annotations

import json
import os
import sys
import threading
import time
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path
from socketserver import ThreadingMixIn

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import kuzu

DB_PATH = ROOT / "kg.kuzu"
PID_FILE = ROOT / ".kg_server.pid"
PORT = int(os.environ.get("KG_PORT", "7475"))

# Module-level state kept alive across requests
_db: kuzu.Database | None = None
_lock = threading.Lock()              # Kuzu connections aren't fully thread-safe
_started_at = time.time()
_request_count = 0


def _open_db() -> None:
    global _db
    if not DB_PATH.exists():
        sys.exit(f"DB not found at {DB_PATH}. Run `py build.py` first.")
    _db = kuzu.Database(str(DB_PATH))
    print(f"[kg_server] DB opened: {DB_PATH}", flush=True)


def _coerce(v):
    if v is None or isinstance(v, (str, int, float, bool)):
        return v
    if isinstance(v, (list, tuple)):
        return [_coerce(x) for x in v]
    if isinstance(v, dict):
        return {k: _coerce(x) for k, x in v.items()}
    return str(v)


class ThreadedHTTPServer(ThreadingMixIn, HTTPServer):
    daemon_threads = True
    allow_reuse_address = True


class Handler(BaseHTTPRequestHandler):

    def log_message(self, *_args):       # suppress per-request stdout noise
        return

    # -------------------------------------------------------------------- write
    def _json(self, code: int, obj: dict) -> None:
        body = json.dumps(obj, default=str, ensure_ascii=False).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        try:
            self.wfile.write(body)
        except (BrokenPipeError, ConnectionResetError):
            pass

    # ---------------------------------------------------------------- routes
    def do_GET(self) -> None:
        global _request_count
        _request_count += 1
        if self.path == "/ping":
            return self._json(200, {"ok": True})

        if self.path == "/info":
            with _lock:
                conn = kuzu.Connection(_db)
                n_nodes = self._scalar(conn, "MATCH (n) RETURN count(n);")
                n_edges = self._scalar(conn, "MATCH ()-[r]->() RETURN count(r);")
            return self._json(200, {
                "nodes": n_nodes, "edges": n_edges,
                "uptime_s": round(time.time() - _started_at, 1),
                "requests_handled": _request_count,
                "db_path": str(DB_PATH),
            })

        if self.path == "/stop":
            self._json(200, {"stopping": True})
            threading.Thread(target=self.server.shutdown, daemon=True).start()
            return

        return self._json(404, {"error": "not found", "tried": self.path})

    def do_POST(self) -> None:
        global _request_count
        _request_count += 1
        if self.path != "/query":
            return self._json(404, {"error": "not found"})

        length = int(self.headers.get("Content-Length", 0))
        try:
            payload = json.loads(self.rfile.read(length))
        except json.JSONDecodeError as e:
            return self._json(400, {"error": f"bad json: {e}"})

        cypher = payload.get("cypher")
        if not cypher or not isinstance(cypher, str):
            return self._json(400, {"error": "missing 'cypher' string"})
        limit = int(payload.get("limit", 1000))

        try:
            with _lock:
                conn = kuzu.Connection(_db)
                t0 = time.time()
                res = conn.execute(cypher)
                cols = res.get_column_names()
                rows: list[list] = []
                while res.has_next() and len(rows) < limit:
                    rows.append([_coerce(v) for v in res.get_next()])
                elapsed_ms = round((time.time() - t0) * 1000, 1)
        except Exception as e:
            return self._json(400, {"error": str(e), "type": type(e).__name__})

        return self._json(200, {
            "columns": cols, "rows": rows, "row_count": len(rows),
            "elapsed_ms": elapsed_ms,
            "truncated": res.has_next(),
        })

    @staticmethod
    def _scalar(conn, cypher: str):
        res = conn.execute(cypher)
        if res.has_next():
            return _coerce(res.get_next()[0])
        return None


def main() -> int:
    _open_db()
    PID_FILE.write_text(str(os.getpid()), encoding="utf-8")
    try:
        srv = ThreadedHTTPServer(("127.0.0.1", PORT), Handler)
        print(f"[kg_server] listening on http://127.0.0.1:{PORT}", flush=True)
        srv.serve_forever()
    finally:
        try:
            PID_FILE.unlink()
        except OSError:
            pass
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
