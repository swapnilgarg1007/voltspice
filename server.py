#!/usr/bin/env python3
"""Tiny dev server: serves web/ and runs the voltspice CLI for /api/simulate. Standard library only."""
import json, os, subprocess, sys
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

ROOT = os.path.dirname(os.path.abspath(__file__))
WEB = os.path.join(ROOT, "web")
EXAMPLES = os.path.join(ROOT, "examples")
BIN = os.environ.get("VOLTSPICE_BIN", os.path.join(ROOT, "build", "voltspice"))
MAX_BODY = 200_000
TYPES = {".html": "text/html", ".js": "text/javascript", ".css": "text/css"}


class Handler(BaseHTTPRequestHandler):
    def _send(self, code, body, ctype="application/json"):
        data = body if isinstance(body, bytes) else body.encode()
        self.send_response(code)
        self.send_header("Content-Type", ctype + "; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        if self.path == "/api/examples":
            out = {}
            for name in sorted(os.listdir(EXAMPLES)):
                if name.endswith(".cir"):
                    with open(os.path.join(EXAMPLES, name)) as f:
                        out[name[:-4]] = f.read()
            return self._send(200, json.dumps(out))
        path = "/index.html" if self.path in ("/", "") else self.path.split("?")[0]
        full = os.path.realpath(os.path.join(WEB, path.lstrip("/")))
        if not full.startswith(WEB + os.sep) or not os.path.isfile(full):
            return self._send(404, "not found", "text/plain")
        with open(full, "rb") as f:
            self._send(200, f.read(), TYPES.get(os.path.splitext(full)[1], "application/octet-stream"))

    def do_POST(self):
        if self.path != "/api/simulate":
            return self._send(404, "not found", "text/plain")
        n = int(self.headers.get("Content-Length", 0))
        if n > MAX_BODY:
            return self._send(413, json.dumps({"error": "netlist too large"}))
        netlist = self.rfile.read(n).decode("utf-8", "replace")
        try:
            p = subprocess.run([BIN, "--json", "-"], input=netlist, capture_output=True, text=True, timeout=15)
            out = p.stdout.strip() or json.dumps({"error": p.stderr.strip() or "simulator produced no output"})
        except FileNotFoundError:
            out = json.dumps({"error": "simulator binary not found; build the project first"})
        except subprocess.TimeoutExpired:
            out = json.dumps({"error": "simulation timed out"})
        self._send(200, out)

    def log_message(self, *a):
        pass


if __name__ == "__main__":
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8000
    print(f"voltspice UI on http://localhost:{port}  (engine: {BIN})")
    ThreadingHTTPServer(("127.0.0.1", port), Handler).serve_forever()
