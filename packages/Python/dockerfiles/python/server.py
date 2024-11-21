from http.server import BaseHTTPRequestHandler, HTTPServer
import json

class CodeExecutionHandler(BaseHTTPRequestHandler):
    def do_POST(self):
        content_length = int(self.headers['Content-Length'])
        post_data = self.rfile.read(content_length)
        response = {"status": "ok", "output": None, "error": None}

        try:
            code = post_data.decode("utf-8")
            exec_globals = {}
            exec_locals = {}
            exec(code, exec_globals, exec_locals)
            response["output"] = exec_locals.get("result", "No 'result' variable returned")
        except Exception as e:
            response["status"] = "error"
            response["error"] = str(e)

        self.send_response(200)
        self.send_header('Content-type', 'application/json')
        self.end_headers()
        self.wfile.write(json.dumps(response).encode("utf-8"))

def run(server_class=HTTPServer, handler_class=CodeExecutionHandler, port=8080):
    server_address = ('', port)
    httpd = server_class(server_address, handler_class)
    print(f"Starting server on port {port}")
    httpd.serve_forever()

if __name__ == "__main__":
    run()