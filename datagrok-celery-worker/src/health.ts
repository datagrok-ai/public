/** Tiny health endpoint (the worker container EXPOSEs 8000). */
import * as http from 'node:http';

import {logInfo} from './logger';

export function startHealthServer(port: number): http.Server {
  const server = http.createServer((_req, res) => {
    res.writeHead(200, {'Content-Type': 'application/json'});
    res.end(JSON.stringify({'status': 'ok'}));
  });
  server.listen(port, () => logInfo(`Health endpoint listening on port ${port}`));
  return server;
}
