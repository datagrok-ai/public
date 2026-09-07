/* Static file server for the gallery-hosted checks (tour, message-input): those behaviours are
   platform-free, so they run against the plain gallery instead of the app — and the gallery is
   ES modules over the built `src/`, which needs a real origin, not file://. */
import {createReadStream} from 'node:fs';
import {stat} from 'node:fs/promises';
import {createServer} from 'node:http';
import {extname, join, normalize, resolve} from 'node:path';
import {fileURLToPath} from 'node:url';

const ROOT = resolve(fileURLToPath(new URL('..', import.meta.url)));
const SERVED = ['gallery', 'css', 'src', 'vendor', 'node_modules/@floating-ui',
  'node_modules/datagrok-api/src/u2core'];
const MIME = {'.html': 'text/html; charset=utf-8', '.js': 'text/javascript; charset=utf-8',
  '.mjs': 'text/javascript; charset=utf-8', '.css': 'text/css; charset=utf-8',
  '.json': 'application/json', '.map': 'application/json',
  '.woff2': 'font/woff2', '.svg': 'image/svg+xml'};

/** Listens on a free loopback port; `close()` resolves once the port is actually released. */
export function startGalleryServer() {
  const server = createServer(async (req, res) => {
    try {
      let path = decodeURIComponent(req.url.split(/[?#]/)[0]).replace(/^\/+/, '');
      if (path === '' || path.endsWith('/'))
        path += 'index.html';
      const rel = normalize(path).split('\\').join('/');
      const file = join(ROOT, rel);
      if (!file.startsWith(ROOT) || !SERVED.some((dir) => rel === dir || rel.startsWith(`${dir}/`)))
        return void res.writeHead(403).end('forbidden');
      try {
        await stat(file);
      } catch (e) {
        return void res.writeHead(404).end(`not found: ${rel}`);
      }
      res.writeHead(200, {'content-type': MIME[extname(file)] ?? 'application/octet-stream'});
      createReadStream(file).on('error', () => res.destroy()).pipe(res);
    } catch (e) {
      res.writeHead(400).end('bad request');
    }
  });
  return new Promise((listening) => server.listen(0, '127.0.0.1', () => listening({
    url: `http://127.0.0.1:${server.address().port}`,
    close: () => new Promise((closed) => server.close(closed)),
  })));
}
