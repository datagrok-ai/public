#!/usr/bin/env node
// Dev-loop control for a LOCAL mcp-server container.
//
// The platform's own mcp-server containers run whatever image the last `grok publish` deployed,
// so iterating on the op registry against them means an image rebuild + redeploy per change. This
// one is built from the working tree and reached directly, the same way dev/runtime.mjs works for
// claude-runtime. Point a turn at it with mcpServerUrl = http://localhost:<port>/mcp — the runtime
// rewrites //localhost to //host.docker.internal on its way out of the container.
//
// Usage:  node dev/mcp.mjs up | down | build | status

import {execFileSync, spawnSync} from 'node:child_process';
import * as path from 'node:path';

const MCP = path.join(import.meta.dirname, '..', 'dockerfiles', 'mcp-server');
const NAME = 'grokky-dev-mcp';
const PORT = process.env.GROKKY_DEV_MCP_PORT ?? '33003';
const IMAGE = 'datagrok/grokky-mcp-server:dev';

// stdio:'inherit' makes execFileSync return null, so only trim when we actually captured output.
const docker = (args, opts = {}) => {
  const out = execFileSync('docker', args, {encoding: 'utf8', stdio: opts.inherit ? 'inherit' : 'pipe', ...opts});
  return out?.trim() ?? '';
};

const exists = () => {
  try {
    return docker(['ps', '-aq', '-f', `name=^${NAME}$`]).length > 0;
  } catch {
    return false;
  }
};

function build() {
  console.log('mcp: building image from the working tree …');
  docker(['build', '-t', IMAGE, '.'], {cwd: MCP, inherit: true});
}

function up() {
  build();
  if (exists())
    down();
  docker(['run', '-d', '--name', NAME, '-p', `${PORT}:3003`,
    '--add-host', 'host.docker.internal:host-gateway', IMAGE]);
  for (let i = 0; i < 40; i++) {
    const r = spawnSync('docker', ['exec', NAME, 'node', '-e',
      'fetch("http://127.0.0.1:3003/health").then(r=>process.exit(r.ok?0:1)).catch(()=>process.exit(1))']);
    if (r.status === 0) {
      console.log(`mcp: ${NAME} up on http://localhost:${PORT}/mcp`);
      console.log(docker(['logs', NAME]).split('\n').find((l) => l.includes('listening')) ?? '');
      return;
    }
    Atomics.wait(new Int32Array(new SharedArrayBuffer(4)), 0, 0, 500);
  }
  console.error('mcp: did not become healthy in 20s — check `docker logs ' + NAME + '`');
  process.exit(1);
}

function down() {
  if (!exists())
    return console.log('mcp: not running');
  docker(['rm', '-f', NAME]);
  console.log('mcp: removed');
}

const cmd = process.argv[2] ?? 'status';
if (cmd === 'up') up();
else if (cmd === 'down') down();
else if (cmd === 'build') build();
else if (cmd === 'status')
  console.log(exists() ? docker(['ps', '-f', `name=^${NAME}$`, '--format', '{{.Names}}  {{.Status}}  {{.Ports}}']) : 'mcp: not running');
else {
  console.error(`mcp: unknown command '${cmd}'`);
  process.exit(1);
}
