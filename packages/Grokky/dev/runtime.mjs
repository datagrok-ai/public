#!/usr/bin/env node
// Dev-loop control for a LOCAL claude-runtime container.
//
// This container is deliberately separate from the ones Datagrok's spawner manages
// (datagrok-dev_cvm_grokky-claude-runtime-*). It exists so runtime code can be iterated on
// without an image rebuild and without touching the platform's containers.
//
// Two bind mounts make the loop fast:
//   dist/    -> /app/dist          `npm run watch` on the host recompiles in ~1s; restart picks it up
//   plugin/  -> /app/plugin        skills are plain .md — edits need no restart at all
//   .claude/ -> /home/grok/.claude subscription creds (see seed-creds.mjs), read-write so the
//                                  CLI can renew its own access token
//
// Usage:  node dev/runtime.mjs up | down | restart | logs | status

import {execFileSync, spawnSync} from 'node:child_process';
import * as fs from 'node:fs';
import * as path from 'node:path';

const DEV = import.meta.dirname;
const RUNTIME = path.join(DEV, '..', 'dockerfiles', 'claude-runtime');
const NAME = 'grokky-dev-runtime';
const PORT = process.env.GROKKY_DEV_PORT ?? '5355';
const IMAGE = process.env.GROKKY_DEV_IMAGE ?? 'datagrok/grokky-claude-runtime:admin';

const docker = (args, opts = {}) =>
  execFileSync('docker', args, {encoding: 'utf8', stdio: opts.inherit ? 'inherit' : 'pipe', ...opts}).trim();

const exists = () => {
  try {
    return docker(['ps', '-aq', '-f', `name=^${NAME}$`]).length > 0;
  } catch {
    return false;
  }
};

// Docker on Windows wants forward slashes and a drive-letter path it can translate.
const mount = (host, container, mode = 'rw') =>
  ['-v', `${path.resolve(host).replace(/\\/g, '/')}:${container}:${mode}`];

function up() {
  const creds = path.join(DEV, '.claude', '.credentials.json');
  if (!fs.existsSync(creds)) {
    console.error('runtime: no dev credentials — run `node dev/seed-creds.mjs` first.');
    process.exit(1);
  }
  const dist = path.join(RUNTIME, 'dist');
  if (!fs.existsSync(path.join(dist, 'server.js'))) {
    console.error('runtime: dist/server.js missing — run `npm run build` in dockerfiles/claude-runtime.');
    process.exit(1);
  }
  if (exists())
    down();
  docker([
    'run', '-d', '--name', NAME,
    '-p', `${PORT}:5355`,
    ...mount(dist, '/app/dist'),
    ...mount(path.join(RUNTIME, 'plugin'), '/app/plugin', 'ro'),
    ...mount(path.join(DEV, '.claude'), '/home/grok/.claude'),
    // host.docker.internal is implicit on Docker Desktop but not on plain Linux daemons.
    '--add-host', 'host.docker.internal:host-gateway',
    IMAGE,
  ]);
  console.log(`runtime: ${NAME} up on ws://localhost:${PORT}/ws`);
  waitHealthy();
  restoreCliConfig();
}

// ~/.claude.json (CLI state, distinct from .claude/.credentials.json) lives OUTSIDE the bind
// mount, so recreating the container loses it and every turn then dies with "Claude
// configuration file not found" — which surfaces as a 90s watchdog kill, not as an auth error.
// The CLI backs the file up into .claude/backups/, which IS mounted, so restore from there.
function restoreCliConfig() {
  const r = spawnSync('docker', ['exec', NAME, 'sh', '-c',
    '[ -f /home/grok/.claude.json ] && exit 0; ' +
    'b=$(ls -t /home/grok/.claude/backups/.claude.json.backup.* 2>/dev/null | head -1); ' +
    '[ -n "$b" ] && cp "$b" /home/grok/.claude.json && echo restored']);
  if (r.stdout?.toString().includes('restored'))
    console.log('runtime: restored ~/.claude.json from backup');
}

function waitHealthy() {
  for (let i = 0; i < 40; i++) {
    const r = spawnSync('docker', ['exec', NAME, 'node', '-e',
      `fetch('http://127.0.0.1:5355/health').then(r=>process.exit(r.ok?0:1)).catch(()=>process.exit(1))`]);
    if (r.status === 0) {
      console.log('runtime: healthy');
      console.log(authLine());
      return;
    }
    Atomics.wait(new Int32Array(new SharedArrayBuffer(4)), 0, 0, 500);
  }
  console.error('runtime: did not become healthy in 20s — check `node dev/runtime.mjs logs`');
  process.exit(1);
}

// The server logs exactly one "Claude auth:" line at startup; surface it so a missing or
// expired credential is obvious immediately rather than as a 401 mid-benchmark.
function authLine() {
  try {
    const log = docker(['logs', NAME]);
    return log.split('\n').find((l) => l.includes('Claude auth:')) ?? 'runtime: no auth line in logs (yet)';
  } catch {
    return '';
  }
}

function down() {
  if (!exists())
    return console.log('runtime: not running');
  docker(['rm', '-f', NAME]);
  console.log('runtime: removed');
}

const cmd = process.argv[2] ?? 'status';
if (cmd === 'up') up();
else if (cmd === 'down') down();
else if (cmd === 'restart') {up();}
else if (cmd === 'logs') docker(['logs', '-f', NAME], {inherit: true});
else if (cmd === 'status') {
  if (!exists())
    console.log('runtime: not running');
  else {
    console.log(docker(['ps', '-f', `name=^${NAME}$`, '--format', '{{.Names}}  {{.Status}}  {{.Ports}}']));
    console.log(authLine());
  }
} else {
  console.error(`unknown command: ${cmd}\nusage: node dev/runtime.mjs up|down|restart|logs|status`);
  process.exit(1);
}
