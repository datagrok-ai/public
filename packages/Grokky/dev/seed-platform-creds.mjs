#!/usr/bin/env node
// Copies the subscription credentials into the PLATFORM-spawned claude-runtime container.
//
// Why this exists: the Dockerfile seeds ~/.claude/.credentials.json with a literal `{}`. Real
// credentials only get in there via the in-app `claude auth login` flow, which writes into the
// *running container's* filesystem — not into the image and not into any volume. So anything that
// recreates the container (rebuilding the image, `docker rm`, the platform respawning an
// on-demand container) silently resets it to `{}`.
//
// A credential-less runtime does not fail loudly. Every full-prompt turn returns in ~1 s with zero
// tokens and no error, which reads like "the model answered instantly and did nothing" — the exact
// symptom that looks like a model or prompt regression and is neither.
//
//   node dev/seed-platform-creds.mjs           seed every running platform claude-runtime
//   node dev/seed-platform-creds.mjs --check    report status only, change nothing
//
// Source is dev/.claude/.credentials.json (written by dev/seed-creds.mjs, which copies only the
// claudeAiOauth block and leaves unrelated MCP tokens behind). Run that first.

import {execFileSync} from 'node:child_process';
import * as fs from 'node:fs';
import * as path from 'node:path';

const SRC = path.join(import.meta.dirname, '.claude', '.credentials.json');
const DEST = '/home/grok/.claude/.credentials.json';
const checkOnly = process.argv.includes('--check');
// --force overwrites containers that already hold credentials. Needed when the stored copy's
// refresh token has been rotated dead (e.g. the containers holding the live session were deleted):
// "has credentials" then just means "has a corpse", and the default skip becomes a dead end.
const force = process.argv.includes('--force');

const sh = (args) => execFileSync('docker', args, {encoding: 'utf8'}).trim();
const fail = (msg) => {
  console.error(msg);
  process.exit(1);
};

if (!fs.existsSync(SRC))
  fail(`${SRC} not found — run "node dev/seed-creds.mjs" first.`);

const oauth = JSON.parse(fs.readFileSync(SRC, 'utf8')).claudeAiOauth;
if (!oauth?.accessToken)
  fail(`${SRC} has no claudeAiOauth.accessToken.`);
if (oauth.expiresAt && oauth.expiresAt < Date.now())
  console.warn(`warning: access token expired ${new Date(oauth.expiresAt).toISOString()} — ` +
    `the container refreshes it itself, but only if the refresh token is still valid.`);

// Platform containers only; the local dev one already bind-mounts its credentials.
const containers = sh(['ps', '--format', '{{.Names}}'])
  .split('\n').map((s) => s.trim())
  .filter((n) => n.includes('claude-runtime') && n !== 'grokky-dev-runtime');

if (!containers.length)
  fail('No platform claude-runtime container is running. Open the AI panel once to spawn it.');

for (const c of containers) {
  let bytes = '?';
  try {
    bytes = sh(['exec', c, 'sh', '-c', `wc -c < ${DEST} 2>/dev/null || echo 0`]);
  } catch { /* container may not expose a shell; fall through to the copy */ }
  const seeded = Number(bytes) > 10;
  if (checkOnly) {
    console.log(`${c}: ${seeded ? 'has credentials' : 'EMPTY — needs seeding'} (${bytes} bytes)`);
    continue;
  }
  if (seeded && !force) {
    console.log(`${c}: already has credentials (${bytes} bytes), leaving it alone (--force to overwrite)`);
    continue;
  }
  // Write via stdin rather than `docker cp` so the plaintext never lands in a host temp file.
  execFileSync('docker', ['exec', '-i', '-u', 'root', c, 'sh', '-c',
    `cat > ${DEST} && chown grok:grok ${DEST} && chmod 600 ${DEST}`],
  {input: JSON.stringify({claudeAiOauth: oauth}), encoding: 'utf8'});
  console.log(`${c}: ${seeded ? 'overwritten (--force)' : 'seeded'}`);
}

if (!checkOnly)
  console.log('\nDone. The runtime reads the file per turn, so no restart is needed.');
