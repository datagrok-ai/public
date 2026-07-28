#!/usr/bin/env node
// Seeds the dev runtime container's Claude subscription credentials.
//
// Copies ONLY the `claudeAiOauth` block out of the developer's real
// ~/.claude/.credentials.json into dev/.claude/.credentials.json, which is what the dev
// container bind-mounts. Two reasons this is a copy and not a mount of the real file:
//
//   1. The real file also holds unrelated MCP OAuth tokens (Slack, etc). The container has
//      no business seeing those.
//   2. The Claude CLI refreshes the access token in place. Mounting the developer's live
//      file would make the container and the developer's own Claude Code session share one
//      credential, and refresh-token rotation can then log one of them out mid-session.
//
// The copy is mounted read-write so the container renews its own token indefinitely. When
// the refresh token finally expires, re-run this script.

import * as fs from 'node:fs';
import * as os from 'node:os';
import * as path from 'node:path';

const SRC = process.env.CLAUDE_CREDENTIALS ?? path.join(os.homedir(), '.claude', '.credentials.json');
const DEST_DIR = path.join(import.meta.dirname, '.claude');
const DEST = path.join(DEST_DIR, '.credentials.json');

function fail(msg) {
  console.error(`seed-creds: ${msg}`);
  process.exit(1);
}

if (!fs.existsSync(SRC))
  fail(`no credentials at ${SRC}. Log in with \`claude auth login\` first, or set CLAUDE_CREDENTIALS.`);

let creds;
try {
  creds = JSON.parse(fs.readFileSync(SRC, 'utf8'));
} catch (e) {
  fail(`could not parse ${SRC}: ${e.message}`);
}

const oauth = creds.claudeAiOauth;
if (!oauth?.accessToken)
  fail(`${SRC} has no claudeAiOauth.accessToken — this is not a subscription login.`);

fs.mkdirSync(DEST_DIR, {recursive: true});
fs.writeFileSync(DEST, JSON.stringify({claudeAiOauth: oauth}, null, 2));
fs.chmodSync(DEST, 0o600);

const dropped = Object.keys(creds).filter((k) => k !== 'claudeAiOauth');
const when = (ms) => !ms ? 'unknown' :
  `${new Date(ms).toISOString().slice(0, 16).replace('T', ' ')}Z (${Math.round((ms - Date.now()) / 36e5)}h)`;

console.log(`seed-creds: wrote ${DEST}`);
console.log(`  subscription : ${oauth.subscriptionType ?? 'unknown'}`);
console.log(`  access token : expires ${when(oauth.expiresAt)}`);
console.log(`  refresh token: expires ${when(oauth.refreshTokenExpiresAt)}`);
if (dropped.length)
  console.log(`  dropped from copy: ${dropped.join(', ')}`);
if (oauth.refreshTokenExpiresAt && oauth.refreshTokenExpiresAt < Date.now())
  fail('the refresh token has already expired — run `claude auth login` and re-seed.');
