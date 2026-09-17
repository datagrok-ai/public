#!/usr/bin/env bash
# Prints a Datagrok session token.
#
# Usage: grok-token.sh <apiUrl> [devKey]
#
# With GROK_PRIVATE_KEY set (a private JWK, raw or base64-encoded) it signs a
# server-issued nonce, so nothing reusable crosses the wire and nothing reusable
# sits in the workflow environment beyond the key itself. Without it, it falls
# back to the developer key, which is deprecated.
#
# Self-contained on purpose: it runs in steps that have node but not yet
# datagrok-tools. See core/docs/KEYPAIR_AUTH.md.
set -euo pipefail

apiUrl="${1%/}"
key="${2:-}"

if [ -z "${GROK_PRIVATE_KEY:-}" ]; then
  exec "$(dirname "$0")/dev-login.sh" "$apiUrl" "$key"
fi

node - "$apiUrl" <<'NODE'
const crypto = require('crypto');
const apiUrl = process.argv[2];

const raw = process.env.GROK_PRIVATE_KEY.trim();
const privateKey = JSON.parse(raw.startsWith('{') ? raw : Buffer.from(raw, 'base64').toString('utf-8'));

const {d, p, q, dp, dq, qi, ...publicKey} = privateKey;
const members = {EC: ['crv', 'kty', 'x', 'y'], RSA: ['e', 'kty', 'n']}[publicKey.kty];
if (!members)
  throw new Error(`Unsupported key type "${publicKey.kty}"`);
// RFC 7638 JWK thumbprint - the key's id, derived identically by the server.
const canonical = '{' + members.map((m) => `${JSON.stringify(m)}:${JSON.stringify(publicKey[m])}`).join(',') + '}';
const fingerprint = crypto.createHash('sha256').update(canonical).digest('base64url');

const post = async (path, body) => {
  const res = await fetch(`${apiUrl}${path}`, {
    method: 'POST',
    headers: {'content-type': 'application/json'},
    body: JSON.stringify(body),
  });
  const text = await res.text();
  try {
    return JSON.parse(text);
  } catch {
    throw new Error(`Unexpected response from ${path} (HTTP ${res.status}): ${text.slice(0, 200)}`);
  }
};

(async () => {
  const {nonce} = await post('/users/login/key/challenge', {fingerprint});
  const options = {key: crypto.createPrivateKey({key: privateKey, format: 'jwk'})};
  if (publicKey.kty === 'EC')
    options.dsaEncoding = 'ieee-p1363';
  // The API root is signed too, so a server cannot relay this signature to a second stand.
  const signature = crypto.sign('sha256', Buffer.from(`datagrok-login:${apiUrl}:${nonce}`), options).toString('base64url');
  const login = await post('/users/login/key', {fingerprint, audience: apiUrl, nonce, signature});
  if (login.isSuccess !== true)
    throw new Error(login.comment ?? 'Key login failed');
  console.log(login.token);
})().catch((e) => {
  console.error(e.message ?? e);
  process.exit(1);
});
NODE
