---
title: "Keypair authentication"
sidebar_position: 5
keywords:
  - keypair
  - authentication
  - developer key
  - grok login
  - CI credentials
---

Keypair authentication lets you and your CI jobs log in to Datagrok by proving you hold a
private key, instead of by sending a reusable secret. It replaces the
[developer key](#migrating-from-the-developer-key), which is deprecated.

## Why keys instead of a developer key

A developer key is one long-lived string that grants your full account. To use it you copy it
out of the UI and paste it into `~/.grok/config.yaml`, a Jenkins credential, a GitHub secret,
a container environment. Every copy is a place it can leak from, every copy has to be rotated
together, and the server cannot tell one copy from another.

With a keypair:

- The **private key never leaves the machine that generated it.** What the server stores is
  the public half, which is useless to an attacker who reads it.
- Logging in is a **challenge–response**: the server issues a single-use nonce, the client
  signs it, and the signature is worthless afterwards. Nothing replayable crosses the wire.
- Each machine gets **its own key**. Revoking your laptop does not disturb your CI job.
- Keys can have an **expiration date**, so an enrollment for a contractor or a one-off build
  agent stops working on its own.

## Log in from the command line

Install [datagrok-tools](../../develop/develop.md#packages) (`npm install -g datagrok-tools`),
then:

```
grok login https://dev.datagrok.ai
```

This generates an EC P-256 key and opens your browser, where you sign in the usual way
(password, SSO, SAML — whatever your deployment uses). The page asks for the **verification
code** the terminal is showing: type it, approve, and the CLI writes:

- the private key to `~/.grok/keys/<alias>.json`, readable only by you
- the server, alias and your login to `~/.grok/config.yaml`

From then on `grok publish`, `grok test`, `grok s` and the rest of the CLI authenticate with
it automatically. Nothing else to configure.

### Without a browser

On a headless box, or when the browser round trip is awkward, use a one-shot enrollment code:

1. In Datagrok, open your profile (click your avatar) and choose **Public keys...**
2. Copy the generated `grok login ... --code ...` command.
3. Paste it into the terminal on the machine you want to enroll.

The code is single-use and expires 15 minutes after it is generated. Generate a new one for
each machine.

### Options

| Option | Meaning |
|---|---|
| `--code <code>` | One-shot enrollment code from your profile. Skips the browser step |
| `--name <name>` | Key name shown in your profile. Defaults to `user-host` |
| `--expires <days\|date>` | `30`, or an ISO date such as `2027-01-31`. Default: never expires |
| `--alias <alias>` | Config alias to write. Defaults to the server's first host label |

## Manage your keys

Profile > **Public keys...** lists every key on your account with its name, fingerprint,
expiration and when it was last used, and lets you:

- **Revoke a key** — anything signing with it stops working immediately, including sessions it
  already opened.
- **Add a public key** you generated elsewhere — paste a JWK, or an RSA public key in PEM.
- **Generate an enrollment code** for another machine.

Keys marked *from config* were declared in the deployment configuration by whoever operates the
stand. You cannot remove them from the UI; ask your administrator.

## CI and automation

A build agent has no browser and no home directory worth trusting, so give it the private key
through the environment instead of a file:

```bash
export GROK_PRIVATE_KEY='{"kty":"EC","crv":"P-256","x":"...","y":"...","d":"..."}'
grok publish dev
```

`GROK_PRIVATE_KEY` accepts the JWK directly or base64-encoded, since some secret stores mangle
multi-line values. It takes precedence over anything in `~/.grok/config.yaml`.

To create a key for CI without enrolling an interactive machine:

1. Log in as the CI user (or ask an administrator to open that user's profile).
2. **Public keys... > Add a public key**, paste the public JWK, set an expiration.
3. Store the matching private JWK as the secret your pipeline exposes as `GROK_PRIVATE_KEY`.

Generate the pair with any tool that emits JWK; for example, with Node:

```js
const {generateKeyPairSync} = require('crypto');
const {publicKey, privateKey} = generateKeyPairSync('ec', {namedCurve: 'prime256v1'});
console.log(JSON.stringify(publicKey.export({format: 'jwk'})));
console.log(JSON.stringify(privateKey.export({format: 'jwk'})));
```

## Deployment-managed keys

An administrator can put keys directly in the deployment configuration, so a freshly deployed
stand already trusts the operator's automation:

```json
{
  "userKeys": [
    {
      "login": "ci",
      "name": "jenkins-build",
      "publicKey": "{\"kty\":\"EC\",\"crv\":\"P-256\",\"x\":\"...\",\"y\":\"...\"}",
      "expires": "2027-01-01"
    }
  ]
}
```

These are re-applied on every deploy and are marked *from config* in the user's profile.

## Migrating from the developer key

The developer key still works, and nothing breaks the day you start using keys. The CLI prefers
a keypair whenever one is configured for the server and falls back to the developer key
otherwise, so the migration is per-machine and reversible:

1. Run `grok login <server>` on the machine.
2. Confirm `grok publish <alias>` still works.
3. Clear the `key:` field for that server in `~/.grok/config.yaml`.

For CI, replace the dev-key secret with `GROK_PRIVATE_KEY` in the job's credentials. Once every
consumer of an account's developer key has moved, rotate the key one last time in
Profile > **Developer key...** to invalidate the copies you no longer control.

## How it works

1. The client asks for a challenge: `POST /users/login/key/challenge` with the key's
   fingerprint. The server stores a single-use nonce that expires in two minutes.
2. The client signs `datagrok-login:<api-url>:<nonce>` with the private key (ECDSA P-256 /
   SHA-256, raw `r||s`; RSA keys use PKCS#1 v1.5) and sends the signature to
   `POST /users/login/key`. The API URL is in there so a server you point the CLI at cannot
   relay your signature to another one where the same key is registered.
3. The server redeems the nonce — which identifies the user — checks that it is the server the
   signature names, checks the signature against the stored public key, and returns a normal
   session token.

A key is identified by its [RFC 7638](https://www.rfc-editor.org/rfc/rfc7638) JWK thumbprint,
so the client and the server derive the same fingerprint independently.

Keypair login needs **Datagrok 1.28 or later**. Against an older server the CLI says so and
names the version; if a developer key is still configured for that server, it is used instead.

## See also

- [Users and groups](users-and-groups.md)
- [Access control](access-control.md)
- [Package development](../../develop/develop.md)
