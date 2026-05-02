---
title: "Server keys"
sidebar_position: 3
keywords:
  - encryption
  - JWT
  - rotation
  - server keys
  - key management
---

Datagrok uses two RSA key pairs to protect data at rest and to sign
sessions. **Encryption keys** seal the credential blobs that back your
data connections (see [Secrets Managers](data-connection-credentials.md)
for the complementary vault-delegation flow). **Signing keys** sign
the JWTs your sessions ride on. The **Keys** page lets administrators
create, rotate, move, and revoke these keys without restarting the
platform.

You'll find it under **Browse > Platform > Keys**. The page is gated
by the `AdminKeys` permission — by default only members of the
**Administrators** group can see it.

## What a key is

Every key has:

| Field        | Meaning                                                                                                                                                  |
|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Kid**      | Unique identifier. Embedded in JWT headers (`kid`) and stamped on each `credentials` / `connections` row, so the right key is used to decrypt each blob. |
| **Usage**    | `encryption` (seals credentials, connection parameters) or `signing` (mints session JWTs).                                                               |
| **Backend**  | Where the private key material lives — Local PEM file, AWS Secrets Manager, or GCP Secret Manager.                                                       |
| **Status**   | Lifecycle state — see below.                                                                                                                             |

### Status

| Status         | What it means                                                                                                                              | Used for new traffic? | Verifies existing traffic? |
|----------------|--------------------------------------------------------------------------------------------------------------------------------------------|-----------------------|----------------------------|
| `primary`      | The active key for new operations. Exactly one per usage.                                                                                  | yes                   | yes                        |
| `active`       | Available but not the current primary. Useful as a rotation target you've pre-staged.                                                      | no                    | yes                        |
| `rotating_out` | Phase-out in flight. Encryption keys: credentials are being re-encrypted to the new primary. Signing keys: existing JWTs verify until the retirement window elapses. | no | yes |
| `retired`      | No longer used. Material may be deleted at any time.                                                                                       | no                    | no                         |
| `revoked`      | Immediately invalidated. Tokens / credentials sealed with it are rejected within the cache TTL.                                            | no                    | no                         |

## The first time you visit

On the first start after upgrade, your existing `auth.pem` and
`datagrok.pem` files are auto-registered as legacy keys (`kid=auth`,
`kid=datagrok`). Other PEM files in settings storage register as
signing keys you can re-classify later via **Edit…**. Nothing changes
about how those keys are stored — they keep working from local PEM
until you choose to move them.

## Common tasks

Each card has a context-actions menu (the small **▾** in the title
bar of the right-side panel) with **Edit…**, **Rotate…**, **Move…**,
**Revoke…**, and **Delete**.

### Create a key

Use the **New key…** toolbar button. Choose:

- **Usage** — `encryption` for new credential storage, `signing` for
  new JWT signing.
- **Location** — `Local`, or any AWS Secrets Manager / GCP Secret
  Manager [credentials connection](data-connection-credentials.md)
  you've already configured.
- **Length** — RSA bit length (2048 is the default).

A freshly-created key is `active`. It only takes traffic once you
rotate the existing primary into it (or create the very first key for
that usage).

### Rotate a key

**Rotate…** transitions the current `primary` to `rotating_out` and
promotes a replacement to `primary`. New traffic immediately uses the
replacement; existing traffic keeps verifying against the rotating-out
key until drain completes.

- **Encryption keys** must rotate to a fresh `kid`. Datagrok re-
  encrypts every credential and connection sealed with the old key,
  in batches, in the background. The old key transitions to `retired`
  the moment the last row moves. The right-side **Rotation** pane
  shows rows remaining and an estimated next-tick run.
- **Signing keys** can either rotate to a freshly-created `kid` or
  re-target an existing `active` signing key. There's nothing to
  re-encrypt; the old key sits in `rotating_out` for a short
  retention window (default 15 minutes — JWT TTL plus a buffer for
  clock skew) so already-issued JWTs stay verifiable, then auto-
  retires.

The Rotation pane is read-only. Rotation always advances on the
internal scheduler (≈ 60 s tick) — surfacing the next-tick timestamp
is the deliberate alternative to a "Run now" button.

### Move a key between backends

**Move…** copies the key material to the chosen target backend
(Local, an AWS Secrets Manager connection, or a GCP Secret Manager
connection), verifies it loads from the new home, updates the
registry, and deletes from the source. The key's `kid` and any
in-flight rotation are unaffected — only its storage location
changes.

If the dialog says "No other storage locations are configured", add
an [AWS or GCP credentials connection](data-connection-credentials.md)
first, then re-open Move.

### Revoke a key

**Revoke…** has two modes:

- **Queued** (default, only meaningful for primary keys). The key
  goes into `rotating_out` and an automatically-promoted replacement
  becomes `primary`. JWTs signed by the source keep verifying until
  the retention window elapses; new JWTs are signed by the
  replacement. This is the safe choice — no user gets logged out.
- **Force** (`?force=true` from the API, or any non-primary key
  goes here directly). Sets the key to `revoked` immediately. Tokens
  signed with it stop verifying within the cache TTL (≈ 60 s).
  Credentials encrypted with it become inaccessible until rotated to
  another key. Service tokens (long-lived dev keys) break at JWT
  expiry; interactive sessions break almost instantly.

The Revoke confirmation surfaces user/service session counts so you
can tell who'll be inconvenienced before you click.

### Delete a key

**Delete** removes the registry row. The server-side endpoint refuses
to delete a key that's still referenced by any credential / connection
unless you pass `force=true` — then it deletes the material as well.

Delete is for tidying up `retired` rows. Don't use it on an `active`
or `rotating_out` key.

:::note Cloud backend retention

AWS Secrets Manager and GCP Secret Manager manage their own retention
windows — the platform doesn't actively delete cloud-stored material on
**Delete** or **Revoke (force)**. To remove the underlying secret
immediately, scrub it in the cloud console after dropping the registry
row. Local PEM files are deleted right away.

:::

## Permissions

Read access (the gallery, Identity / Usage / Rotation panes) is
admin-only via the standard auth check. Writes — create, rotate, move,
revoke, delete — additionally require the `AdminKeys` global
permission. New deployments grant it to the **Administrators** group
automatically; on existing deployments the upgrade migration grants
it to **Administrators** and the **Admin** group.

`AdminKeys` is also what gates the **Browse > Platform > Keys** node
itself.

## Storage backends

| Backend                       | Pros                                                              | Cons                                                                          |
|-------------------------------|-------------------------------------------------------------------|-------------------------------------------------------------------------------|
| **Local PEM**                 | Zero external dependencies. Survives anywhere Datagrok runs.      | Material lives in settings storage; you carry the backup story.               |
| **AWS Secrets Manager**       | Cloud-managed encryption, IAM-scoped access, AWS audit log.       | Requires an AWS Secrets Manager [credentials connection](data-connection-credentials.md) and IAM access to the secret. |
| **GCP Secret Manager**        | Same idea, on Google Cloud.                                       | Requires a GCP credentials connection.                                        |

Material on cloud backends is stored as the PEM string under a secret
named `datagrok/server-keys/<kid>`. You can inspect it from the cloud
console; you should not edit it there — round-trip everything through
the Keys page.

## When to rotate

- **Routine**: rotate on a cadence that matches your security policy.
  Quarterly is a reasonable starting cadence for encryption keys;
  signing keys can rotate more aggressively because the retention
  window is small. The platform doesn't auto-trigger rotation —
  schedule a job (cron, CI, etc.) that calls **Rotate…** or hits
  `POST /admin/keys/manage/{kid}/rotate`.
- **Compromise suspected**: revoke immediately. For a primary key,
  start with queued revoke unless the threat model demands instant
  invalidation.
- **Backend migration**: rotate when you're moving between Local /
  AWS / GCP, or when changing the AWS / GCP credentials connection
  the key uses.

## Troubleshooting

| Symptom                                                            | Likely cause                                                                                                                | Fix                                                                                              |
|--------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------|
| Key stuck in `rotating_out` for an encryption rotation             | Background re-encryption hasn't finished. Each row migrates one at a time on the scheduler tick.                            | Wait for the next tick (the **Rotation** pane shows the timestamp). The key auto-retires when remaining hits zero. |
| Key stuck in `rotating_out` for a signing rotation                 | The retention window hasn't elapsed yet (default 15 min).                                                                   | Wait. The **Time left** countdown shows when the next tick will retire it.                       |
| Move dialog says "No other storage locations are configured"       | No AWS / GCP credentials connection exists, or the only one is the key's current location.                                  | Add a connection under [Secrets Managers](data-connection-credentials.md), then re-open Move.    |
| Some users get "Session expired" right after a force revoke        | Their JWT was signed by the just-revoked key.                                                                               | Expected. Refresh logs them back in with a token signed by the replacement primary.              |
| Trying to delete a key throws "still referenced"                   | Credentials or connections still have `key_kid = <this-kid>`.                                                               | Rotate first. Force-delete only when you're prepared to lose access to those rows.               |
