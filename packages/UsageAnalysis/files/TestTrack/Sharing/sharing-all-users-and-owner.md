---
feature: sharing
target_layer: playwright
coverage_type: regression
priority: p2
realizes: []
realized_as:
  - sharing-all-users-and-owner-spec.ts
related_bugs: []
---

# Sharing & Permissions — All users (everyone) & owner-retains-access

Covers two boundary cases of the permission system: sharing with the built-in **All users**
group (share with everyone), and the **owner override** — the owner always retains access
even after every grant, including their own, is removed. Requires two browser sessions: the
owner and a second, non-admin recipient account.

## Setup

- `Recipient` — a second, **non-admin** account
- An entity **owned by the owner** (any standard shareable entity, e.g. a project created under the demog table).

## Scenarios

### Block A — Share with "All users" reaches every authenticated user

1. As **owner**, open the Share dialog for the project and add the **All users** group at **View and use**; click OK.
   * Expected result: the Sharing pane lists **All users — can view and use**.
2. Log out from the `owner` session and log in as the `recipient`. As **recipient**, open **Shared with me** and search for the entity.
   * Expected result: the entity is **accessible** to the recipient at **view and use**.
3. As **recipient**, attempt to **delete / re-share**.
   * Expected result: **denied**.

### Block B — Downgrade / revoke All users

1. Log out from the `recipient` session and log in as the `owner`.
2. As **owner**, remove the **All users** grant.
3. Log out from the `owner` session and log in as the `recipient`.
   * Expected result: the recipient **loses access** — the project is no longer accessible.

### Block C — Owner always retains access

1. Log out from the `recipient` session and log in as the `owner`. As **owner**, open the Advanced editor and remove **every** grant — including the owner's own rows — and SAVE (or attempt to).
   * Expected result: the **owner can still view and edit** the project. Verify the owner can re-open the project and re-share afterwards.
   * Unresolved: it is unclear whether the server blocks the removal of the owner's grant row, silently re-adds it after SAVE, or another mechanism preserves the owner's access.
2. Confirm a **non-owner** with no grant cannot access it at this point.
