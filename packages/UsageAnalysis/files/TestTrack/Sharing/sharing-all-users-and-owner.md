---
feature: sharing
sub_features_covered:
  - sharing.permission-types
  - sharing.share-dialog
  - sharing.advanced-editor
  - sharing.server.privileges-router
  - sharing.browse-shared-with-me
target_layer: manual-only
coverage_type: regression
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/sharing/sharing-all-users-and-owner.md
migration_date: '2026-06-11'
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities:
  - owner-grant-removal-mechanism-unspecified
scope_reductions: []
related_bugs: []
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-06-11-sharing-migrate-02
    timestamp: '2026-06-11T10:00:00Z'
    failure_keys: []
    review_round: 1
    claims:
      - check: A-STRUCT-MECH-01
        status: PASS
      - check: A-STRUCT-MECH-02
        status: PASS
      - check: A-STRUCT-MECH-03
        status: PASS
      - check: A-STRUCT-MECH-04
        status: PASS
      - check: A-STRUCT-MECH-05
        status: PASS
      - check: A-STRUCT-MECH-06
        status: PASS
      - check: A-STRUCT-03
        status: PASS
      - check: A-STRUCT-04
        status: PASS
      - check: A-LAYER-ALIGN-01
        status: PASS
      - check: A-CONT-01
        status: PASS
      - check: A-BUG-01
        status: PASS
      - check: A-MERIT-01
        status: PASS
      - check: A-MERIT-02
        status: PASS
  d:
    verdict: PASS
    cycle_id: 2026-06-11-sharing-migrate-02
    timestamp: '2026-06-11T08:00:00Z'
    failure_keys: []
---
# Sharing & Permissions — All users (everyone) & owner-retains-access

Cross-cutting suite for the two boundary targets of the permission system: the built-in **All users** group (share-with-everyone) and the **owner override** (the owner always retains access even after all grants are removed). Two-actor: owner + recipient (`selenium1`, non-admin), second window.

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
   * Unresolved: it is unclear whether the server blocks the removal of the owner's grant row, silently re-adds it after SAVE, or another mechanism preserves the owner's access. See `owner-grant-removal-mechanism-unspecified` in `unresolved_ambiguities`. This maps to the atlas edge_case: `sharing.advanced-editor` + `sharing.server.privileges-router` — "Owner always retains full access".
2. Confirm a **non-owner** with no grant cannot access it at this point.

## Notes

- `target_layer: manual-only` per atlas `manual_only[]`: the core assertion (recipient's access) requires two independent authenticated browser sessions (owner + non-admin recipient), which cannot be reduced to deterministic single-session UI automation.
- Block A step 1 exercises `sharing.share-dialog` (open Share dialog, add All users group) and `sharing.permission-types` (SHARE_WITH_EVERYONE global permission required for owner to share with All users).
- Block A step 2 exercises `sharing.browse-shared-with-me` (recipient opens Shared with me and finds the entity).
- Block C step 1 exercises `sharing.advanced-editor` (remove all grants via PermissionsView matrix) and `sharing.server.privileges-router` (server-side enforcement of owner retention).
- Source: atlas `critical_paths[cp-share-with-all-users]` (p2), `edge_cases[]` entry "Owner always retains full access" (`coverage_type: edge`), `derived_from: core/server/datlas/lib/src/routers/privileges.dart#L4`.
- Atlas `dep_lifecycle_ops` covering this scenario: `grant_permission` (Block A), `revoke_permission` (Block B), `open_share_dialog` (Block A step 1), `open_advanced_editor` (Block C step 1). All ops are `source_agnostic: true` — the project entity used as test vehicle is incidental; behavior is identical across all source classes.

---
{ "order": 13, "datasets": [] }
