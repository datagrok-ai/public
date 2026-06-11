---
feature: sharing
sub_features_covered:
  - sharing.entity-types.files-spaces
  - sharing.share-dialog
  - sharing.context-panel-pane
  - sharing.permissions-editor
target_layer: manual-only
coverage_type: regression
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/sharing/share-files-and-spaces-permissions.md
migration_date: '2026-06-11'
source_text_fixes:
  - typo-log-out-form-to-from-block-a
  - typo-log-out-form-to-from-block-b
  - typo-log-out-form-to-from-block-c
  - typo-log-out-form-to-from-block-d
candidate_helpers: []
unresolved_ambiguities:
  - block-c-step-1-cascade-notice-confirm-live
scope_reductions: []
related_bugs: []
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-06-11-sharing-migrate-02
    timestamp: 2026-06-11T15:00:00Z
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
    timestamp: 2026-06-11T12:00:00Z
    failure_keys: []
---

# Sharing & Permissions — File shares & Spaces

## Setup

- `Owner` — the tester's own account, logged into the primary browser session
- `Recipient` — a second, **non-admin** account
- **Files** — a **user file share** (e.g. My Files)
- **Spaces** — a **Space the tester owns** (e.g. a newly created space: Browse > Spaces; right-click the Spaces item and select Create Space...; Enter a unique name and click OK) containing a dataset.

## Scenarios

### Block A — Files are shared via the file-share connection (two-actor)

1. As **owner**, share the **file-share connection** behind the tester's share with
   `Recipient` at **View and use**.
   * Expected result: the connection is shared; **ListFiles** is part of view-and-use.
2. Log out from the `owner` session and log in as the `recipient`.
3. As **recipient**, browse the shared file share.
   * Expected result: the recipient can **list and open files** in the share, but
     cannot edit the connection or write files.

### Block B — File-share negatives, then revoke (two-actor)

1. Recipient attempts to **edit the connection / delete the share / re-share it**.
   * Expected result: all **denied** — view-and-use grants browse/read only.
2. Log out from the `recipient` session and log in as the `owner`.
3. As **owner**, open the Share dialog and click the **×** next to the recipient's grant; click **OK**.
   * Expected result: the recipient grant is removed; the Sharing pane shows
     owner-only again.

### Block C — Share a Space (namespace project) cascades to its contents (two-actor)

1. As **owner**, set the **Space** as the current object and open its **Sharing**
   surface; click **SHARE...** (a Space is a Project, so it uses the standard Share
   dialog).
   * Expected result *(confirm live)*: the **Share <space>** dialog opens with the
     **View and use / Full access** dropdown and the cascade notice **"This entity will
     also be shared to view and use:"** listing the space's contained entities.
2. As **owner**, add `recipient` name and click **OK**.
   * Expected result: the space's grant + cascaded view-and-use on its contents are
     created.
3. Log out from the `owner` session and log in as the `recipient`.
4. As **recipient**, open **Spaces** / **Shared with me**.
   * Expected result: the shared space appears for the recipient and its contained
     entities are **viewable/usable** (tables load, queries run, layouts apply) — the
     namespace-level project cascade. (Same mechanism as the Project/Dashboard case.)

### Block D — Space negatives, then revoke (two-actor)

While the recipient holds **View and use** on the space (before revoke):
1. As **recipient**, attempt to **edit** the space or its contents, **delete**, and
   **re-share**.
   * Expected result: all **denied** — view-and-use grants no Edit / Delete / Share on
     the space or its members.
2. Log out from the `recipient` session and log in as the `owner`.
3. As **owner**, **revoke** the space share last.
   * Expected result: the space and its cascaded contents disappear from the
     recipient's browse.

## Notes

- `target_layer: manual-only` per atlas `manual_only[]` entry for `sharing.entity-types.files-spaces`: two-actor test (owner + non-admin recipient) requires two independent browser sessions; cannot be reduced to deterministic single-session UI automation.
- Block C step 1 is marked *(confirm live)* — live verification required to confirm the cascade notice actually appears for a Space entity in the current build. See `unresolved_ambiguities`: `block-c-step-1-cascade-notice-confirm-live`.
- Chain: `sharing` (revision 1). `pyramid_layer: manual` per chain analysis. `ui_coverage_delegated_to: null` — flow-open-share-dialog and flow-sharing-context-panel are orphaned in this chain pending the sharing-ui-smoke.md consolidation file.
- Atlas source: `sharing.entity-types.files-spaces` — `core/client/xamgle/lib/src/commands/file/share_dataset.dart#L7`. A file share is a DataConnection (dataSource: Files); Space sharing follows the project permission model.
- Source-text fixes applied: four occurrences of "Log out form" corrected to "Log out from" (Blocks A, B, C, D step 2).

---
{
  "order": 10,
  "datasets": []
}
