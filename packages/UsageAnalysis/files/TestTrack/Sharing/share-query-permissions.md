---
feature: sharing
sub_features_covered:
  - sharing.context-panel-pane
  - sharing.share-dialog
  - sharing.permissions-editor
  - sharing.notification
  - sharing.advanced-editor
  - sharing.entity-types.query
  - sharing.browse-shared-with-me
target_layer: manual-only
coverage_type: regression
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/sharing/share-query-permissions.md
migration_date: '2026-06-11'
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
related_bugs:
  - GROK-19403
gate_verdicts:
  d:
    verdict: PASS
    cycle_id: 2026-06-11-sharing-migrate-02
    timestamp: '2026-06-11T12:00:00Z'
    failure_keys: []
  a:
    verdict: PASS
    cycle_id: 2026-06-11-sharing-migrate-02
    timestamp: '2026-06-11T12:30:00Z'
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
---

# Sharing & Permissions — Query

## Setup

- `Owner` — the tester's own account, logged into the primary browser session
- `Recipient` — a second, non-admin account
- `Custom group` — a dedicated group containing the recipient, for the group-sharing block
- **Entity under test** — a DataQuery owned by the owner whose run depends on a
  System:Datagrok connection (go to Browse > Databases > Postgres, right-click the Datagrok
  connection, select New Query..., insert `select * from public.groups` and click SAVE).

## Scenarios

### Block A — Open the Share dialog from the Sharing pane (owner)

1. As **owner**, set the query as the current object (Browse > Platform > Functions > Queries
   and click the query, so the Context Panel shows it). Expand the **Sharing** pane in the
   Context Panel.
   * Expected result: the Sharing pane lists the current grant — "`owner` can view and use"
     and "You are the owner" — and shows a **SHARE...** button.
2. Click **SHARE...**.
   * Expected result: the **Share `query name`** dialog opens with: a **"User, group, or
     email"** autocomplete input, a level dropdown defaulting to **View and use**, the
     existing owner grant row (with a **×**), an **Advanced editor...** link, and
     **OK** / **CANCEL**.

### Block B — Share-dialog controls (owner)

1. Type a few letters of the recipient into the **"User, group, or email"** input.
   * Expected result: an autocomplete list suggests matching users/groups; selecting one
     adds it as a pill with its own level dropdown and a **×** to remove.
2. Observe the dependent-entity notice and the notification controls.
   * Expected result: because the query depends on a connection, the dialog shows **"This
     entity will also be shared to view and use:"** followed by the connection name
     (cascade). When adding a new recipient, a **Send notifications** checkbox and a
     message box ("Type in message here") are available.
3. Click **CANCEL**.
   * Expected result: the dialog closes; no grant is changed (Sharing pane still shows
     owner-only).

### Block C — Advanced editor per-permission matrix (owner)

1. Open the Share dialog again and click **Advanced editor...**.
   * Expected result: a permissions view opens (route `/permissions/<id>`) showing a
     **Group × Object** matrix. Columns are grouped under **Common** — **View**, **Edit**,
     **Delete**, **Share** — plus the entity-specific **Execute** column for a DataQuery.
     The owner's existing rows show **View** (and **Execute**) checked.
2. Inspect the bottom **"Type in user, role or group to add..."** row.
   * Expected result: a new group/user can be added and individual permission checkboxes
     toggled independently, with a **SAVE** action. Close without saving.

### Block D — Share "View and use" to the recipient; recipient gains read + execute (two-actor)

1. As **owner**, open the Share dialog, add the `Recipient` name (or the **`Custom Group`**),
   leave the level at **View and use**, leave **Send notifications** as desired, and
   click **OK**.
   * Expected result: the dialog closes with no error; the Sharing pane now lists the
     recipient with **"can view and use"** in addition to the owner.
2. Log out from the `owner` session and log in as the `recipient`. In the **recipient**
   browser session, open **My stuff > Shared with me** (or search for the query in
   Browse > Platform > Functions > Queries view).
   * Expected result: the recipient has access to the shared query.
3. As **recipient**, open the query and **Run** it.
   * Expected result: the query opens and executes successfully.

### Block E — Negative: with "View and use", recipient cannot edit, delete, or re-share (two-actor)

1. As **recipient**, attempt to **edit** the query — change its SQL / parameters and try
   to **save**.
   * Expected result: editing is **not permitted**; the query content is unchanged.
2. As **recipient**, attempt to **delete** the query (Delete action / context menu).
   * Expected result: delete is **not permitted**; the query still exists.
3. As **recipient**, open the query's **Sharing** pane and attempt to **re-share** it
   (the way you would as owner).
   * Expected result: there is **no usable SHARE... affordance** (absent / disabled), or
     any attempt to grant is **rejected**. The grant list is unchanged.

### Block F — Revoke access from the recipient; access lost (two-actor) — do this LAST

1. Log out from the `recipient` session and log in as the `owner`.
2. As **owner**, open the Share dialog and click the **×** next to the recipient's grant;
   click **OK** (alternatively, use the Advanced editor to uncheck the grants and SAVE).
   * Expected result: the recipient grant is removed; the Sharing pane shows owner-only
     again.
3. Log out from the `owner` session and log in as the `recipient`.
4. In the **recipient** session, refresh and look under **Shared with me** / try to open
   the query.
   * Expected result: the query no longer appears under Shared with me and the recipient
     can no longer open or run it (access revoked).

## Notes

- Atlas: `sharing.entity-types.query` is in `manual_only[]` — two independent browser
  sessions (owner + non-admin recipient) are required for Blocks D, E, F. Automated
  Playwright cannot satisfy this constraint.
- The dependent-entity cascade notice in Block B step 2 is query-specific: a DataQuery
  always depends on a DataConnection, so "This entity will also be shared to view and use:"
  always appears. Source: `sharing.permissions-editor` (display co-shared dependent
  entities for non-Project targets).
- The **Execute** column in the Advanced editor (Block C) is entity-specific to DataQuery
  (DataQuery.Execute permission). Source: `sharing.advanced-editor`,
  `sharing.permission-types`.
- GROK-19403 (project dependency cascade failure) is the known bug for this feature section.
  This scenario does not reproduce it directly (query, not project), but the dependency
  cascade notice behavior in Block B is structurally related.
- See `sharing.yaml` chain for consolidation proposal: `sharing-ui-smoke.md` is proposed
  to own single-actor (owner-only) coverage of flow-open-share-dialog,
  flow-sharing-context-panel, and flow-advanced-editor without a second browser session.

---
{ "order": 1, "datasets": [] }
