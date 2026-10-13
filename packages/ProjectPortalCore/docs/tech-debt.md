# Tech debt

What we chose to live with for now, and what removes it. Platform issues themselves are in
[platform-findings.md](platform-findings.md) (numbers below refer to it); open product decisions are in
[build-log.md](build-log.md). Add an item when a shortcut is taken; strike it when it is paid.

Status as of 2026-10-09.

## Startup and hooks

| # | Debt | Why we took it | Cost | Pays it off |
|---|---|---|---|---|
| 1 | **The package autostarts** (`portalAutostart`, `meta.role: autostart`). The whole package (about 1.5 MB) loads at every platform start, for every user. | The Model card pane and the program-record handler (icon, Open → program page) only exist once the package code runs. Without it they work only after something else loaded the package (usually My programs on Home). | Startup download and execution for users who never use the portal; our hooks run on every context panel and function run; a startup bug affects everyone. | Platform: honour `meta.role: panel` for a func input (finding 8) and let a package declare handlers/icons for its domain rows without loading (finding 37). Interim if the platform fix is not coming before a customer rollout: a tiny hooks-only package that autostarts and loads the core on demand. Kept for now by decision (2026-10-09). |
| 2 | **Model card pane is added by hooking the context panel's build**, the way Compute2 adds REST, not as a declared info panel. `modelCardPanel` keeps `meta.role: panel` for when the platform honours it. | Finding 8. | Depends on accordion events and pane names; breaks quietly if the platform changes them. | Finding 8 fixed: delete the hook, keep the panel function. |
| 3 | **Program records are recognised by name** in the object handler (`portal.program:` entity records), so the flask icon and Open work in Workspace and Favorites. | Favorites hand over a bare record, not the domain row (finding 37). | Fragile to naming; no icon in a Space tree for linked rows. | Finding 37 fixed. |
| 4 | **The PowerPack hook role `workspaceModelPanel` contains "model"**, so Model Hub's `like "%model%"` filter lists `workspaceModelCard` as a model ("Model card" tile). | Named on 2026-10-02 without Model Hub in mind. | A stray tile in Model Hub. | Rename the role in PowerPack and here (e.g. `workspaceEntityPanel`), one coordinated change. |

## Security and permissions

| # | Debt | Why | Cost | Pays it off |
|---|---|---|---|---|
| 5 | **Rows are shared by linking them into workspace folders** (program, approval, sync history), not by row grants. | GROK-21041: grants on promoted rows are refused (finding 1). | A row without a folder is invisible to non-admins; access follows folder sharing, not the portfolio's groups directly; `_Sync history` / `_Approvals` folders exist only to carry reads. | GROK-21041 fixed: grant per row to the portfolio's groups; folders go back to holding work only. |
| 6 | **Only a row's author and admins can write decisions, evidence and documents** on someone else's program row. | GROK-21041. | Approvers who are not admins cannot decide in the app. | GROK-21041 fixed. |
| 7 | **Still readable by every user:** `program_role`, `portfolio`, the model catalogue tables. | Table-secured tables cannot be scoped per row and their mode cannot change after deploy (finding 36). | Users see other portfolios' role assignments, names and catalogue entries (not their data). | Move them to row-secured tables like `approval` / `sync_log` and copy the rows, as `moveToSecuredTables` did. |
| 8 | **Some checks run only in the browser.** Package functions run in the user's session, so a role check in them is the app's. Workspace folder building and export have no server-side guard beyond what the platform enforces on the objects. | No server-side service identity for packages. | A determined user calling the API directly is limited only by platform permissions. | A server-side function identity, or write rules on the tables involved. |
| 9 | **Grants are per table, not per portfolio** where tables stay table-secured. | Finding 36. | Coarse. | Same as 7. |

## Data and migrations

| # | Debt | Why | Pays it off |
|---|---|---|---|
| 10 | **Old tables kept after the security move:** `approval_request`, `sync_run`, `follow`, the compound values on `compound`. Closed to all users; old follows were not copied. | Additive schema; no data loss while verifying. | Drop them (schema migration) once every deployment has run `moveToSecuredTables`. |
| 11 | **One-off migration functions kept in the package:** `moveToSecuredTables`, `moveRecordsOutOfRoot`, `securePortalTables`. | Needed once per deployment. | Retire them after the last deployment has run them; keep a note in the CHANGELOG. |
| 12 | **Folder favorites still count as pins** next to program-record favorites; they were not migrated. | Backward compatibility. | A one-off migration of folder favorites to record favorites, then drop the folder path in `pinning.ts`. |
| 13 | **`model_run` keeps no reference to the platform's call** where GROK-21042 blocks the soft reference (finding 20). | GROK-21042. | GROK-21042 fixed: write the reference; backfill. |
| 14 | **Sync runs only when someone runs it.** The package cannot schedule its own functions (finding 9). | Platform. | Platform scheduling for package functions, or a server-side job. |
| 15 | **No notification API from JS** (finding 10); in-platform notices go through what the stand offers. | Platform. | A notifications client in the js-api. |
| 16 | **Soft-deleted program rows can leave live entities** (EMS); the example stand has some. | Platform behaviour. | Platform fix, then a clean-up of orphans. |
| 28 | **Each table of a saved run gets an entity project** named after the table ("Fed-batch time course", a program code), because `uploadTables` shares it with the viewers and a share wraps a standalone entity in a project (finding 46). Deleting the table leaves that project behind, and no one can delete it then. The package never deletes run tables; its tests do, through `deleteRunTable` (project first). | Sharing the table is how a viewer reads a run's data; there is no share without the project. | Finding 46 fixed (the entity's project goes with it). Any future code that deletes run tables or shared connections must call `deleteRunTable` / `deleteEntityProjectsOf` first. |

## Discussion and Ask SME (on the platform's Chats)

| # | Debt | Why | Pays it off |
|---|---|---|---|
| 22 | **Chats are called through raw `fetch` on the server's `/chats` routes** (`src/discussion.ts`), not a js-api client. | No chats client in the js-api (finding 38). | A js-api chats client; swap the calls. |
| 23 | **A thread's anchor (criterion, run, document) is written into the first line of each message**, and a decision cites a thread by its address. | Chats have no threads or references (finding 40). | Platform threads/references on chats; migrate the first-line anchors. |
| 24 | **Ask SME from a non-admin notifies nobody**; the question is posted, but only user managers may notify others. Questions put to a person do not show in Waiting for me. | Finding 41. | Platform @mentions / notifications for other users; then add questions to Waiting for me. |
| 25 | **The app offers no edit or delete of posts.** | Delete is guarded by the server (author, or the `DELETE_COMMENTS` permission; checked 2026-10-12), but edit is not: anyone who can read a chat can rewrite any comment and it still shows the original author (finding 42, core#768). Offering only delete would be odd, so neither is exposed. | core#768 fixed; then author-only edit and delete with an "edited" marker (F-60). |
| 26 | **Every chat post emails participants and watchers, and the first post emails all admins**, wherever email is configured (finding 39). Not an issue on a stand without mail. | Platform behaviour. | Decide per deployment before going live with email; ask the platform for a per-chat or per-package switch. |
| 27 | **Seeded example threads carry the seed date.** | Comments cannot be back-dated (finding 44). | Cosmetic; leave. |

## Process and tooling

| # | Debt | Pays it off |
|---|---|---|
| 17 | **`grok publish` makes an admin-only debug version; users run the last `--release`** (finding 28). Forgetting to release leaves users on an old build. | A release checklist step, or a CI publish. |
| 18 | **`grok publish` never deletes files removed from the package** on the stand; stale copies had to be removed by hand. | Platform tooling fix; until then, note removed files in the release steps. |
| 19 | **Comments are parsed as function headers** ("Unknown header tag" warnings on publish). Harmless but noisy. | Reword the comments that start with `//` + word + colon. |
| 20 | **Tests that depend on stand state:** scratch programs left by aborted BDD runs, a flaky pinned-programs test, BDD timeouts when the stand is slow. `grok test` sometimes hangs on the preloader (finding 33). | Clean-up hooks that run even on abort; isolate pinned-programs; CI on a fresh stand. |
| 21 | **The package is named `ProjectPortalCore`** so it can sit beside an older draft package on the same instance. | Rename to `ProjectPortal` once the older draft package is retired everywhere it is deployed. |
