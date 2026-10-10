# Platform findings

Platform issues found while building Project Portal Core, with the workaround the package uses. Ticket numbers
are given where one was filed. Items 1 to 19 were first listed in the package README.

## Issues and workarounds

1. **Sharing a domain row is refused (GROK-21041).** A promoted domain row cannot be shared through
   `grok.dapi.permissions.grant`, even by the row's author holding Share (`~can_share` true).
   `privileges_service.savePermission` checks Share against the projects repository, which never holds domain rows.
   Program rows therefore stayed private to their author: a portfolio viewer saw no programs. **Workaround since
   2026-10-06:** each program row is linked into its workspace folder. The server lets whoever may view a space view
   the domain rows linked into it (View only; it clamps that path, so no write follows), so the row reaches the
   portfolio's groups and the program's roles through the folder (`src/program-rows.ts`). A program without a folder
   stays private, and Edit on a row still needs the fix: a contributor cannot write the master tables (gate log,
   evidence, documents) of a program they did not author. One package test pins the grant itself. Workspaces are
   not affected: a space takes the grant, so each portfolio's workspace is isolated today.
2. **`string_list` columns can't be written from JS.**
   - `Types.fromValue` reports any `List` as `list`, so a JS array never validates against a `string_list` column.
   - The workaround is to hold stage spellings in `stage.spellings` as a `|`-separated string.
3. **The npm-installed `grok` lags the repo.** It rejects current `domain-schema.schema.json` manifests, Grit's included. Use the workspace build at `node_modules/.bin/grok`.
4. **`dapi.connections.filter('name = "X"')` matches the friendly name.** `findConnection` in `src/portal.ts` falls back to an exact match.
5. **A range filter over an all-zero numeric column drops every row.** A column of `[0, 0, 0]` in the filter panel leaves 0 rows passing; `[5, 5, 5]` keeps all three. The portfolio table leaves a numeric column out of the filter panel until its values vary.
6. **The bar chart sorts categories by name, not by the column's category order.** Sorting by category ignores `setCategoryOrder`, while the box plot respects it. The table's bar chart splits on a hidden `stage in order` column (`01 D1`, `02 D2`).
7. **Listing a space's linked runs failed once (GROK-21043, closed as cannot reproduce).** In a test run on 2026-10-01,
   `GET /spaces/{id}/children?includeLinked=true` failed with `Converting object to an encodable object failed` once
   the space linked a run. Rechecks the same day did not reproduce it: existing program and stage folders holding linked
   runs, and a fresh space with a run just linked, list them (with and without a `FuncCall` types filter), and the
   folder view shows them. The ticket was closed; reopen it with a repro if it recurs. The tests still read
   `grok.dapi.projects.find(id).links`.
8. **A panel that takes a function is not offered for one.** A `meta.role: panel` function whose input is a `func`
   (or an `object`), with a condition that holds, does not appear in the context panel of a model script. The
   model card is therefore reached from the program page, through `modelCard(model)` and in the Workspace tab on
   Home (`workspaceModelCard`, a Power Pack role), not as an info panel.
9. **A package cannot schedule its own functions.** The server's scheduler (`#schedule:` with an optional
   `#schedule.runAs:` group) runs only server-side scripts and queries, and the JS API has no way to create a
   schedule. A package's JavaScript function does not run headless either. So `syncDuePortfolios` exists but is not
   registered; it runs when someone, or something with a browser session, calls it.
10. **No JS API to send a notification.** `POST /users/notifications` creates one, but only for the caller unless they
   hold Edit users, and the stand's `grok.dapi` has no notifications client at all. The portal calls the route with
   `fetch`.
11. **A column named `user` breaks the schema deploy.** The deployer does not quote column names in its DDL, so
   `"user": {"type": "user"}` fails with `syntax error at or near "user"`. `follow` names it `user_id`.
12. **A row of a table without soft delete cannot be deleted.** On a table declared `"softDelete": false`, both
   `delete(id)` and a transaction's delete answer `not-found` for a row that `query` returns. Changing `softDelete`
   afterwards is refused without a migration. `follow` keeps rows and turns them off instead.
13. **`FuncCall.edit()` runs a call whose inputs are all set** instead of showing its dialog. The portal's Run builds
   the dialog with `buildEditor`.
14. **A package function registers under its TypeScript name.** A different `//name:` header is ignored. The
   function documented as `programData` was therefore registered as `programDataTable`; since 0.3.0 the TypeScript
   function is named `programData`, and `programDataTable` stays as a second function for scripts that use it.
15. **Deleting a user setting that is not there throws.** `grok.userSettings.delete(name, key)` fails with
   `Cannot read properties of null` when the key is absent. The portal deletes only a key it has read.
16. **`in (...)` over a string column of uuid-shaped values fails.** `space_id in ("…", "…")` on `program.space_id`
   (a string) answers `operator does not exist: text = uuid`; `space_id = "…"` works. The pinned programs are read
   with `=` joined by `or`, 50 at a time.
17. **A `json` column holds a map, not a list.** Writing a JSON array to a domain `json` column fails with
   `Cannot convert '[...]' to map`. `observation.evidence` keeps its list under `items`.
18. **A table view opened from a click in the app view does not take the focus.** `grok.shell.addTableView` called
   in a click handler on the program page adds the table, but the app view stays current. The DMTA section's
   Open status defers the call with `setTimeout`.
19. **A domain schema redeploy never revokes a grant.** Removing a group or a permission from a table's `grants`
   in `schema.json` leaves the grant on a server that has it; changing `securityMode` or `defaultRowVisibility`
   of a deployed table is refused. Narrowing access is a runtime `revoke` by someone with Share on the table.
20. **A soft reference to `Core.func_calls` cannot be written (GROK-21042).** Inserting a `model_run` with
   `func_call` set fails with `column p.func_id does not exist`: the server's existence check for the soft
   reference aliases its target `p`, and the permission subquery it embeds declares its own `permissions p`, so
   `p.func_id` resolves against the wrong table. Saving runs and linking them into spaces works; only the record's
   reference waits. `saveModelRun` writes `run_id` as text and tries the soft reference first, so it switches over
   by itself once fixed. One package test pins it. Migration debt: backfill `func_call` from `run_id` with SQL.
21. **`Property.units` is null for a script parameter that declares `units:`.** The value is only in
   `property.options.units`.
22. **A domain `query({filter})` without `limit` returns the server default of 100 rows, with no warning.** Every
   list read in the package passes a limit (`src/rows.ts`).
23. **A datetime column comes back as a dayjs object whose `String()` is `toUTCString()`** (weekday first), so
   string sorts are wrong. Times are compared as times (`src/model/time.ts`).
24. **Every space lists its own storage connection** (`Files`, a `DataConnection`) among `children`, so "is this
   folder empty" has to skip it. Probably by design.
25. **A space's files have no list call.** A space file's full path (`Space:Program:Files/name`) reads through
   `grok.dapi.files`, but `SpaceFilesClient` has no list, so a folder's files come from `children.ofTypes('FileInfo')`.
26. **A property set on an app view object does not come back** through `grok.shell.views` or `grok.shell.v`,
   which hand out new wrappers. The portal keeps its views by root element.
27. **A click inside Home makes Home the current view again after the handler runs,** so a handler that switches
   views must do it after the click (the My programs tile waits 50 ms).
28. **`grok publish` without `--release` makes a per-user debug version.** Only the publishing user sees it; a
   test run as another user tests the released version.
29. **A function description holding `;` inside `{description: …}` breaks publishing** (`Missing ":" in option`).
   Known parser behaviour.
30. **The schema type for a floating-point column is `float`;** `double` fails `grok api` validation, although
   `double` is what the platform's property types call it.
31. **A column that is `immutable` cannot also declare `onDelete: "setnull"`** (manifest validation,
   `immutable-setnull`). By design.
32. **Inside `asUser` (impersonation) `DG.User.current()` stays the session's user,** while the server answers as
   the impersonated user. Tests that check identity read it from the server.
33. **`grok test` against a local server can stop starting** (Puppeteer "frame got detached", "Waiting failed", or a
   hang on the preloader), for any package. The in-browser `test` function still works and was used to count
   results when it happened.
34. **Json columns return keys in the server's order,** not the order written. Profile diffs sort keys.
35. **An admin does not read another user's row in a row-secured table with no row grants.** A `following` row
   written by a plain user is invisible to the admin session (query and count). Row security on a plugin table
   has no admin bypass in practice, so "only its author" really is only its author. The portal's consequence: a
   sync no longer updates other users' follows; each follower's news comes from their own visit.
36. **A table-secured table cannot be scoped per row, and its mode cannot change.** Reads of a `table` table check
   one grant on the table entity; only `row` tables support promotion, and so only their rows can be linked into a
   space. With item 19 (no `securityMode` change on a deployed table), narrowing who reads a table's rows means a
   new row-secured table and moving the rows (0.2.33: `approval`, `sync_log`, `following`).
37. **A package cannot give its domain rows an icon of their own, except through a JS object handler.** The domain
   schema has no `icon` per table; the platform draws every domain row with `fa-table`. A favorited row comes back as a
   bare `EntityRecord`, which Favorites and the Power Pack Workspace do not turn into a domain row, so a handler for
   the table's type never sees it there. The portal's handler claims the bare record by its name and gives programs the
   flask; it works only once the package has loaded (Home loads it through My programs). In a Space's tree a linked
   row always shows `fa-table`: the tree calls the platform's own domain-row meta, not the registered handler. Asks:
   an `icon` per table in the domain schema; favorites and the Workspace resolving a domain `EntityRecord` to its row
   (as Spaces do); the Space tree using the registered handler.

Items 38 to 45 were found building the project discussion and Ask SME on the platform's chats (F-51, 2026-10-09).
Each says what is missing, what the portal does instead, and what the platform would need.

38. **No js-api for chats.** The server has `/chats` routes (list by entity, comments, post, watch, delete), but
   `grok.dapi` has no chats client and `dapi2.chats` is generated empty. The portal calls the routes with `fetch` and
   the session token (`src/discussion.ts`), relying on the server's JSON shape. A new chat's id is the client's only
   when the comment carries its own id too; without one the server makes an empty chat. Need: `grok.dapi.chats`
   (find by entity, read, post, watch) and a `ui` chat widget a package can embed.
39. **A chat post emails people.** Every post to a public chat makes the server email its participants and watchers,
   and the first post of a chat emails every administrator, wherever email is configured (`email.dart`,
   `COMMENT_POSTED`). The portal must notify in the platform only. On the stand email has no server, so nothing goes
   out; on a deployment with email, a discussion post or an Ask SME question would email. The portal cannot turn this
   off per chat. Need: a chat (or a caller) that notifies in the platform only, or a per-user setting that the
   portal can rely on.
40. **No threads, anchors or recipients in a chat.** A chat hangs on entities and is one flat list; the Chats panel
   shows the first chat of an entity only. The portal keeps one chat per program folder and writes the anchor
   (criterion, run, document) and the asked users into a first line with `pp:` tokens, which the platform shows
   as text. Need: a reply-to or topic on a comment, and a reference to another entity (or a domain row) on a comment.
41. **No mentions, and no notification to the people asked.** Mention notifications are commented out in the server.
   `POST /users/notifications` creates a notification for someone else only for a user manager ("Permission denied
   to access notifications of another user" for a viewer). So a viewer's or a contributor's question is posted, and
   the people asked are not told; an admin's question tells them. A user cannot make someone else watch a chat
   either. Need: @mentions that notify in the platform, or a notification API that a package may call for the
   members of a portfolio's groups (item 10).
42. **Anyone who reads a chat can rewrite or delete any comment in it.** Posting a comment with an existing id
   replaces its text and keeps its author (proved: a viewer rewrote the admin's comment); `DELETE
   /chats/comments/{id}` checks nothing. A discussion that a decision cites can be changed after the decision. The
   portal offers no edit or delete, and a decision cites the thread's address, not a copy. Need: edit and delete by
   the author (and admins) only, and an edit history.
43. **A post is not checked against the chat's entities.** Reads are filtered by View on an attached entity, but a
   post only checks a private chat's group: a user who knows a public chat's id may post to it without seeing it.
   Ids are not listed to such a user, so the risk is small. Need: posting checked like reading.
44. **A comment cannot be back-dated.** The server sets `posted` to now. The example threads are dated when seeded,
   not in the example history's past. Need: an admin import with a date (as the gate log allows).
45. **A chat on a domain row needs the row promoted.** Criteria, runs and documents are domain rows; a chat on one
   would need it to be an entity, and promotion is per table. The portal anchors to them inside the program
   folder's chat instead (item 40). Need: chats on domain rows without promotion, or comment references (item 40).
46. **Deleting a shared entity leaves its entity project behind, and no one can delete that project.**
   `grok.dapi.permissions.grant(entity, group)` on a standalone entity (an uploaded table, a connection) goes
   through `Entity_Share` → `GET /projects/for_entity/{id}`, which creates an entity project named after the entity
   (`is_entity`, `is_dashboard`, linked to it) and shares that. When the entity is deleted, its relation goes, the
   project stays. It is listed in Home search under the entity's name and fails to open ("Can't open project
   EXSKIN1"). `grok.dapi.projects.delete` refuses it to everyone, admin included ("Insufficient privileges to delete
   Project"), since the permission check finds nothing left to check against. Repro, as admin:
   ```js
   const df = DG.DataFrame.fromColumns([DG.Column.fromList('int', 'x', [1])]); df.name = 'TEST-ORPHAN';
   await grok.dapi.tables.uploadDataFrame(df);
   const t = df.getTableInfo();
   await grok.dapi.permissions.grant(t, someGroup, false);       // one project "TEST-ORPHAN" appears
   await grok.dapi.tables.delete(await grok.dapi.tables.find(t.id)); // the project stays
   await grok.dapi.projects.delete(project);                       // "Insufficient privileges to delete Project"
   ```
   Deleting the project first, while the entity exists, works and leaves the entity and its grant in place. A
   platform bug, not our misuse: the share is the documented call. Our tests deleted run tables after saving runs,
   so every full core run left about ten (by 2026-10-09: 1,071 "Fed-batch time course", 72 "EXS-KIN1", 72
   "Estimated solubility (example)", 30 scratch connections). **Workaround:** `deleteRunTable` and
   `deleteEntityProjectsOf` (`src/model-runs.ts`) delete the project first, found by
   `isEntity = true and relations.entity.id = "<id>"`; a full core run now leaves none. On 2026-10-09 the 73 "EXS-KIN1"
   orphans were removed from the stand's database (the same deletes the server's own delete runs); the other names
   were left for a decision. Need: deleting an entity deletes its entity project (or the
   project is deletable by its author and admins once empty).

## What the server enforces on writes

Read in the server source (domain repository, predicate builder, deployer) while securing the portal's workflow
tables. Package functions run in the user's browser, so a check in them is the app's, not the server's.

- **Table mode:** insert and update need Edit on the table entity; delete needs Delete. Per group, not per row.
- **Master mode:** a write needs Edit on the referenced row (the program), or Edit on the securing table.
- **Row mode:** the author always reaches their own rows; others need a row grant (refused today, item 1) or the
  table grant.
- **Column security:** a write to a column needs Edit on its property schema (`shareColumn` / `restrictColumn`).
- **Keeping actions:** sync, a workspace build, export and advance by data check for the portfolio's approvers in
  the function (`requireKeeper`), which is the app's check. The server refuses a non-author's program and work-unit
  writes on its own; building folders (contributors hold Edit on the root) and export (the portfolio table is
  readable by all) are guarded by the app alone (`src/model/keeper-actions.ts`).
- **Not enforcement:** custom table permissions (`permissions: ["approve"]`) are answered to the app, never
  checked on a write, so "approver of one portfolio but not another" is the app's check. Check constraints see the
  row, not the caller.
- **No service identity:** server-side scripts run as the caller; the scheduler's `runAs` covers scripts and
  queries, not package JavaScript. There is no supported "write as the package" path.
- **Deploys are additive:** a schema redeploy never revokes a grant, and refuses a `securityMode`,
  `defaultRowVisibility` or `delegate` change on a deployed table (item 19).

## What the server enforces on reads

- **Table mode:** View on the table, or nothing. One grant covers every row, so it cannot keep one portfolio's
  rows from another's.
- **Row mode:** the author; direct row grants (refused today, item 1); and whoever may view a space the row is
  linked into (View only). With `defaultRowVisibility: table`, View on the table shows every row, so such a table
  must carry no View grant.
- **Master mode:** whoever reads the referenced row (the program).
- **Not enforcement:** the app's filter to the user's portfolios (`ofMyPortfolios`, `myPortfolioIds`) and the
  catalogue scoping run in the browser. They hide rows the server would also hide, and the catalogue, which the
  server shows to all users.

### Platform asks that would close the remaining gaps

- Row grants on promoted domain rows (GROK-21041).
- Soft references to `Core.func_calls` (GROK-21042).
- Row-scoped custom permissions, or a per-row write rule for table-mode tables (for example, Edit only where
  `portfolio_id` is one of the caller's portfolios).
- A supported server-side write identity for a package.
- Scheduling a package function (item 9) and a notification API (item 10).
- Chats a package can use for project discussion: a js-api, in-platform-only notification, threads and references,
  mentions, author-only edits (items 38 to 45).
