# Projects features — what is still missing

Open project tickets that break data or leave a project that cannot be opened, and have no test yet.

The UI cases of this group are features in this folder: core#531
(`projects-regressions-closed-source.feature`), core#602 (`projects-regressions-selection-column.feature`),
core#546 (`projects-regressions-duplicate-name.feature`),
core#563 (`projects-regressions-move-same-name.feature`), core#296 (`projects-regressions-shared-tag.feature`)
and core#626 (`projects-regressions-rename-from-status-bar.feature`), all as `@known-failure`s, and the fixed
core#945 and core#567 as regression guards (`projects-regressions-quick-save.feature`,
`projects-url-parameters-list.feature`); core#628, core#627 and core#463
already are (`projects-complex`, `projects-lifecycle-spaces`, `projects-regressions-query`). The rest of the
GROK-20511 group (core#446) is covered too: core#569 and core#630 (the link after the parameter is unticked)
in `projects-url-parameters.feature`, core#446 in `projects-function-url-run.feature`, and the fixed core#566,
core#568 and core#630 (the alias right after the save) by the Source checks after a save under a new name
(`projects-url-parameters.feature`).

## API and server tests requested

Each of these is server or JS API behaviour: the check is a few lines in an API or server test and does not
depend on the UI. Every case was reproduced on localhost, core master `d9ff84a962`, 2026-10-20.

1. **core#976, a space deleted with a project linked into it.** `core/server/datlas/test/spaces/spaces_test.dart`:
   save a project, `spaces.id(space).addEntity(project.id, link: true)`, delete the space, then the project is
   still found by its id and keeps its namespace. Today it is deleted with the space. The suite has
   "remove linked entity from space" and "delete space with child spaces cascade", but no case for a linked
   entity surviving the delete of its space.
2. **core#546, two entities with one name.** `ApiTests/src/dapi/projects.ts` and `queries.ts`: saving a second
   project (and a second query) of the same name in the same namespace is refused or gets a distinct nqName.
   Today both save as `<namespace>:<Name>`. The UI side for projects (the `/p/` address opens neither, "Unable
   to get project asset") is `projects-regressions-duplicate-name.feature`; a query of a taken name, and the
   project built on it that reopens on the other query, have no test.
3. **core#607, a connection deleted with its wrapper project.** A datlas connections test: deleting a
   connection through `DELETE /connectors/connections/<id>` removes its credentials row and its wrapper
   project with the permissions on it; `DELETE /projects/<id>` on a wrapper project either removes it or is
   refused. Today the credentials row and the wrapper project stay, and the project delete answers 200 and
   removes nothing.
4. **core#509, a share right after the wrapper project is created.** A datlas privileges test: a call to
   `GET /projects/for_entity/<id>` that creates an entity's wrapper project, followed at once by
   `POST /privileges/permissions/user` sharing it, succeeds every time (run in a loop). Today the share is
   intermittently refused with "You don't have a permission to share this object"
   (`scheduled_function_test.dart` hits it).

## Steps requested

1. **The signed-in account is the one the feature started with, and the frame of the "Embed viewer" dialog
   has loaded (core#564).** The ticket's case is all UI, and every other step of it exists: demog.csv opened
   from Browse > Files > Demo, a scatter plot added from the toolbox, the table saved as a project through the
   ribbon's Save dialog, `user picks "General > Embed..." from the viewer menu of scatter plot viewer`, the
   "Embed viewer" dialog closed with CANCEL, `user reloads the page`. Two steps are missing:
   - a wait for the page framed in the dialog to load. The frame opens the project under the guest session
     minted for it, and that session replaces the user's only once the frame has loaded: CANCEL right after
     the dialog shows closes it before then, the reload keeps the user, and the scenario passes with the bug
     present (localhost `3896258354`, 2026-10-21; with the dialog left open for a few seconds the reload
     signs in Guest User);
   - a claim that the account behind the page is the one the feature started with. `the running account
     should be signed in` refuses when the feature has not signed in as another account first.
2. **A project not tagged on the server (core#296).** `the project named {string} should be tagged {string} on
   the server` has no negative form. `projects-regressions-shared-tag.feature` claims the error balloon only;
   with the step it would also claim that the tag the second account added was not saved (today the tag even
   shows in Context Panel > Details).
