# General features — what is still missing

The features translated from TestTrack `General/` (this folder, and Helm, Bio and EDA for their own cases).
The first list of gaps was reviewed on 2026-10-07 against the library's rules (`public/libraries/bdd/CLAUDE.md`:
a feature claims what only the browser shows, what an API test can check goes there, and a feature is
fast). Section numbers are the first list's.

## Done

- **§5, a signed-out page.** `user is logged in` on a page a failed scenario left on the login form signs
  the context's own session back in (`bindings/common/session.ts`) instead of waiting 180 s for a shell.
- **§6, the fixture user's name.** `a user {string} is on the server` puts a fixture user back to its own
  name (the login, no last name) now and at feature end; `the user {string} should have the name {string}
  on the server` claims a save. `profile.feature` edits the name once and claims it on the page and the
  server, with no reloads.
- **§9, the Table Manager's rows.** `ItemsGrid` answers `getWidgetStatus` with its grid's status
  (`core/client/d4/lib/src/common/items_grid.dart`), so `table-manager.feature` claims the listing, a row
  click, Open as table and a closed table leaving the list.
- **§10, the order of the view tabs.** `the view tabs should be in the order {string}` reads the tab strip
  and counts only the named views.
- **§14b, choosing a model by name.** The Apply dialog keys each model by its id
  (`predictive_modeling_core.dart`) and `user selects the predictive model {string} in {element}` picks it,
  so `apply-and-delete` no longer walks positions and both EDA model features left the serial lane.
- **§15, the ad-hoc PDB round trip.** Covered by `biostructure-viewer/data-file-persistence.feature` and
  the Mol* viewer's status.
- **§18, the search syntax.** The cases are ddt matcher tests (`core/shared/ddt/test/data_frame/
  matcher_test.dart`, "search box: …"); the quoted value and the year alone are skipped there until
  GROK-20229 is fixed.

## Still open

- **"and" / "or" in the table view's Search box (GROK-20229).** The box parses its whole text with one
  matcher per column type (`xamgle features/search.dart`); the fix would route it through the row matcher
  the expression filter uses, with its test beside the fix.
- **§2, a valid sign-in through the login form.** Needs an account whose password the run knows; CI signs
  in by dev key. A lead's decision: a per-stand fixture account with a password, or none.
- **§7, the profile picture.** An image fixture, the file chooser step for a project file, a reading of the
  profile's own avatar (`user_profile.dart`), and a signal from the core when a non-image is dropped, which
  is silent today (`if (img == null) return;`).
- **§8, Change password.** The mismatched confirmation is the dialog's own check and is UI; it needs a gate
  on the signed-in account having a password (`User.hasPassword`, which the JS `DG.User` does not expose).
  The refusal of a wrong current password is the server's answer — a datlas test.
- **§13a, a query drawn in the HELM substructure filter.** A name on the filter's clickable panel and an
  element for its editor dialog (`Bio/src/widgets/bio-substructure-filter-helm.ts`).

## Not gaps

- **§1, the credential matrix.** Fourteen wrong, empty, boundary and injection pairs are one server answer
  each — the server's authentication, a datlas test; `login-logout.feature` claims the form's own part once.
- **§3, the page title; §4, an alert from a typed script.** No behaviour a person sees: the login text is
  never rendered back.
- **§5, the rest.** "The session ended on the server" is a datlas claim; the `.disabled` ancestor of the
  login form is a library refinement nothing needs yet.
- **§11, a folder's lifecycle on a file share.** The old spec used only `grok.dapi.files` (ApiTests
  `dapi/files`); the Browse commands for files would be a feature of their own, not a gap of this one.
- **§12, the file share's cache mappings.** Server configuration, read back through an endpoint the JS API
  does not wrap — a datlas test.
- **§13b.** Superseded: `bio-menu.feature` claims no error or warning balloon at all.
- **§13c, §19, what a renderer drew.** The package tests the renderer (the 2026-10-06 ruling); a feature
  claims the cell type and the value.
- **§14a, §14c, the accelerometer table and a random walk.** iris and an inline table show the same UI.
- **§16, the Notebooks menu state.** The menu item is registered by the core whether or not the capability
  exists; a new gate kind is not worth it.
- **§17, the PyKNN engine.** A stand fact, not a gap: the Train Model view asks every engine whether it
  applies and logs an error for an engine whose container is not started — a product finding for the
  view, not for the feature.
- **§20, the check mark of a menu item.** The server claim is the proof of Add To Favorites; the account
  administering a group is how every stand the suite runs on is set up.

## Not translated by rule

- Login with Google (login-ui.md 1-4): Google's consent screen is an outside service.
- First login (first-login-case-ui.md 1-6): needs an account that has never signed in; a run spends it and
  users cannot be deleted.
- Inactivity (inactivity-response-ui.md): a 20-minute wall-clock wait, then a query and a server script.
- Startup time (startup-time-ui.md): a performance threshold on cold cache, not a UI behaviour.
- network.md: no assertion ("inspect the network for redundant queries").
- profile-settings-spec.ts: a JS API round trip, nothing UI-specific (ApiTests).
- Chemprop (chemprop-spec.ts): trains and predicts in the chem-chemprop Docker container.
- Notebooks lifecycle (notebooks-lifecycle-jupyter-container): every step past the gallery runs the
  Jupyter container or is API-only.
