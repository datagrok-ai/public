# General features — what is missing

The features translated from TestTrack `General/` (this folder, and Helm, Bio and EDA for their own cases) with the vocabulary the library has today
(no library, binding or core change). They pass on a local stand built from master. This file lists
what they could not say, what each missing piece must check, and where its signal would come from.
Every item is a gap of the library or of the core, not of the features.

## 1. Names on the login form (login.test.ts, 27 of 31 tests)

`loginInput` and `passwordInput` (`core/client/xamgle/lib/src/signup_login/signup_login.dart:61,63`)
are bare `<input>`s with a placeholder and nothing else; no kind reaches them (the `placeholder`
strategy looks for a host that contains the placeholder). Wanted: an `aria-label` ("Login or Email",
"Password") or a `name=` on each — then the `text input` kind resolves `"Login or Email" text input` —
or registered elements `login field` / `password field` (`#signup-login-fields input[placeholder=…]`)
in `bindings/platform/elements.ts`. The status label (:119, `label.signup-status`) wants a name too,
so an empty label can be claimed (`login status should have text ""`); `login-logout.feature` claims
it today as no "Login failed" shown.

```gherkin
Scenario Outline: Wrong credentials fail and keep the form
  When user signs in as the sharing user
  And user opens the address "/u"
  And user clicks on "Logout" link
  And user enters "<login>" into "Login or Email" text input
  And user enters "<password>" into "Password" text input
  And user clicks on "Login" button
  Then "Login failed" text should be visible
  And browse tab should be absent
  When user signs in as themselves again
  Then the running account should be signed in

  Examples:
    | login                 | password                    |
    | no_such_user_xyz      | wrongpass123                |
    | bddsigner             | wrongpass_123!              |
    | bddsigner             |                             |
    |                       | anypassword                 |
    | "   "                 | "   "                       |
    | 1234567890            | 1234567890                  |
    | notreal@notexist.xyz  | somepassword                |
    | x                     | y                           |
    | ' OR '1'='1           | ' OR '1'='1                 |
    | admin'--              | anypassword                 |
    | <b>admin</b>          | password                    |
    | пользователь          | пароль                      |
    | 😀user                | password                    |
    | no_such_user          | !@#$%^&*()_+-=[]{}\|;:,.<>? |

Scenario: The fields are empty when the form shows
  When user signs in as the sharing user
  And user opens the address "/u"
  And user clicks on "Logout" link
  Then "Login or Email" text input should have value ""
  And "Password" text input should have value ""
```

The 500-character login and password go in the same Outline. Three different wrong pairs in a row
replace today's three empty-field clicks (the label is then also cleared by `onInput`, :521, so each
"Login failed" is the new answer). The old null-byte test accepted any non-empty label and should
claim "Login failed" like the rest. The server answers "Login failed" for every internal failure
(`datlas/.../user_management_service.dart:1245`), and there is no lockout, so the matrix is safe.

## 2. A valid sign-in through the form, and Enter in the password field (login.test.ts L01, L02, L31)

```gherkin
Given the sharing user can sign in with a password
When user enters the sharing user's credentials into the login form
And user clicks on "Login" button
Then browse tab should be visible
And the sharing user should be signed in

When user enters the sharing user's login and presses Enter in the password field
```

A step that types `DATAGROK_SHARING_LOGIN` / `DATAGROK_SHARING_PASSWORD` and never prints them; the
Given skips where the password is unset (CI mints the second user's token; a dev key alone gives no
password). The session the form makes is ended at feature end (its own Logout, or `signs in as
themselves again`). The same step after a failed attempt covers L31.

## 3. The page title (login.test.ts L03)

```gherkin
Then the page title should contain "Datagrok"
```

`document.title`; the login page's comes from `core/client/xamgle/web/login.html:9`.

## 4. No browser alert from a typed script (login.test.ts XSS)

```gherkin
Given browser alerts are recorded
When user enters "<script>alert(1)</script>" into "Login or Email" text input
And user clicks on "Login" button
Then "Login failed" text should be visible
And the browser should not have shown any alert
```

The negative of `the browser should have shown the alert {string}` (`common/steps.ts`), read after a
positive in the same scenario so it is not read on the first tick. Needs section 1 too.

## 5. `user is logged in` on a signed-out page

A Background that runs after a scenario failed on the login page waits 180 s for the shell and
fails: `loggedIn` (`bindings/common/session.ts`) signs the started session back in only when the page
is in the shell. Wanted: when the page is not in the shell and the context has a started session, sign
it back in (`signInWithSession`) instead of loading `/`. Until then every login scenario signs the
running account back in as its own last steps (as `login-logout.feature` does).

Also for the login form: `{element} should be enabled` reads `aria-disabled`, native `disabled` and a few
d4/u2 classes, but the form disables itself by putting `.disabled` (`pointer-events: none`) on
`#signup-form` (`signup_login.dart:436`, `signup_login.css:395`), so a Login button stuck unusable reads
"enabled". Wanted (library): the state check notices an ancestor with `.disabled` or `pointer-events: none`.
And a server claim that a logout ended exactly its own session:
`Then the session the sharing user signed in with should be ended on the server` (the server ends only the
request's session, `user_management_service.dart:1722`), plus a reload after a sign-in made through the
form itself (first-login-case-ui.md 7), which needs section 2.

## 6. A user's name put back at feature end (profile-settings.md)

`profile.feature` puts the fixture's name back through the same dialog in its last steps; a run killed
in between leaves `Prof<time> Ile<time>` on `bddprofile`. Wanted:

```gherkin
Given the first and last name of user "bddprofile" are put back at feature end
Then the user "bddprofile" should have the name "Prof{time} Ile{time}" on the server
```

Read `firstName`/`lastName` with `grok.dapi.users.find` at the Given, write them back with
`grok.dapi.users.save` at feature end and read them back. The server claim pairs the UI one read after
the reload.

Also wanted: a name on the profile's own name element (`div.grok-user-profile-name`,
`core/client/xamgle/lib/src/views/user_profile.dart:353`), so the shown name is claimed there rather than
anywhere on the page (the view's tab title carries the same display name).

## 7. The profile picture (profile-settings-ui.md, photo)

```gherkin
Given the picture of user "bddprofile" is put back at feature end
When user signs in as "bddprofile"
And user opens the address "/u"
And user uploads "fixtures/avatar.png" through first "Edit" icon
Then the profile picture should show an uploaded image
And the user "bddprofile" should have a picture on the server
When user uploads "fixtures/avatar.jpg" through first "Edit" icon
Then the user "bddprofile" should have a different picture on the server than before
When user uploads "fixtures/not-an-image.txt" through first "Edit" icon
Then the user "bddprofile" should have the same picture on the server as before
```

- No UI removes a picture, and the picture file stays on the server (as a project's does): the
  put-back clears `user.picture` through the API and deletes the stored file.
- The avatar with a picture is `Icons.image('user', pictureUrl)` (`name="icon-user"`,
  `core/client/xamgle/lib/src/meta/grok_user_meta.dart:306`), without one an initials avatar; a claim
  on the profile's own avatar needs a name on `div.grok-user-profile-picture`
  (`core/client/xamgle/lib/src/views/user_profile.dart:724`), since the sidebar shows the same icon.
- A non-image is dropped without any signal (`if (img == null) return;`, `user_profile.dart:728`):
  the negative needs a balloon or an event from the core to be claimable at all.
- The file is chosen through `htmlOpenBinaryFile`: `user uploads … through {element}` exists for
  the downloaded file only; a fixture image needs the same chooser step with a project file. Fixture
  images (png, jpg, gif) and a text file in `bdd/fixtures/`.

## 8. Change password — a gate before the link (profile-settings-ui.md, password)

```gherkin
Scenario: A mismatched confirmation is refused
  Given the signed-in account has a password
  When user opens the address "/u"
  And user clicks on "Change password..." link
  Then "Change password" dialog should be visible
  When user enters "bdd-not-the-password-{run}" into "Current password:" input in "Change password" dialog
  And user enters "Bdd-New-{run}-a" into "New password:" input in "Change password" dialog
  And user enters "Bdd-New-{run}-b" into "Retype password:" input in "Change password" dialog
  And user clicks on OK button in "Change password" dialog
  Then an error balloon containing "Password confirmation doesn't match" should have been shown
  And "Change password" dialog should be visible

Scenario: A wrong current password is refused
  Given the signed-in account has a password
  When user opens the address "/u"
  And user clicks on "Change password..." link
  And user enters "bdd-not-the-password-{run}" into "Current password:" input in "Change password" dialog
  And user enters "Bdd-New-{run}" into "New password:" input in "Change password" dialog
  And user enters "Bdd-New-{run}" into "Retype password:" input in "Change password" dialog
  And user clicks on OK button in "Change password" dialog
  Then an error balloon containing "Failed to change password" should have been shown
  And "Change password" dialog should be visible
```

Every phrase exists except the gate. Without it, an account with no password gets `resetUserPassword`
and a reset mail instead of the dialog (`user_profile.dart:480-485`); the fixture users have none. The
gate reads `hasPassword` of the signed-in user (Dart `User.hasPassword`; the JS `DG.User` does not
expose it) and skips where it is false. The sharing user (`bddsigner` locally, `test2` on CI) has a
password, so the scenarios would sign in as it. The new password must never be the right one to change
to: a success revokes every other session of the account
(`datlas/.../user_management_service.dart:1764`).

## 9. The Table Manager's rows (table-manager.md, table-manager-ui.md)

What the manager publishes today: the grid on `.grok-tables-manager` (`name="viewer-Grid"`) is registered
as the **TableManager** widget (`DG.Widget.find(root).type === 'TableManager'`, a DartWidget), and its
`getWidgetStatus()` returns `{values: {}, hitAreas: {}}`. So `the "rows" reading of Grid viewer in "Tables"
dock panel` says "the viewer reports: no readings", and no `cell N of name` area exists to click,
right-click or drag over. Wanted (core, not a step): `TableManager` / `ItemsGrid`
(`core/client/d4/lib/src/common/table_manager.dart`, `items_grid.dart`) answers `getWidgetStatus` with its
inner grid's status (the grid's `rows`, `text of cell N of <col>`, `cell N of <col>`, `header <col>`), or the
widget registry gives the root to the grid. With it, these scenarios are writable with existing steps:

    Scenario: The manager lists every open table, in the order they opened
      When user presses Alt+T
      Then the "rows" reading of Grid viewer in "Tables" dock panel should be 3
      And the "text of cell 1 of name" reading of Grid viewer in "Tables" dock panel should be "demog"
      And the "text of cell 2 of name" reading of Grid viewer in "Tables" dock panel should be "smiles"
      And the "text of cell 3 of name" reading of Grid viewer in "Tables" dock panel should be "spgi-100"
      When user switches to the "smiles" table view
      And user closes the current view
      Then the "rows" reading of Grid viewer in "Tables" dock panel should be 2
      And the "text of cell 2 of name" reading of Grid viewer in "Tables" dock panel should be "spgi-100"
      When user presses Alt+T
      Then "Tables" dock panel should be absent
      And no errors should have been logged

    Scenario: A click on a row makes its table's view current and the table the current object
      Given the context panel is open
      When user presses Alt+T
      And user clicks on the "cell 2 of name" area of Grid viewer in "Tables" dock panel
      Then the "smiles" view should be current
      And the context panel should show "smiles"
      When user clicks on the "cell 1 of name" area of Grid viewer in "Tables" dock panel
      Then the "demog" view should be current
      And the context panel should show "demog"
      When user presses Alt+T
      Then "Tables" dock panel should be absent
      And no errors should have been logged

    Scenario: Open as table makes a table of the manager's list
      When user presses Alt+T
      And user picks "Open as table" from the context menu of the "cell 1 of name" area of Grid viewer in "Tables" dock panel
      Then the "rows" reading of Grid viewer in "Tables" dock panel should be 4
      And the table should have 3 rows
      And the table should have a column "name"
      And the table should have a column "rowCount"
      And the value of "name" column in row 2 should be "smiles"
      When user presses Alt+T
      Then "Tables" dock panel should be absent
      And no errors should have been logged

    Scenario: Several selected tables get one submenu for all of them
      When user presses Alt+T
      And user drags from the "cell 1 of name" area to the "cell 3 of name" area of Grid viewer in "Tables" dock panel holding Shift
      And user right-clicks on the "cell 2 of name" area of Grid viewer in "Tables" dock panel
      Then the open menu should list "3 tables"
      When user closes the context menu
      Then the context panel should contain text "3 tables"
      When user presses Alt+T
      Then "Tables" dock panel should be absent
      And no errors should have been logged

    Scenario: Show > All adds every attribute as a column and takes them away again
      When user presses Alt+T
      And user picks "Show > All" from the context menu of the "cell 1 of name" area of Grid viewer in "Tables" dock panel
      Then Grid viewer in "Tables" dock panel should have a "header rowCount" area
      And Grid viewer in "Tables" dock panel should have a "header colCount" area
      And the "text of cell 1 of rowCount" reading of Grid viewer in "Tables" dock panel should be "5850"
      And the "text of cell 3 of rowCount" reading of Grid viewer in "Tables" dock panel should be "100"
      When user right-clicks on the "cell 1 of name" area of Grid viewer in "Tables" dock panel
      And user hovers over "Show" menu item in context menu
      Then "rowCount" menu item in context menu should be checked
      When user closes the context menu
      And user picks "Show > All" from the context menu of the "cell 1 of name" area of Grid viewer in "Tables" dock panel
      Then Grid viewer in "Tables" dock panel should not have a "header rowCount" area
      And Grid viewer in "Tables" dock panel should have a "cell 1 of name" area
      When user presses Alt+T
      Then "Tables" dock panel should be absent
      And no errors should have been logged

Unverified until the status exists: the Shift-drag gesture (the manager has no row header) and what the
context panel shows for a list of tables. Checked by a probe on localhost (5 Oct, core a014812cc2): the
second Show > All leaves the manager with no column at all, `name` included (`showAll` sets
`currentProps = []` when every attribute was shown, `items_grid.dart:68`), while the md expects only the
attribute columns to go — GROK-17558 (reopened 5 Oct 2026); the last scenario above claims the md's
behaviour.

The manager's identity is not claimed either: `Grid viewer in "Tables" dock panel` resolves to any grid
docked under a panel titled Tables. The old spec checked the `.grok-tables-manager` class, which no kind
or registered element reaches. Wanted: a registered element `table manager` (`.grok-tables-manager
[name="viewer-Grid"]`, `core/client/d4/lib/src/common/table_manager.dart`), so the feature claims that
Alt+T shows the Table Manager, not merely a grid.

## 10. The order of the view tabs (tabs-reordering-ui.md; tabs-order-and-projects.feature)

Translated with today's steps (Pass 2, 5 Oct): the drag is `user drags "<name>" view tab to "<name>" view
tab`, and the order is read position by position as `Nth view tab should have text`. Measured on the
local stand: a tab dropped on the middle of another lands right after it, and dropping it on a tab of
another dock group (Browse, Toolbox) moves the view into that group. What is still wanted:

    Then the view tabs should be in the order {string}

The table tabs of the document strip, left to right, by their titles. The ordinals the feature uses count
every visible view tab in page order — the Browse and Toolbox pane tabs, then Home — so they hold only
because the Background pins the shell (simple mode off, Browse open, Toolbox shown): a stand or a later
setting that docks those panes elsewhere shifts every position. Signal: the DOM order of
`.tab-handle[name^="view-handle: "]` inside the document manager's tab list (dock_spawn
`tab/tab_host.dart`; `tab_handle.dart` `changeTabPosition` moves the element).

## 11. A folder's lifecycle on a file share (files-cache.md, files-cache-spec.ts)

The UI walk needs nothing new; the hard cleanup rule and the server reads are missing (no step removes a
folder now and at feature end, none reads a server path back). files-cache-spec.ts itself used only
`grok.dapi.files` — an ApiTests matter.

    Given no folder {string} is on the server, now and at feature end
    Then the file {string} should be on the server
    Then the file {string} should not be on the server
    Then the file {string} on the server should hold the text {string}

- the sweep: `grok.dapi.files.delete` on the path (recursive), now and at feature end, then poll
  `grok.dapi.files.exists` until false; plus the `{run}` family over an hour old in the parent
  (`isStaleFixture`, runtime/server.ts; list the parent with `grok.dapi.files.list(parent, false)` taking
  directories — the user's-files sweep in platform/workspace.ts skips them).
- the reads: `grok.dapi.files.exists` / `readAsText`, polled (the txt editor's SAVE does not await its
  write: `FileInfoMeta.writeText`, `features/file_editors.dart`).

    Feature: A folder and a file on a file share are created, edited, renamed and deleted
      Background:
        Given user is logged in
        And the browse panel is open
        And no folder "System:AppData/UsageAnalysis/BDD Folder cache test {run}, System:AppData/UsageAnalysis/BDD Folder cache test1 {run}" is on the server, now and at feature end

      Scenario: A folder is created, gets a file with text, both are renamed and the folder deleted
        When user picks "Create folder..." from the context menu of "Files > App Data > UsageAnalysis" tree node inside browse tree
        And user enters "BDD Folder cache test {run}" into Name input in "Create folder" dialog
        And user clicks on OK button in "Create folder" dialog
        Then the "Create folder" dialog should close
        And the file "System:AppData/UsageAnalysis/BDD Folder cache test {run}" should be on the server
        When user picks "Create file..." from the context menu of "BDD Folder cache test {run}" tree node inside browse tree
        And user enters "test.txt" into Name input in "Create file" dialog
        And user clicks on OK button in "Create file" dialog
        Then the file "System:AppData/UsageAnalysis/BDD Folder cache test {run}/test.txt" should be on the server
        When user double-clicks on "test.txt" tree node inside browse tree
        And user replaces the code of code editor with "Hello world!"
        And user clicks on SAVE button
        Then the file "System:AppData/UsageAnalysis/BDD Folder cache test {run}/test.txt" on the server should hold the text "Hello world!"
        When user picks "Rename..." from the context menu of "test.txt" tree node inside browse tree
        And user enters "test1.txt" into "File name" input in Rename dialog
        And user clicks on OK button in Rename dialog
        Then the file "System:AppData/UsageAnalysis/BDD Folder cache test {run}/test1.txt" should be on the server
        And the file "System:AppData/UsageAnalysis/BDD Folder cache test {run}/test.txt" should not be on the server
        When user picks "Rename..." from the context menu of "BDD Folder cache test {run}" tree node inside browse tree
        And user enters "BDD Folder cache test1 {run}" into "Folder name" input in Rename dialog
        And user clicks on OK button in Rename dialog
        Then the file "System:AppData/UsageAnalysis/BDD Folder cache test1 {run}/test1.txt" should be on the server
        And the file "System:AppData/UsageAnalysis/BDD Folder cache test {run}" should not be on the server
        When user picks "Delete..." from the context menu of "BDD Folder cache test1 {run}" tree node inside browse tree
        And user clicks on DELETE button in "Are you sure?" dialog
        Then the file "System:AppData/UsageAnalysis/BDD Folder cache test1 {run}" should not be on the server
        And "BDD Folder cache test1 {run}" tree node inside browse tree should be absent
        And no errors should have been logged

UI facts from the core (not probed): Create folder... / Create file... on a connection or folder
(`data_connection_meta.dart:379,409`), Rename dialog inputs "File name" / "Folder name"
(`file_info_meta.dart:287`), Delete... → "Are you sure?" with DELETE (`file_info_meta.dart:306`).

## 12. The file share's cache (files-cache-ui.md)

Ticking Cache on the shared Demo connection changes a connection every stand and feature uses, and nothing
reads the mappings back. Wanted: a Files connection the feature owns, and a read of its cache settings.

    Given a Files connection named {string} on the folder {string} is on the server
    Then the {string} connection should cache files
    Then the {string} connection should have a cache mapping for {string} invalidated on {string}
    Then the {string} connection should have no cache mapping for {string}

- the fixture: `DG.DataConnection.create(name, {dataSource: 'Files', dir: …})` over a folder from §3's sweep,
  deleted at feature end like `a {string} connection named {string} is on the server` (platform/steps.ts).
- the mappings: `dapi.connect.connections.getCacheMappings(connection.id)` (Dart,
  `data_connection_meta.dart:245`); the JS API has no wrapper (runtime/server.ts endpoint or a js-api method).
- md correction: "Cache..." is on a folder (or file) of a connection whose Cache is on
  (`file_info_meta.dart:33`), not on the connection; the connection's own Cache checkbox is in Edit... (Cache tab,
  `data_connection_meta.dart:791`).

    @serial
    Feature: A cached file share keeps a mapping per folder and drops it with the folder
      Background:
        Given user is logged in
        And the browse panel is open
        And no folder "System:AppData/UsageAnalysis/BDD Cache root {run}" is on the server, now and at feature end
        And a Files connection named "BDD-Cache-{run}" on the folder "System:AppData/UsageAnalysis/BDD Cache root {run}" is on the server

      Scenario: Cache is switched on in the Edit dialog
        When user picks "Edit..." from the context menu of "BDD-Cache-{run}" tree node inside browse tree
        And user clicks on "Cache" tab in "Edit connection" dialog
        And user checks "Cache" checkbox in "Edit connection" dialog
        And user clicks on OK button in "Edit connection" dialog
        Then the "BDD-Cache-{run}" connection should cache files

      Scenario: A folder's mapping with a cron is saved and goes when the folder is deleted
        When user picks "Create folder..." from the context menu of "BDD-Cache-{run}" tree node inside browse tree
        And user enters "Folder cache test" into Name input in "Create folder" dialog
        And user clicks on OK button in "Create folder" dialog
        And user picks "Cache..." from the context menu of "Folder cache test" tree node inside browse tree
        And user enters "*/2 * * * *" into text input in "Cache settings for Folder cache test" dialog
        And user clicks on OK button in "Cache settings for Folder cache test" dialog
        Then the "BDD-Cache-{run}" connection should have a cache mapping for "/Folder cache test" invalidated on "*/2 * * * *"
        When user picks "Delete..." from the context menu of "Folder cache test" tree node inside browse tree
        And user clicks on DELETE button in "Are you sure?" dialog
        Then the "BDD-Cache-{run}" connection should have no cache mapping for "/Folder cache test"
        And no errors should have been logged

## 13. Helm × Bio menu (Helm/bdd/features/integration/bio-menu.feature)

### 13a. A query drawn in the HELM substructure filter (General/helm-bio-menu-integration, Subsequence Search)
The feature claims only that Subsequence Search adds a filter on the HELM column. On a HELM column the
filter is `HelmBioFilter` (Bio/src/widgets/bio-substructure-filter-helm.ts): a bare
`ui.div('', {style: {cursor: 'pointer'}})` that opens the HELM editor in `ui.dialog({showHeader: false})`
on click — not the full-screen editor the Helm project's `HELM editor` element names
(`.d4-dialog-full-screen:has([data-testid="app-root"])`), so neither the click target nor the editor's
palette/canvas elements resolve. Wanted, in Helm's or the library's elements:

    When user clicks on HELM query of "HELM" filter in filters viewer
    Then HELM filter editor should be visible
    When user clicks on Peptides palette tab in HELM filter editor
    And user clicks on A monomer tile in HELM filter editor
    And user clicks on an empty spot of editor canvas in HELM filter editor
    And user clicks on OK button in HELM filter editor
    Then fewer than 53 rows should pass the filter
    And the filter should pass exactly the rows where "HELM" contains "A"

- a name on the filter's clickable panel (`_filterPanel`, e.g. `name="div-helm-filter-query"`);
- an element for the filter's editor dialog: the dialog that contains `[data-testid="app-root"]` and is
  not full-screen, with the same in-editor parts as `HELM editor`.

### 13b. The per-command "no Helm-related balloon" filter of the md (H17)

    Then no error balloon matching "helm|jsdraw|pseudo-?molfile|getmolfiles|cell\s*render" should have been shown

The negative of the existing `an error or warning balloon matching {string}` (viewers tier): no balloon
of kind error since the scenario's floor whose text matches. The md tolerates a UMAP warning on the
showcase; on the local stand no such balloon appeared, so the feature claims "no errors logged" only.

Since Pass 2 the feature claims every scenario stricter instead: no error or warning balloon of any
kind (`no error or warning balloon should have been shown`); the showcase showed no tolerated warning,
so the filtered negative is not needed yet.

### 13c. The cells of the search viewers' result grids

Similarity Search and Diversity Search copy the column's renderer to their own result grid
(`similar (HELM)`, `diverse (HELM)`; `Bio/src/analysis/sequence-similarity-viewer.ts:83`), which is the
Helm side of those two leaves. The viewers report no canvas of their own, so a pixel claim reads the
result grid's chrome, and no step reaches the inner grid's cells. Wanted:

```gherkin
Then the "cell type of similar (HELM)" reading of the result grid of "Sequence Similarity Search" viewer should be "helm"
Then the "cell 1 of similar (HELM)" area of the result grid of "Sequence Similarity Search" viewer should be painted in at least 3 colors
```

The inner grid is a `DG.Grid` inside the viewer root; a scope phrase that resolves a grid inside a viewer
(`DG.Widget.find` on the nested `[name="viewer-Grid"]`) would make the existing grid readings work there.

## 14. Predictive models (EDA/bdd/features/models/apply-and-delete.feature)

### 14a. The old spec's dataset
    Given user opens accelerometer dataset keeping the first 2000 rows

`sensors/accelerometer.csv` (System:DemoFiles, 80,616 rows, accel_x/accel_y/accel_z/time_offset) is not
registered in `libraries/bdd/bindings/platform/datasets.ts`, and `{dataset}` does not take a platform path
at compile time (`compile.ts:107` "dataset … is not registered", although the parameter's regexp and the
survey suggested it would). The feature trains on iris instead.

### 14b. Choosing a model in the Apply dialog by name
    When user selects the "BDD-Iris-PLS-{run}" model in "Apply predictive model" dialog

The Model choice lists `"<createdOn>: <friendlyName>"` cut at 40 characters
(`core/client/xamgle/lib/src/features/predictive_modeling/predictive_modeling_core.dart:31-36`), so
`user selects {string} in …` (exact label, `gestures.ts` selectNative) can never name a `{run}` model.
The feature walks the list with ArrowDown/ArrowUp and claims the chosen model in the context panel
(`fillOptions` sets `AppEvents.currentObject`). Wanted: a select step matching an option whose label
ends with / contains the name (or a name/value on the option, e.g. the model id, from core).
Also: the dialog's default choice does not become the current object until the choice changes
(the panel reads `Dialog "Dialog"` while it opens) — `fillOptions()` runs before the dialog shows and
the dialog takes the current object after; worth a core look, not a defect a user sees.

Which model OK applied: each prediction column carries the tag `predictive.model` = the model's markup,
which holds its friendly name and id (`core/client/xamgle/lib/src/features/predictive_modeling/engines.dart:96-103`).
The library compares a tag only for equality, so the feature claims the two prediction columns differ
instead (PLS trained with two components). Wanted:

```gherkin
Then the "predictive.model" tag of "Petal.Width (3)" column should contain "BDD-Iris-PLS-{run}"
```

### 14c. A generated table (old spec step 3, grok.data.testData random walk)
    Given user opens the "random walk" test dataset with 1000 rows and 10 columns
    Then the Inputs of "Apply predictive model" dialog should map "Sepal.Length" to "<column>"

`grok.data.testData(name, rows, cols)`; and a name on the dialog's ColumnsMapInput rows to claim the
mapping (predictive_modeling_core.dart:43). The feature applies to an inline table with the input
columns instead.

## 15. BiostructureViewer GROK-17485 — ad-hoc PDB survives a project round trip (not written)

Source: General/biostructureviewer-bug-grok-17485-spec.ts (copy in TestTrack/BiostructureViewer). The
PDB is a local package file (`System:AppData/BiostructureViewer/samples/1bdq.pdb`), no RCSB call. The old
spec's core claim only `console.warn`-ed on failure (`KNOWN_REGRESSION_SOFT_WARN = true`), so it never
failed. BiostructureViewer has no bdd project (needs `grok-bdd init` there — the lead's call) or the
feature goes to UsageAnalysis/bdd/features/general.

```gherkin
@journey
Feature: An ad-hoc PDB in the Biostructure viewer survives a project save and reopen (GROK-17485)
  Background:
    Given user is logged in
    And the "BiostructureViewer" package is installed
    And the user's own project "bdd-bsv-17485-{run}" is removed now and at feature end

  Scenario: The structure goes into the viewer
    Given user opens a table "bsv-17485" with:
      | sentinel |
      | 1bdq     |
    When user adds a Biostructure viewer
    And user sets "pdb" property of Biostructure viewer to the text of the "System:AppData/BiostructureViewer/samples/1bdq.pdb" file
    Then "representation" property of Biostructure viewer should be "cartoon"
    And the "structure loaded" reading of Biostructure viewer should be "true"
    And the "chains" reading of Biostructure viewer should be at least 1
    And no errors should have been logged

  Scenario: The project keeps the structure
    When user saves the current view as project "bdd-bsv-17485-{run}"
    Then 1 project named "bdd-bsv-17485-{run}" should be on the server
    When user closes all views
    And user opens the "bdd-bsv-17485-{run}" project and waits for its table
    Then Biostructure viewer should be visible
    And "pdb" property of Biostructure viewer should hold the text of the "System:AppData/BiostructureViewer/samples/1bdq.pdb" file
    And the "structure loaded" reading of Biostructure viewer should be "true"
    And no errors should have been logged

  Scenario: The project shared with the second account opens with the structure
    # share via Browse > Dashboards card menu "Share...", Send notifications off (share-model pattern)
    When user signs in as the sharing user
    And user opens the "bdd-bsv-17485-{run}" project and waits for its table
    Then "pdb" property of Biostructure viewer should hold the text of the "System:AppData/BiostructureViewer/samples/1bdq.pdb" file
    When user signs in as themselves again
```

Missing:
- `When user sets {string} property of {widget} to the text of the {string} file` — `grok.dapi.files.readAsText(path)`
  in the page, `viewer.props[name] = text`, then the viewer's render signal.
- `Then {string} property of {widget} should hold the text of the {string} file` — equality, first
  differing offset in the message (the old spec compared lengths).
- Package signal (Pass 0): `MolstarViewer` (BiostructureViewer/src/viewers/molstar-viewer/molstar-viewer.ts)
  publishes no `getWidgetStatus`/`isRenderPending` (only `onRendered`, ~1613-1623): readings
  `structure loaded`, `chains`, `atoms`, `representation` from the Mol* plugin's loaded structure.
  Without them "painted" reads a WebGL canvas the library cannot read.
- A step to open a project's Share dialog from its Dashboards card (share-model uses a gallery label;
  check the phrase for a project card).

## 16. Notebooks (General/notebooks-lifecycle-jupyter-container, Notebooks/playwright) — not written

Everything past the gallery runs the Jupyter container (initContainer, convertNotebook, JupyterLab) or is
API-only; the Notebooks browser view is already claimed in UsageAnalysis browse-platform-and-databases.
A menu-listing scenario was dropped as vacuous: `ML | Notebooks | New Notebook...` is registered by core
itself (`core/client/xamgle/lib/src/features/jupyter_notebook/jupyter_notebook_plugin.dart:7-14`) whether
or not the Notebooks package or the capability exists. What would mean something:

    Given the stand advertises the "NOTEBOOKS" capability
    Then the "ML > Notebooks > New Notebook..." top menu item should be enabled
    # and on a stand without it / without the package:
    Then the "ML > Notebooks > New Notebook..." top menu item should be disabled with the reason "This feature is not available on this server"

- menu item enabled/disabled state (+ the reason from `checkEnabled`, shown as its tooltip) — a top-menu
  reading of `aria-disabled` on the item;
- a capability gate on `grok.shell.startupData.fleetCapabilities` — a new gate kind, the lead's decision
  ("Nothing else skips").
- (only if wanted) `Given a notebook {string} is on the server` — `DG.Notebook.template({})` +
  `grok.dapi.notebooks.save`, deleted at feature end — for the gallery's search/Share legs.

## 17. Stand facts met (not gaps of the features)

- The Train Model view asks EVERY registered engine's `isApplicable` on each features change
  (`predictive_modeling_view.dart:537-546`) and `log.error`s a failure. Samples registers a Docker engine
  (PyKNN, `Samples/dockerfiles/app.py`); on the local stand Samples is installed without its container
  (`/api/docker/containers` lists only dbtests, bio, api-tests…, notebooks, chem-chem, chem-chemprop), so
  the server throws "Container is not started" (`action_logger_service.dart:72`) and the first scenario of
  apply-and-delete fails its error floor there. With Samples published with its container (5 Oct) the
  feature is green on localhost and on dev. (train-on-cars and share-model, also EDA, do not hit it on
  the same stand — why their feature changes do not reach PyKNN was not established.)

## 18. Toolbox Search: what a year alone means (toolbox-search.feature)

Not a missing step but an open expectation for GROK-20229: the feature claims `STARTED > 1990` reads the
year as its first day (5573 rows on demog, the same as `STARTED > 1/1/1990`). Reading it as "after the
whole year" would give 2674. Both differ from today's 0, so the known failure flips whichever way it is
fixed; if the fix chooses "after the year", the expected count in the feature changes to 2674.

## 19. How the reopened molecule column is drawn (molecule-in-exported-csv.md)

The md asks that the column be "correctly visualized as SMILES". The feature claims what the column is
(semantic type Molecule, units smiles, 100 distinct non-empty values, no MOLBLOCK in the file), not how the
grid draws it: no step reads a grid cell's renderer or what it painted for a molecule. Wanted:
`Then the "Structure" column of grid should be drawn by the "Molecule" renderer` (the grid's cell type
reading for the column) and a claim that a cell painted a structure (ink in the cell area beyond text).

## 20. Favorites after GROK-21108 (browse, spaces, users-groups-roles features)

GROK-21108 made "Add to favorites" a check: a plain item for an account that administers no group,
otherwise a submenu of "Only for me" and the groups it administers (`Favorites.buildMenu`,
`core/client/xamgle/lib/src/features/favorites.dart`); "Remove from favorites" is gone, and files are
favorites too. The four features that used it (browse-context-panel-and-menus, browse-my-stuff,
spaces-favorites, users-manage Users-21) now pick "Add to favorites > Only for me" to add and pick it
again to take out, and claim the result on the server or in Browse > My stuff > Favorites.

Missing:
- A step that reads whether a menu item is checked, e.g. `the open menu item "Add to favorites > Only
  for me" should be checked` / `... should not be checked` (the check mark of the d4 menu item). Without
  it the check the menu shows after adding is not claimed, only its effect.
- The path assumes the running account administers a group (admin does on dev, the local stand and CI).
  An account without one gets the plain item, and the submenu path fails. A step choosing the path by
  the account, or a precondition step `the account administers a group`, would make it explicit.

## Not translated by rule

- Login with Google (login-ui.md 1-4): Google's consent screen is an outside service.
- First login (first-login-case-ui.md 1-6): needs an account that has never signed in; a run spends
  it and users cannot be deleted. Step 7 (a reload keeps the account) is in `login-logout.feature`.
- Inactivity (inactivity-response-ui.md): a 20-minute wall-clock wait, then a query and a server
  script.
- Startup time (startup-time-ui.md): a performance threshold on cold cache, not a UI behaviour.
- network.md: no assertion ("inspect the network for redundant queries").
- profile-settings-spec.ts: a JS API round trip, nothing UI-specific (ApiTests).
- Chemprop (chemprop-spec.ts): trains and predicts in the chem-chemprop Docker container; the Train Model
  UI before it is covered by EDA's model features.
- Notebooks lifecycle (notebooks-lifecycle-jupyter-container): every step past the gallery runs the
  Jupyter container or is API-only (section 16).
