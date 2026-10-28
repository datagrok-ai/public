# Sharing: what the features could not say yet

Gap round of 2026-10-21 over `TestTrack/Sharing`, `playwright-public/Sharing` (10 specs) and the sharing
steps of other old specs (Projects, Chat, Browse, Widgets, EDA, PowerPack). The features were written with
today's steps only; what they need beyond that is listed here, as Gherkin where a scenario is ready.

Covered by this round: `sharing-connection` (its Sharing pane after a share), `sharing-layout` (the
second account applies a shared layout), `sharing-access-levels` (Full access, All users), `sharing-chat`
(two accounts in one chat), and the second account applying a shared model in EDA `models/share-model`.
Who may see, change, re-share or lose a shared entity is not a UI claim: it is
`ApiTests/src/dapi/permissions.ts`, which this round extended (§7). A script or a space shared and
taken back has no BDD feature for that reason.
Projects (Olena's area) keeps its own parked items: raising a grant to Full access and the
"Inherited from" row (`projects-lifecycle-files`, `projects-move`).

## 1. The "will also be shared" notice of the Share dialog

The dialog lists the other entities of the entity's project ("This entity will also be shared to view and
use:", `permissions_browser.dart`). A query saved through the JS API on the feature's own connection
(`a query "X" on "Admin:ShareConn1" reads …`) showed no such notice on localhost (2026-10-21): the query
has no project relations. The old specs only logged the notice as a remark.

Needed: a fixture whose project holds another entity (a query saved through the UI, or a script inside a
project), and then:

```gherkin
When user picks "Share..." from the context menu of "ShareConnQ1" gallery card
Then "Share ShareConnQ1" dialog should contain text "will also be shared to view and use"
And "Share ShareConnQ1" dialog should contain text "ShareConn1"
```

Question: which entities are expected to show it — a query and its connection, a script in a project?

## 2. The level the Sharing pane names for a connection

core#977 is fixed: right after a View and use share the Sharing pane of a connection says "can view and
use". `sharing-connection` guards the label (`@realizes:core#977`, no longer a known failure).

## 3. Steps the Share dialog lacks

- `When user removes {string} from {element}` (added) — the Remove of any grant row (a group such as
  All users), beside `user removes the sharing user from` which matches the sharing user only.
  `sharing-access-levels` now removes the All users grant through the dialog and checks the Sharing pane no
  longer lists it.
- `When user sets the access level of the sharing user in {element} to {string}` — open the selector of an
  existing grant row and check a level; only the new-grant selector (`share access selector`) is reachable.
- A level is picked by its **check box** (`user checks "Full access" tree node inside privilege tree`); a
  click on the name only selects the row. The note in `spaces/spaces-share-and-edges` that "clicking a
  level did not change the selector" describes a click on the name, not a platform failure.

## 4. Chat readings (`sharing-chat`)

The md's steps 2, 6-10 need readings of the Chats pane that no step has:

- who wrote a message (`the message "X" in Chats section should be by the sharing user`);
- the order of messages (`"Reply from B" should come after "First comment from A" in Chats section`);
- a message's own menu (edit, delete) and that another account's message offers neither;
- a mention highlighted and its notification.

## 5. Facts found on the way

- A layout saved through the JS API (`user saves the layout of the current table view as …`) is not
  offered in the Layouts pane; `sharing-layout` saves through the pane's Save and renames in the panel.
- A bare table shared on its own (not inside a project) and opened by the second account from Shared with
  me is not covered; a table inside a shared project is (`projects-lifecycle-*`).

## 6. What every run leaves on the server, and no step can remove

A query, script or connection saved through the API gets its own entity project (`projects.is_entity`);
deleting the entity leaves that project behind, empty, with the grants the Share dialog made on it. The JS
API lists none of them, and `DELETE /api/projects/<id>` answers "ok" and keeps them (localhost
2026-10-21: 140 such projects since 2026-09-23, from every feature that saves a query, script or
connection, sharing or not, e.g. ShareDialog1, BddScriptBrowser*, BDD-Conn-Browser-*). Grants on deleted
projects and spaces stay as well (rows in `permissions` whose entity is gone). Candidate platform bug,
awaiting a hand check; until the platform removes an entity's project with the entity, the features
cannot leave the stand as they found it.

## 7. API tests

Added to `ApiTests/src/dapi/permissions.ts` by this round (localhost 2026-10-21, the category 15/15):

- the owner keeps every privilege after every grant is revoked (old `sharing-all-users-and-owner` C,
  `sharing-edge-cases` 2);
- a space shared for viewing cannot be changed, deleted or shared further by the grantee (old
  `share-files-and-spaces-permissions` D);
- a connection shared on its own shares none of its queries (deleted TestTrack
  `share-connection-permissions` D.3).

Still to write: a folder share gives read but not write, and files outside the folder stay closed (old
`share-files-and-spaces-permissions` A/B). A script has no Execute privilege of its own (its type lists
only View, Edit, Delete and Share), so View is what lets the grantee run it, and the existing View and
use test covers it.

## 8. Edge cases (`share-dialog-edge-cases`) — still open

Written for the edge-case feature, with no library, binding or code change. That feature claims a removal
cancelled, and OK disabled until a recipient is accepted (core#980, a known failure: an email becomes a
recipient only once Enter accepts it).

### Steps requested

1. **A disabled user is not offered in the suggestions** (decided). Needs a set-up step that disables an account
   made for the feature (`a disabled user named {string} is on the server`) and a reading of the suggestions
   (`the user suggestions should not offer {string}`, beside the existing `the user suggestions should not offer
   the current user`). Then:

   ```gherkin
   When user types "<disabled user's name>" into "User, group, or email" input in "Share ShareEdgeA{time}" dialog
   Then the user suggestions should not offer "<disabled user's name>"
   ```
2. **A grant revoked while the recipient has the entity open** (decided): the revoked privilege is closed on the
   recipient's next save, refresh or reopen of the entity, whichever privilege it was. Needs the recipient's session kept open while the owner revokes,
   which one browser page cannot do (`user signs in as the sharing user` replaces the running session): a second
   page for the sharing user, or a step that revokes a grant through the API while the recipient's page stays open.

### Won't do

- **What a View-and-use recipient may do with a shared query or space** beyond re-sharing and putting it into a
  view-only space: the server's answer is ApiTests `permissions.ts` (§7); `sharing-recipient` claims the UI side.
- **Escape on the Share dialog**: it closes as CANCEL does (`Dialog._closeOnEsc`); `sharing-dialog` claims CANCEL
  and the one branch of its own, Escape closing the list of levels first.
