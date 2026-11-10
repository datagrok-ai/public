# Projects features — what is still missing

## Still open

- core#546 — two projects or queries of one name in a namespace share an nqName (the `/p/` side is
  `projects-regressions-duplicate-name`); claim when fixed: ApiTests, the second save refused or renamed.
- core#607 — a deleted connection leaves its wrapper project and credentials row; claim when fixed: a datlas test.
- core#509 — a share right after `GET /projects/for_entity/<id>` made the wrapper project is refused at times;
  claim when fixed: a datlas privileges test repeating the pair.
- core#564 — after "Embed viewer" a reload signs in as Guest once the dialog's frame has loaded; needs a wait for
  the frame and a claim on the signed-in account.
- core#1162 — a View-and-use recipient's card menu offers Rename... when the permission check lands after
  `Menu.beforeShowTimeout` (d4 `menu.dart`); claim when fixed, as projects-move's recipient: the menu lists no Rename...

## Won't do

- Two sessions both saving the original project: the last save wins, as `dashboard.md` ("Versioning") documents.
- A long project name: the only bound is the 255-character `friendly_name` column.
- The owner deleting a project a recipient has open: two sessions at once, which one page cannot hold.
- Special characters in a URL parameter value: by decision.
- Data sync when the source changes (query deleted or replaced, a column gone, file deleted): the Data loading error
  is `projects-lifecycle-script`; the rebuild on reopen and GROK-19725's warning are ApiTests `Dapi: projects: data sync`.
- A viewer split by a column the data-synced source lost: the viewer's generic missing-column message.
- Name rules of the Save and Rename dialogs (empty, blank, digits only, "/"): prop_gen 'FriendlyName test'; each
  dialog's invalid mark is claimed once (`projects-regressions-names`, `spaces-rename`).
- CANCEL on the Save, Rename and Delete Project dialogs: every dialog writes on OK only.
- The owner's Save dialog of a reopened project holds its name: the GROK-17700 and GROK-20006 scenarios save under it.
- A wrong-typed value in a dashboard address (?lim=abc): xamgle 'URL params bind list elements'.
- Open local file before a save: the no-creation-script switch rule is claimed on a JS-API table
  (`projects-regressions-names`); the upload is browse's subject.
- The Rename validator's message text: cosmetic; the invalid mark and the disabled OK are claimed.
- Save original project disabled for a recipient (a Dart radio option): the recipient's "Copy of <name>" shows the
  mode (`projects-move`).
