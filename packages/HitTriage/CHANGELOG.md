# HitTriage changelog

## v.next

* Hit Triage, PepTriage: A cloned template's name and key are free, as in Hit Design: "<name> (copy)" ("(copy 2)" and on when that is taken) and the key's next free number
* Hit Triage: Fixed the grid and its filters not showing what "Calculate additional properties" wrote into columns the campaign already had
* All apps: Fixed "Calculate additional properties" dropping the filter on an integer column (MW) it fills with fractions, rebinding the viewers that showed it and moving it to the end
* [core#1146](https://github.com/datagrok-ai/core/issues/1146): Campaigns: breadcrumbs are added when the pick view is ready, without an extra 300 ms
* Hit Design, PeptiHit: A new campaign's molecule column holds molblocks, so a sketched molecule is kept exactly as drawn; OK on an empty sketcher no longer calculates the properties of an empty molecule
* All apps: The landing header is u2's `AppHeader` (moved out of `@datagrok-libraries/utils`): the same icon, bullets and Learn more
* All apps: The campaigns table is u2's `CatalogTable`: every campaign is listed (no 20-row cap), the most recently saved first while no column is sorted, groups cover every campaign, start open and show their counts, a Template group is labelled by the template's name, a Date field shows its date, and Refresh reads the list from the server again; the eye menu stays open while columns are toggled
* All apps: Edit Permissions is one grid of users and groups with Can View and Can Edit, and deleting a campaign or a template asks through a confirmation that names it; a deleted campaign leaves the table in place
* All apps: The template picker offers to create the first template when there is none, says why the last template cannot be deleted, and no longer reports "Failed to delete template null" after a delete
* All apps: A link to a campaign that cannot be opened (its template is gone) lands on a working landing and says why, without logging that expected state as an error
* All apps: The landing's help is the app's README, kept as the help of the app's view (it lost to the platform's start page when the app opened from its address)
* All apps: Fixed two sessions starting a campaign from one template at the same time getting the same code: the first save takes the next free code under the campaign's lock
* All apps: Fixed a campaign that failed to open (its template or table gone) switching the app half over to it, so the campaign still open saved over the wrong campaign: what the open reads is checked first, and a failed open leaves the open campaign in place
* All apps: Fixed deleting a campaign leaving its last save, lock and V-iD links in the database, and a session letting go of another session's campaign lock
* All apps: Edit Permissions on the landing writes under the campaign's lock over the campaign as the server holds it, so another session's save made meanwhile is kept
* Hit Triage, PepTriage: Fixed a save dropping campaign.json keys it does not know, and a closed app still rewriting the address when a view of the same name became current
* Hit Design, PeptiHit: The new-campaign form is a u2 form: one label column, START disabled with the reason while a required field is empty (spaces only count as empty), and a line saying what START makes
* PeptiHit: Fixed a new campaign taking over the campaign opened before it, ignoring the Campaign Name, and failing after a view-only campaign
* Hit Triage, PepTriage: START is disabled with the reason while the file, a required field or a data-source parameter is missing; a new campaign no longer takes over the template and fields of the campaign opened before it
* Hit Triage, PepTriage: Fixed saved campaigns opening with a blank view: a campaign saved without filters no longer passes `{}` as the filter state
* All apps: Fixed the campaigns list hiding campaigns shared with All users from everyone but their author: they are listed for everyone, and editing, deleting and sharing a campaign follow its stored edit permissions (All users by default); an unreadable campaign is skipped with a warning naming it
* Hit Design, PeptiHit: The campaign view is one session owned by the view (`shared/campaign-session.ts`): the crumbs "Hit Design › <campaign>" are the platform's, shown at once (no second's delay, never on another view), the project Save is left off through the platform, and the package's ribbon group is pinned to the ribbon's end with its small icons folding into "…" at 1366 px; the progress tracker, dialogs, grid handlers and the check for other sessions' saves end with the view
* Hit Design, PeptiHit: Fixed a save over another session's save losing that session's edits to the rows both hold: a cell this session left as it read it takes the other session's value
* Hit Design, PeptiHit: Fixed a view-only campaign's refused save telling every open session of new changes: the save is refused before the lock, and the status bar says the campaign is view only
* Hit Design, PeptiHit: Fixed a save dropping campaign.json keys it does not know, and a save whose table cannot be written hanging with the campaign locked for 30 s: it lets the lock go and says why, and logs no save, so no other session is told of changes the share never got
* Hit Design, PeptiHit: Submit, Merge, Progress tracker and Modify Stages are u2 dialogs, all but the tracker with CANCEL: Submit is always there, disabled while no submit function is chosen, Enrichment Methods by friendly name, SAVE waits for the save; the tracker's lanes scroll at 1366, it closes with its campaign view; Modify Stages says before OK where the rows of a removed stage go
* Hit Design, PeptiHit: Fixed the Submit dialog's Name and File Path pencils and the download icon doing nothing: the pencil opens the value with Save and Cancel (Enter and Escape alike), a problem shows under the field, and the download saves the rows as CSV; a Molecule field's structure sits in the value column on a row as tall as the picture (it hung over the next row)
* All apps: Edit Permissions shows each user's picture or initials and each group's icon, on roomier rows across the dialog
* Hit Design, PeptiHit: Merge with table picks its source file from Datagrok Files again ("Add file from Files" beside Open file, as the platform's file input had it): a file from a share is remembered with the campaign and, with Auto-merge on, merged on every open
* Hit Design, PeptiHit: A new campaign opens with its template's layout; a layout's viewers bound to columns the table lacks (the first calculation adds them) start empty instead of failing, and such a column is saved only once it holds values
* Hit Design, PeptiHit: Calculate additional properties opens the template editor's compute box in a u2 dialog; the campaign's compute is stored as before
* Hit Design, PeptiHit: Fixed a save asked for while another campaign opened, or queued behind one, writing one campaign's rows or keys into the other: a save writes the campaign it was asked for, and an open takes the new campaign at once
* Hit Design, PeptiHit: Fixed a save over another session's save mixing up the rows of a molecule several rows hold, and reading the other session's table from our path instead of the one its campaign.json names
* Hit Design, PeptiHit: Fixed the progress tracker's cards reading another open campaign's table (a peptide's V-iD on a Hit Design card, cards missing) when the campaign's card layout named a table of the same name
* Hit Design, PeptiHit: The progress tracker shows a row whose stage the template does not have (a typo in the grid) in a lane of that stage instead of leaving it out
* Hit Design, PeptiHit: Fixed number columns showing every digit the campaign's table holds after a reload (446.9100036621094 for 446.91), in the grid and on the tracker's cards; a format the campaign chose stays
* Hit Design, PeptiHit: A save that cannot take the campaign's lock because a query fails says why, instead of "Someone else might be editing the campaign"
* Hit Design, PeptiHit: Fixed Modify Stages saving the progress tracker's Tile Viewer into the campaign's layout, a campaign whose layout has no grid no longer saving its edits, and two quick edits of a row leaving the first molecule's V-iD
* PeptiHit: Fixed the previous campaign view staying open when another campaign is opened, and the view being named after the campaign's code until its first save
* Hit Design, PeptiHit: Merge with table says why each option it keeps off is off ("Pick a source file first", "Used only when merging by V-iD", "Turn on Add unmatched rows first", …), as its tooltip and as its description for screen readers and automation
* Hit Design, PeptiHit: Merge with table, Calculate additional properties, the designer's form and the dialog for a dropped file are pinned to their campaign view, as Submit and the progress tracker are: hidden while another view is current, closed with the view
* Hit Design, PeptiHit: The template editor is a u2 page (`shared/template-editor.ts`) with its sections in the order they had, its fields and stages rows with "Add field" / "Add stage" (stages reorderable), the crumb back to the app the shell's (it closes the editor, asking first once something was typed), the help on the Templates section, and CANCEL / CREATE in a footer that stays in view; CREATE is disabled while the name or key is empty, blank or taken, with the reason, and waits for the write: a failed write, or a name or key another session took meanwhile, keeps the editor open and says why
* Hit Design, PeptiHit: The template editor's Compute is one u2 master-detail box (`shared/compute-editor.ts`): the functions on the left with their checkboxes (↑/↓ move the current one, Space checks it; a sash between the list and the settings), the descriptors as a searchable tree whose groups say how many of theirs are picked (they showed unchecked), a function's settings with "Re-run when campaign opens", dimmed with the reason while the function is not checked, and "Reorder…" to order and hide the functions; picking a descriptor turns Descriptors on
* Hit Design, PeptiHit: Fixed a cloned template keeping the name and key (now "<name> (copy)" and the key's next free number), losing its submit function, scripts, queries and layout
* Hit Triage, PepTriage: Fixed the template editor's CREATE never blocking on its own errors (an empty name wrote `templates/.json`, a taken name overwrote that template): it is disabled with the reason while the name or key is empty or taken, or PepTriage's sequence column is empty; a failed write is reported instead of "Template created successfully"; a clone keeps the template's layout
* Hit Design, PeptiHit: Fixed "Calculate additional properties" showing "Cannot read properties of undefined (reading 'label')" for every function but the descriptors while a table without a name was open
* All apps: Fixed a calculation over the whole table naming its columns after a hidden column ("Acceptor (~Molecule.canonicalSmiles)", beside the "Acceptor (Molecule)" a row's own calculation makes) and adding "(2)" copies when run again: the functions get the campaign's columns, the molecule column holding the canonical SMILES under its own name, and what they add or change fills the columns of that name; a hidden `~` column is never taken for the molecule column, never written to campaign.json, and dropped when a campaign's table holds one
* Hit Design, PeptiHit: Fixed a new campaign whose template layout has a Filters viewer showing "Error adding filter" for every filter on a column the campaign calculates later, and losing those filters: their columns are added empty with the layout's other columns, and a calculation's counts fill them without dropping the filters
* Hit Design, PeptiHit: Which columns a layout's viewers are bound to is the platform's answer (`DG.ViewLayout.boundColumnNames`) instead of the package's own reading of the layout, so a trellis plot's inner viewers are covered too
* Hit Design, PeptiHit: Fixed "Add new row" in the ribbon logging "No such column: then" on every press
* Hit Triage, PepTriage: Template, Campaign Name and the data source are one label column again, and START sits at the form's end, as before the u2 landing
* Hit Triage, PepTriage: START over a data source that answers no rows creates no campaign and says so ("<source> returned no rows, so no campaign was created"); it made an empty campaign silently
* All apps: Boolean campaign fields read true / false in the campaigns table again; a group by a Date field is labelled as its cells write the date (2026/6/29); "Apply template layout" has its platform icon again
* Hit Design, PeptiHit: A Stage cell is a dropdown of the template's stages, set when the campaign opens and when Modify Stages changes them, never saved with the campaign; a row in a stage the template lacks keeps it (in its dropdown too) until another is picked; the template editor's stages stay free text
* Hit Design, PeptiHit: Closing the progress tracker never puts a move back: its one action is CLOSE, and CLOSE, ✕, Escape and closing the campaign view all keep the moves and save them
* Hit Design, PeptiHit: Fixed a campaign view whose first read of the last save failed never noticing other sessions' saves: the check every 3 s starts all the same; a view closed while it read starts none
* Hit Design, PeptiHit: Merge with table no longer offers hidden columns (Chem's `~<column>.canonicalSmiles` cache, which a table saved as CSV keeps): the molecule and V-iD column choices list the columns the platform's pickers do, and pick their defaults among them; no merge (the dialog, a dropped table, Auto-merge on open) copies an incoming table's hidden columns into the campaign any more
* All apps: Hidden `~` columns are no longer shown to the AI assistant (a function's table result lists only the columns the platform's pickers do), and Hit Triage's campaign.json no longer names them in its column types; a campaign.json an older version wrote with them opens as before

## 1.10.9 (2026-10-15)

* Hit Design, PeptiHit: Fixed the "Campaing Name" typo in the new-campaign form
* Hit Design, PeptiHit: Fixed column tags (such as `.choices`) and calculated columns' formulas being lost when a campaign was saved and reopened; a reopened calculated column follows its source again
* AI: Exposed AI view functions on all four Hit apps — the landing views list/open/create/delete campaigns and templates, and the campaign table views get a hit-specific vocabulary (state, molecules, stages, status, save, submit) plus a search/call bridge to the standard table commands.

## 1.10.6 (2026-07-16)

* Hit Design: Added a "Re-run when campaign opens" option per compute function/script/query, available both in a campaign's compute dialog and at template creation. When enabled, the function is recomputed over the whole molecule column every time the campaign is opened — useful for database lookups that resolve after a molecule was first drawn.
* Hit Design: Recomputing a function/script/query now overwrites the existing result column in place instead of appending a duplicate `<name> (2)` column.

## 1.10.2 (2026-05-15)

* Auto merging of data file if configured in campaign settings.

## 1.10.0 (2026-04-09)

### Features

* Introduced PepTriage — a peptide-focused workflow alongside Hit Triage and Hit Design
* Hit Triage:
  * Added campaign grouping by template, status, author, last modified user, or custom fields
  * Added sortable campaigns table columns
  * Added column visibility editor for the campaigns table
  * Migrated the compute function engine to a shared implementation now used together with Chem MPO

### Bug Fixes

* Browse panel compatibility restored

## 1.9.1 (2025-02-03)

* HD: Parallel running functions, fix table joining
* Hit Design: Custom comparer for vid
* GROK-19523: AppTreeBrowser decorator: Fixed behavior, so it finds the app properly, fixed usages

## 1.9.0 (2025-12-17)

* Global Vid system and molecule tracking across multiple campaigns.

## 1.8.0 (2025-12-09)

* Improved collaboration support
* Automatic sygnaling of changes in campaigns between multiple users.

## 1.7.7 (2025-09-19)

* Sortable campaigns table columns via doubleclick.
* Add campaingPropts to both visible columns and groupby options.

## 1.7.5 (2025-09-03)

* Add Ability to hide or show custom columns in campaigns table.
* Improved styling of icons in the info view.

## 1.7.3 (2025-05-01)

* Fix descriptors tree not scrolling

## 1.7.2 (2025-04-09)

* Ability to reorder/hide functions in new template or functions dialog.

## 1.7.0 (2025-04-01)

* Fix not saving campaigns

## 1.5.0 (2025-03-01)

* Add caching to loaded campaigns.
* Parallel loading of campaigns.
* Lazy loading of function dialogs.
* Fixes to url logic.
* Working with new browse version.

## 1.4.4 (2025-02-07)

* Fix permission check for admin members
* Add Default file storage for campaign files

## 1.4.3 (2025-02-06)

* Fix campaign auto-naming when creating new campaigns

## 1.4.2 (2025-02-03)

Correctly set semantic types for all columns in joining dataframes

## 1.4.1 (2025-01-31)

* Support current molecule sending to designer functions
* Support multiple molecule columns in joining dataframe

## 1.4.0 (2025-01-16)

* Support adding tables to hit design campaigns
* Support molecular generation functions such as reinvent

## 1.3.11 (2024-12-11)

* Add loader to campaigns app view

## 1.3.10 (2024-12-04)

* Bump dependencies versions

## 1.3.9 (2024-11-29)

* Enable working with molblocks without overwriting original data with smiles.
* Template layouts, data sync.

## 1.3.8 (2024-11-20)

* Add Campaign name setting and changing.

## 1.3.7 (2024-11-17)

* Add Package level management of default sharing permissions.

## 1.3.6 (2024-11-08)

* Correct grouping of campaigns table, Added author/last modified user info.

## 1.3.5 (2024-11-07)

Correct readme links

## 1.3.4 (2024-11-07)

* Groupping of campaigns
* Status change of campaigns

## 1.3.3 (2024-10-28)

* Modification of HD Stages within campaign.

## 1.3.2 (2024-10-17)

* Fix HT Campaigns and templates

## 1.3.1 (2024-10-10)

* Remove project save button from the ribbon panel

## 1.2.2 (2024-08-26)

* Fix progress tracker layout loading

## 1.2.1 (2024-08-23)

* Fix issues related to progress tracker
* Fix layout grid related issues

## 1.2.0 (2024-07-30)

Permission management and sharing of campaigns. Deleting templates. optimization of function loading

## 1.1.16 (2024-07-29)

Fix re-running calculations on the first index.

## 1.1.15 (2024-07-29)

Fixed duplicate detection on fingerprints. added rerunning calculations

## 1.1.14 (2024-07-26)

Added detection of duplicate molecules in campaign and vid coloring of duplicates.

## 1.1.13 (2024-07-23)

Fix chem descriptors breaking loop.

## 1.1.12 (2024-07-18)

Fix broken demo campaign file.

## 1.1.11 (2024-06-03)

* Add ability for collaboration on same campaign.
* Fixed bugs with tiles viewer forms crashing in case of removing columns.

## 1.1.10 (2024-04-09)

Fixed bugs with tiles viewer forms and functions loading.

## 1.1.9 (2024-03-27)

Bug and style fixes.

## 1.1.8 (2024-02-22)

Add python script support.

## 1.1.6 (2024-01-26)

* Fixes to adding new functions to the campaign.
* Conversion to canonical smiles format before calling compute functions.

## 1.1.5 (2024-01-25)

* Ability to use queries with `HitTriageFunction` tag directly as compute functions.
* Fixes and improvements to performance.

## 1.1.4 (2024-01-17)

* Add ability to use queries directly as source functions.
* Fixed bugs with hit triage, including incorrect saving and addition of calculated properties.

## 1.1.3 (2023-12-27)

Fix incompatibility with old api version

## 1.1.2 (2023-12-26)

Hit Triage application improvements

## 1.1.1 (2023-11-20)

Hit Design application improvements

### Bug Fixes

Hit Design:

* New campaign loading incorrectly after closing old one.
* Correct Validation of template name, key and campaign fields.
* Last Campaign field / Stage not being added to template correctly.
* Campaigns not being saved correctly after calculations.
* Progress tracker still being created if no stages defined in the template.
* Fixed date format in campaign summary form.

### Features

Hit Design:

* Progress tracker view now is available through the button on the ribbon panel and not added automatically.
* Ability to add/remove calculated functions to saved campaigns (Using the `🔧` icon in the ribbon panel before
  the `Progress tracker` button).
* Ability to add new rows from the ribbon panel using `+` icon.
* Campaigns fields can now support molecule inputs.
* Template cloning allows to create a new template based on an existing one.
* Ability to change the location of saved dataframe (could be users own connection) throught the `Submit` button.
* Ability to delete campaigns from `Campaigns table`.
* Ability to add/duplicate molecules and edit them from the context menu of any molecule cell.
