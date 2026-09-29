# Tutorials changelog

## v.next

* Scripting: Fixed the step counter (11 declared for 10 steps)
* Dashboards: Fixed: the connection step named the menu item "Add connection..." — it is "New connection..."; the Dashboards step looked its tree row up before "Close all" rebuilt the tree and could miss the click
* Data Connectors: Fixed: the counter was 12 for 11 steps; the connection step named the menu item "Add connection..." — it is "New connection..."
* Sticky Meta: Fixed: the Types and Schemas steps looked their tree nodes up before the Sticky Meta group had loaded and listened on the caption only, so they could never complete; the entity checkbox was looked up before the dialog's rows existed
* Tutorials: A step that waits for a click sees it even when the control stops the click's propagation (a Browse tree row); the buttons of the dialog tours are named like the other tours'
* Differential equations: Fixed: the counter was 15 for 14 steps; the Lotka-Volterra model was taken as the third library card; "Set Delta to 0.1" watched the gamma input (the added eta parameter shifted the form) and could never complete — inputs are found by caption now
* Activity Cliffs: Fixed: the Show only cliffs step completed on a click on any switch of the page, and its hint was captured when the step began; it now waits for the plot's own filter
* R-Groups Analysis: Fixed: the counter never reached its last step (14 declared for 15 actions); the gear step completed on a click on any gear of the page; the Category and Value steps watched the context panel's text and now read the trellis's inner viewer; the hover step threw when the Distributions pane was not built yet; a stale Jupyter prerequisite is gone
* Similarity and Diversity Search: Fixed: the counter never reached its last step; the gear step completed on a click on any gear of the page; the Follow Current Row and Molecule Properties steps waited for texts in the context panel and now read the viewers' properties; the Edit hint looked for a class the viewer does not carry, so it showed nothing
* Substructure Search and Filtering: Fixed: the counter never reached its last step; switching the card off completed on any class change of the card and now waits for the card to say it is off; three hints were captured when their step began; typos "flter" and "filers" in the steps
* Multivariate Analysis: Fixed: the counter never reached its last step; the viewer tour took the charts two seconds after RUN by position (a slower first run left it short or wrong) and now waits for them by title; the tour's texts were one off from the charts, since the Scores text was split in two
* Filters: Fixed: the row-count and bin-selection steps completed on any selection event; the indicator and reset hints were captured when their step began and could point at nothing
* Grid Customization: Fixed: the row-selection step asked for a plain click on a row number, which makes the row current and selects nothing; it asks for a Shift+click now. The new-row and remove-rows hints follow their controls instead of the nodes captured when the step began
* Viewers: Fixed: the counter never reached its last step (20 declared for 20 actions); three selection steps completed on any selection event, so the tail of one gesture could tick the next; the gallery's search box was captured when its step began
* Calculated Columns: Fixed: the expression step was declared twice — the first copy was skipped or shown depending on the dialog's layout, shifting every step number after it — and the Edit step pointed at a button called "Edit" while the pane says "Edit in dialog" (no highlight); both fixed
* GROK-20419: Embedded Viewers: The last step pointed at a gear icon no element carries, so it showed no highlight; it now points at the inner plot's Color selector in the Trellis strip and says so, and the Trellis plot is found even when it is added after the menu click
* Added a BDD project (`bdd/`): the tutorials walked from their cards to the congratulations, each step claimed both as ticked and as done on the platform
* Demo app: A demo started from the Browse tree fires the `demo-loaded` custom event with the demo's `path` once it has run and its view is named
* Tutorial cards carry `role="button"`, `aria-label` (the tutorial's name) and `data-status` (`done` / `not done`)
* Scatter Plot: Fixed: "Select points" and "Deselect points" completed on the same selection event, so the trailing event of a Shift-drag could tick "Deselect points" with no Escape; they now wait for a selection and for an empty one
* R-Groups Analysis: The Distributions pane step resolves the pane when it appears instead of when the step is built — it exists only while a RowGroup is the current object
* GROK-20416: Steps that point at a rebuilt element (the Add viewer ribbon icon, the aggregation tag) now resolve it per tick instead of caching it, so both the click listener and the highlight follow the new node; the icon is found by its `aria-label`, not by its position in the ribbon
* Sticky Meta: Fixed the tutorial dying at "Save schema" — step 11 looked up the schema section and the property input with strings that disagreed in case with what the steps tell the learner to type, and the miss threw one line later
* Activity Cliffs: Fixed two steps that could never complete — the scatter plot is added before the step subscribes to `onViewerAdded`, and the tooltip check needed a second mouse move; both now read the state that is already there
* Compute: Sensitivity Analysis and Parameter Optimization address form inputs by the caption the model declares instead of by child index, and the tour wizard buttons carry `name` / `aria-label` even when the caller passes an empty label
* Dashboard: Dropped two 1.5 s "UI generation delay" sleeps — the engine now waits for the control it is about to drive
* Tutorials: Dashboards tutorial follows the dashboard it saved by id — on a shared server its last steps accepted any project named "Coffee sales dashboard", so a project from someone else's run could be opened instead and fail with a connection-permission error

## 1.11.4 (2026-09-11)

* Tutorials: Links to the Tutorials app and to individual tutorials can now be shared and opened directly
* Demo app: Correlation Plot is offered as a demo again - its `demoPath` was missing a space, so consumers that rebuild the path from the hierarchy never matched it
* Demo app: The second viewer pane in a viewer demo is titled after the viewer it holds, not always "histogram"
* Tutorials: A track with no tutorials renders a 0% progress bar instead of `NaN%`
* Fixed the package build failing on `TS2416` — the u2 `Component` base introduced `run(fn)`, which every widget now inherits, so `Tutorial`'s own `run()` no longer matched; it is now `start()`
* GROK-20602: BREAKING regen — `grok api` codegen v2 for the Northwind demo schema: datetime fields are dayjs, typed expand/transaction surface, lazy db.ts clients
* Demo app: Added a Domain Databases demo (Data Access | Domain Databases) — ships the classic Northwind schema and data as a plugin-declared domain database (databases/northwind) and walks through browsing, security, audit history, and the JS API
* GROK-19337: Peptides-SAR tutorial: Fixed callouts rendering body text as bold headings (`paragraphs()` now emits a bold title only when one is explicitly passed), removed spurious mid-sentence `<br>` line breaks, and narrowed the action panel (~18% narrower dock ratio)
* Demo app: Always clear the "Updating..." overlay when a demo finishes — moved `setUpdateIndicator(false)` into a `finally` so a failing demo no longer leaves the indicator covering the home page
* Demo app: Fixed clipped/misplaced "Updating..." overlay — place it on the tab-content of the active view instead of `tab-content[0]`, so it stays bounded (tab-content is `position: relative`) and lands on the pane that's actually being replaced rather than a sibling (e.g. the home page) when the dock is split
* Demo app: Made the "Updating..." overlay clearly readable during demo loading — scoped CSS (`.demo-app-loading > .d4-update-shadow`) gives the overlay a solid translucent white background with a subtle backdrop blur and renders the label/loader at full opacity, so widget text underneath no longer bleeds through
* Demo app: Fixed view leak when reopening a multi-table demo (e.g. Table Linking) — newly opened views are tagged via `view.temp['demoApp']` and all of them are closed on the next demo switch, instead of only `currentView`
* GROK-15442: Tutorials: Added a Sunburst viewer step to the Viewers tutorial
* Tutorials: Fixed the Radar step highlight in the Viewers tutorial landing on the hidden "Recently used" card instead of the visible gallery card
* Tutorials: Updated the Diff Studio tutorial with respect to the main app modifications
* GROK-20145: Demo | Form viewer: Error on URL direct opening
* GROK-20408: Tutorials: Fixed `waitForElementClick` dropping the first 500ms of clicks, hanging forever on a not-yet-rendered element, and leaking listeners with no cancellation — replaced by `elementClick` (immediate listener, polling getter with timeout, Observable so `firstEvent` cancels on close)
* GROK-20648: Multivariate Analysis tutorial: Improved the step and viewer descriptions, and fixed the fourth viewer hint describing Regression Coefficients while pointing at Variable Importance
* GROK-20424: Tutorials: Fixed the Calculated Columns tutorial hanging on the "Enter the expression" step — the MutationObserver watched `.cm-line`, which is null before CodeMirror renders, throwing `observe: parameter 1 is not of type 'Node'`; now watches the always-present dialog root subtree

## 1.11.2 (2026-04-10)

* Added demo image for Retrosynthesis

## 1.11.1 (2026-03-23)

### Bug Fixes

* GROK-17970: Demo app no longer closes user-opened content when running demos
* GROK-19746: Demo app browse-toolbox no longer flickers when switching between demos

## 1.11.0 (2026-03-13)

## 1.10.6 (2025-12-15)

* Tutorial: Updated the scientific computing tutorials

## 1.10.5 (2025-10-31)

* Tutorial: Updated the parameter optimization tutorial

## 1.10.4 (2025-10-30)

* Fitting: New layout related update

## 1.10.3 (2025-10-08)

* GROK-18139: Activity dashboard: Retired Learn widget and added it to the activity dashboard instead

## 1.10.2 (2025-10-07)

* Demo app: Fixed the issue that names with `-` in name sometimes can't open

## 1.10.1 (2025-10-07)

* GROK-18767: Fixed opening Demo in learn widget opens empty toolbox
* Fixed Demo > Cheminformatics is opened even for old users

## 1.10.0 (2025-09-30)

* GROK-18139: Added demo app hierarchy as function in Tutorials
* GROK-18718: Expanded Demo > Cheminformatics by default for new users

## 1.9.6 (2025-09-25)

## 1.9.5 (2025-09-19)

* GROK-18879: Fixed Data aggregation tutorial
* GROK-17845: Fixed each click opens additional app menu in Tutorials widget
* GROK-18791: Fixed toolbox isn't shown in Scatterplot tutorial
* Added Peptides tutorial

## 1.9.3 (2025-09-17)

* Dashboards tutorial: Fixed UI inconsistency in names
* GROK-17832: Fixes in sync with new Browse UX

## 1.9.2 (2025-09-15)

Fixed the Parameters Optimization tutorial

## 1.8.5 (2025-04-25)

Updated UI in the Multivariate Analysis tutorial

## 1.8.4 (2025-04-22)

Fixed the behavior of hints in the Parameters Optimization tutorial

## 1.8.3 (2025-04-11)

Updated the Parameters Optimization tutorial

## 1.8.2 (2025-04-04)

Updated Demos

## 1.8.1 (2025-03-31)

Updated the Scientific Computing tutorials with respect to new Model Hub UI

## 1.8.0 (2025-03-29)

Updated demo features and the Scientific Computing tutorials with respect to new Browse UX

## 1.7.3 (2025-01-10)

Added scientific computing tutorials:

* Differential equations
* Sensitivity analysis
* Parameters optimization

## 1.3.27 (2024-07-19)

### Bug fixes

* Fixed Cheminformatics tutorials

## 1.3.26 (2024-06-04)

### Features

* GROK-15447: Tutorials: Add radar as a part of viewers tutorial

### Bug Fixes

* GROK-15228: Tutorials | Dashboards: Not running
* GROK-15620: Tutorials | Scripting: Cannot read properties of undefined (reading 'content')
* Fixed connection to starbucks

## 1.3.25 (2024-04-16)

Updated the Multivariate analysis tutorial.

## 1.3.24 (2024-04-15)

Fixed demo cards' images for Docking and Bio Docking conformations

## 1.3.23 (2024-04-11)

Fixed a couple of bugs in Cheminformatics tutorials

## 1.3.22 (2023-12-21)

### Features

* Added card for Diff Studio demo app.

## 1.3.21 (2023-12-07)

### Bug Fixes

* Fixed dock size.

## 1.3.20 (2023-12-07)

### Features

* Added card for SimPKPD demo app.

## 1.3.19 (2023-12-04)

### Bug Fixes

* Fixed demo app not working.

## 1.3.18 (2023-12-04)

### Features

* Improved Chem tutorials: Similarity and Diversity Search.

### Bug Fixes

* Fixed demo app not working.

## 1.3.17 (2023-11-09)

### Features

* Improved Chem tutorials: R-Group Analysis.

## 1.3.16 (2023-11-07)

### Features

* Improved Chem tutorials: Activity Cliffs, R-Group Analysis, Similarity and Diversity Search, Substructure Search and Filtering.

## 1.3.15 (2023-11-02)

### Bug Fixes

* Improved Chem tutorials: Activity Cliffs, R-Group Analysis, Similarity and Diversity Search, Substructure Search and Filtering.
* Improved Scatterplot and Dashboard tutorials.
* Changed tutorials order in the EDA track.
* Fixed links in tutorials.

## 1.3.14 (2023-10-23)

### Features

* Added tests for all demo functions.
* Added PK-PD modeling demo.
* Added Chem tutorials: Activity Cliffs, R-Group Analysis, Similarity and Diversity Search, Substructure Search and Filtering.

### Bug Fixes

* Harmonized menu names.
* Moved Bioreactor demo to Compute section.
* Fixed multivariate analysis tutorial.

## 1.3.13 (2023-08-28)

### Bug fixes

* [GROK-13826](https://reddata.atlassian.net/browse/GROK-13826): Tutorials:
  * Temporarily commented the calculated columns tutorial.
  * Temporarily commented the predictive modelling tutorial.
* [GROK-13824](https://reddata.atlassian.net/browse/GROK-13824): Tutorials:
  * Fixed doc links.
  * Replaced PostgreSQL by Postgres.

## 1.3.12 (2023-08-24)

### Bug fixes

* [GROK-13357](https://reddata.atlassian.net/browse/GROK-13357): Tutorials: Fixed indicator on current view change.

## 1.3.11 (2023-08-24)

## 1.3.10 (2023-08-07)

## 1.3.9 (2023-08-07)

### Features

* Updated steps description in Multivariate Analysis tutorial

### Bug fixes

* Recognition steps in Multivariate Analysis tutorial

## 1.3.8 (2023-07-21)

### Features

* Added [Bioreactors demo](https://public.datagrok.ai/apps/Tutorials/Demo/Bioreactors).
* Added demo for [heatmap](https://public.datagrok.ai/apps/Tutorials/Demo/Visualization/General/Heatmap) with dendrogram.
* Added demo for [Chem](https://public.datagrok.ai/apps/Tutorials/Demo/Cheminformatics).
