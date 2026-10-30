# Compute2 changelog

## v.next

* Help: Fixed the help panel staying empty for every later step after one help file failed to load
* Runs: Added the Save run dialog: saves to My runs or a space with a name, description and tags; names are unique in a space, where a taken run name asks to replace it (the space link moves to the new run) and a name taken by another item is refused; its spaces open like folders, and clicking a run of a space offers to overwrite it. My runs lists all of the user's runs, including those saved to spaces but not the step runs saved as part of a workflow, and names can repeat there
* Runs: The Runs panel, Model Hub's Tools > Browse Runs... and the Compare Runs run picker browse runs of all models or one in My runs and spaces like a file browser: subspaces as folders with `..` and a path, tree counts of the runs in a space and all its subspaces, linked subspaces under each parent, and filters that search a space and its subspaces and show where each run is; they offer tag, author and favorites filters, platform favorites for runs and spaces (the Save dialog too), a tree filter that shows only favorite spaces and the spaces above them, Edit details, Add to space for one or several runs and Delete with a confirmation for one or several runs (a run that fails to delete is named and stays selected, and the others are still deleted); they replace the history panel
* Runs: Removed Share to workspace and the `sharingMethod` setting
* Custom function views: Save opens the Save run dialog and History the Open run dialog instead of a docked history panel, so the whole view is the custom content; saved runs are listed under the custom view's function in the run browsers and open in that view again
* Run compare: Several comparisons of the same runs open as tabs in one view; each tab has its own selected values, Hide equal values, chart and export, while index columns, mappings, matched columns and Hide outdated outputs are shared
* Workflows: A workflow that was never saved shows its root `title` link in the view name
* Workflows: Fixed Run ready steps running a step before its debounced or async validators checked its new inputs, or before an async link set them; the workflow stays busy until they finish, and validator debounce applies to edits only
* Workflows: Added `RunWorkflow` to run a workflow without UI from an initial configuration, as `StartWorkflow` takes it; it returns every step's values, status and validations, and a summary. `SaveWorkflowRun` saves the result to history
* RFV: In workflows, inputs with the `disabled` or `restricted` restriction are hidden behind a Show linked inputs row and show a link icon while the value is the one the workflow set; `meta.inputToggles` picks the inputs (by category or name, linked or not) and the row captions, the `keepVisible` input annotation exempts an input, and inputs with errors, warnings or an inconsistent value are always shown
* RFV: In workflows, `meta.inputsHidden` hides every input behind the Show linked inputs row instead of collapsing the form; inputs with errors, warnings or an inconsistent value stay visible. The standalone RFV ignores it and `meta.inputToggles`
* RFV: In workflows, a dataframe input with the `mutateOnly` annotation can't be replaced under the `restricted` or `info` restriction, only edited in its viewers
* RFV: The `standalone` prop replaces `localValidation`, `skipInit`, `skipLogic`, `showRunButton` and `keepExportsVisible`; autorun is read from `runOnInput`
* RFV: Malformed JSON in `outputCategoryGroups` or an output's `hidden` annotation logs a warning instead of breaking the view
* Workflows: Nested pipelines have two history loads: Load completed run loads an editable copy, Load completed run (read-only) a locked one; read-only steps and pipelines show a lock on the left of their tree row; both have a branch icon at the end of the row whose menu opens the original run; steps name the source run in the tooltip and menu; a loaded pipeline also has a line in its view to open the original run in a new view
* Workflows: An editable step or pipeline based on a completed run compares its current state with that run (Compare with original in the branch icon menu or the pipeline view), or adds it to a comparison (Add to comparison...). The compare icon next to Save adds the current state of a workflow or a standalone model to a comparison. These and Compare in the run browsers open one dialog: each run has an editable name, and when comparisons are open you pick one or a new one (Compare and Compare with original start on a new one, the rest on the most recent). The current state is a snapshot that is not saved to history; the dialog warns about outdated outputs, and the compared model joins the comparison's model list
* Run compare: Snapshots of the current state are marked unsaved and ask before removal, and a save icon saves one to history; runs with outdated outputs are marked, and Hide outdated outputs leaves those outputs out
* Workflows: Steps with history have a lock ribbon icon next to history that loads a completed run read-only
* Workflows: Fixed OpenWorkflowRun opening the default input dialog for runs of functions without a Compute2 editor, such as workflow steps; they open in the Rich Function View, nested pipeline runs in the Tree Wizard
* Workflows: Read-only subtrees can no longer be duplicated, and duplicating a pipeline that contains one keeps it read-only
* Workflows: Actions hidden by a node meta link (`setActionVisible`) are not shown
* Run compare: Selected runs and workflow steps load in parallel; a run that fails to load no longer stops the others, and runs with steps that could not be loaded show a warning
* Run compare: Fixed runs with the same name merging in charts: a repeated name gets a ` (2)`, ` (3)`... suffix
* Run compare: Fixed an added run not taking the index picked for the same function's runs, which cleared the merged index row and left the run out of the charts; split and axis mode picks carry over too
* Run compare: Empty index cells no longer show up as an `undefined` category on a categorical index
* Run compare: Fixed a blank chart when the first checked value disappears while others stay checked; a charted value stays listed when Hide equal values would hide it
* Run compare: Chart and export columns carry their units in the name (`Biomass (g/L)`, `Time (h)`), and the Compare list shows each value's units
* Run compare: A run that has the value in other units shows `units differ (mg/L)` instead of `no similar data`
* Run compare: Clicking the selected value again clears the selection; the expand chevron takes its whole cell, so clicks near it no longer change the selection
* Run compare: Unit mappings (`1 g/L = 1000 mg/L`) in the Mappings dialog let values in other units join a row; they are multiplied into the first run's units when the chart is built
* Run compare: Units compare case-sensitively (`MPa` ≠ `mPa`), folding only spacing, micro signs, superscripts, `l`/`L` and time unit spellings; the time units picker reads the same spellings, so `H` (henry) is no longer taken for hours
* Run compare: The Multiple values toggle sits next to the value filter, and the name mappings icon got a label
* Run compare: Charts of several values show every run any of them has, and the partial badge no longer depends on which value was checked first
* Run compare: Fixed repeated workflow steps (e.g. one per day) overwriting each other's tables and showing as identical rows; each occurrence is now its own row labeled with its step, and runs pair by step
* Run compare: Several dataframe columns from different workflow steps chart together, like scalars; each step of a run is its own line
* Run compare: Several columns in the same units chart in one panel on a shared y axis; the Shared y axis switch brings back the stacked panels
* Run compare: The Mappings dialog lists values that are in no compare row, with shortcuts to map their name or units
* Run compare: Raw tables without units no longer raise the missing-units warning; an info icon says their values are taken in the row's units
* RFV: Dataframe outputs annotated `hidden: true` are hidden, like scalar and graphics outputs
* TreeWizard: Step tooltips no longer name inputs that were only changed among those with warnings
* TreeWizard: Removing a step or workflow from the tree asks for confirmation
* Replaced hardcoded colours with the platform design tokens
* TreeWizard: The selected step shows its actions without hovering, and tree icons are labelled for assistive tools
* TreeWizard: Each step status has its own icon, so statuses no longer differ only by colour
* RFV: Added the Feedback menu (Report a bug, Request a feature) to the standalone view
* Inspector: Links tab shows each link as written, empty aliases, targets dropped by `$linked` and action visibility, and no longer mixes up the spec of links that share an id across steps; Config tab shows the config as written with resolved step io; links generated from annotations are shown in the Links and Log tabs only with the Annotation links toggle, and missing-value checks are never shown
* TreeWizard: Added duplicating a dynamic workflow item (step or nested workflow) from the tree; the copy keeps input and state values, and the outputs of a run step, and does not rerun init hooks
* TreeWizard: The workflow Excel export adds a `000_summary.xlsx` with each step's status, run error, messages and inconsistent inputs, and each workflow's validations and step counts; custom exports get the same data from `getExportSummary` and `reportSummaryExcel`, and the step callback gets the `status`
* TreeWizard: Step forms only render and bind (`skipLogic`): the workflow owns values, choices, validation, visibility and enabled state, and the form no longer overrides them
* TreeWizard: Opt-in compact view for single-step workflows (`compactView: true`; no tree or navigation; Save run, the Runs panel and export on the step ribbon, the Runs panel always available)
* RFV: Default Excel export includes validation and consistency
* Run compare: Fixed a fuzzy-named value taking a run's slot from an exact match declared after it
* RFV: Graphics outputs (R/Python plots) render as images in their own panel instead of base64 text, and export to Excel
* TreeWizard: Fixed validation advice actions running without their `additionalParams`
* TreeWizard: The failed step tooltip shows the run error

## 1.6.0 (2026-09-23)

* TreeWizard: Fixed saving a step to history overwriting the previously saved run (a loaded or re-saved step now gets a new id; saved runs are immutable)
* Ribbon panels and menus migrated to the injected per-view service from webcomponents-vue; removed view prop drilling, busy-guarded ribbon actions use the built-in disabled API
* TreeWizard: Fixed the subtree save icon saving the currently selected step instead of the clicked one
* TreeWizard: Stopped re-creating the pipeline FuncCall on every re-render (help panel churn and repeated context help fetches)
* RFV: "Copy link with inputs" reads entity ids off the current table/file values instead of ids recorded at load time, so a replaced input no longer resurfaces through the copied link
* Failed run loads, saves, shares, viewers hooks, and exports now show a shell error instead of failing silently
* TreeWizard: Fixed step deep links encoding only the first path segment separator (broke depth-3+ links pasted as text)
* RFV: Fixed dataframe param subscriptions accumulating on every tab rebuild (leaked across runs, step switches, and unmount)
* Excel exports detach the ad-hoc viewers they create (one leaked set per export click before)
* History: Fixed the Version column always showing; edit/delete dialogs cleanly abort on cancel and Esc; failed run loads on row selection show a shell error
* RFV: Fixed the skipInit prop default typo; report/export/SA/fitting failures show a shell error; scalar table rows align with and without validation icons; restored the scalar row separator
* Run compare: Fixed int64 columns corrupting the chart outside the elapsed-time axis; unchecking the anchored value moves the anchor to a charted one; export survives the comparison changing behind the dialog; the add-table input is released on close
* Added optional run sharing (workflow, RFV model run, and per-step when the step has enableHistory) via the sharingMethod setting: workspaces (links the run into an existing space, saving it to history first when needed)
* Saved Compute2 runs (e.g. linked into a space) now render and open through their own editor from Browse instead of the generic funccall view
* Added OpenWorkflowRun for opening a saved workflow or function run by its FuncCall id
* Fixed view pinning failing depending on the platform's js-api version (pin vs isPinned)
* Run compare: Fixed Save & Share failing to save and producing projects that crash on open

* RFV: Added a default Grid viewer for dataframe outputs without a viewer annotation (e.g. query results)
* RFV: Added a "Save as project" export (visible tables and viewers; data sync when inputs are literals)
* RFV: URL parameters as input overrides with auto-run, and a "Copy link with inputs" export

## 1.5.10 (2026-08-07)

* Fixed context help not showing for functions annotated with `readme` instead of `help`
* Added a general-purpose run comparison tool (Model Hub > Tools > Compare Runs)

## 1.5.8 (2026-07-07)

- Model Hub: add the `roleOnlyModelFilter` package setting to use the faster role-only model filter (off by default; the legacy tag-or-role filter is unchanged)
- Pick up compute-utils 1.46.6

## 1.5.7 (2026-07-02)

- TreeWizard: fix ribbon flicker on step runs and link updates by keeping the Actions, Run-ready, Save and Confirm controls visible and guarding their clicks instead of hiding them on transient locks

## 1.5.6 (2026-06-30)

- Use the G7 float format for form inputs on platforms with js-api >= 1.27.7; rendering and export keep `#0.###` and honor per-field format overrides
- Skip empty dataframes when exporting to Excel instead of producing blank tabs
- Enable help files for action steps via `nqName`
- Pick up compute-utils 1.46.5

## 1.5.5 (2026-06-25)

- RFV: don't duplicate the Diff Studio lookup "Default" scenario in Sensitivity Analysis and Fitting (pass `disableLookupDefault`)
- Pick up compute-utils 1.46.4

## 1.5.4 (2026-06-16)

- RFV: Help panel is dockable via `dockSpawnConfig['Help']`, so models can anchor it instead of it landing on the first chart
- RFV: multithreaded fitting for Diff Studio models
- RFV: propagate Diff Studio mapped (lookup) input to Sensitivity Analysis and Fitting views
- Update the standalone RFV URL with the run id on save
- Fix custom export running the model instead of the declared function
- GROK-19187: Fix back navigation duplicating the Model Hub tab
- GROK-19188: Fix models not having a proper URL
- Track the RFV focused tab separately for inputs and outputs
- Use type-only imports/exports to silence rspack warnings
- Fix custom view test code
- Pick up compute-utils 1.46.3

## 1.5.3 (2026-06-05)

- Fix navigation tree desync after a drag reorder: clear stranded `dragNode` state in the after-drop handler and treat the driver tree as the single source of truth
- Repair tree selection when the chosen step is removed: climb to the nearest surviving ancestor instead of leaving a dead uuid selected
- Stream pipeline-validator results to the navigation tree so the validator icon updates as soon as the validator resolves
- Pick up compute-utils 1.46.2

## 1.5.2 (2026-06-04)

- Default focus to the Inputs tab when the form is shown as a tab
- Pick up compute-utils file-input and `pipelineValidator` changes

## 1.5.1 (2026-06-04)

- Fix float-display test for the `#0.###` default mask

## 1.5.0 (2026-06-03)

### Features

- Per-step funccall history in the tree wizard
- `inputsHidden` with inline toggle; `formAsTab` replaces `formOnly`
- Inspector: parsed link IOs rendered structurally, clickable link badges on error-log items, actions vs. links distinguished, unified inputs/outputs shape across tabs
- Honor RTD action visibility in TreeWizard; render `nodeMeta` body on the pipeline action page
- Output category groups for scalar outputs; hide empty outputs with tab persistence
- Render boolean scalars in the scalars table

### Build

- Default build now uses rspack (webpack kept as `build-webpack`); added rspack + tsgo option, Babel for `.tsx`

### UI

- Default float mask `#0.###` across outputs and grids
- Hide ribbon actions and input viewer panels when UI is blocked / `meta.hidden` is set

### Bug fixes

- Fixed latent type errors surfaced by tsgo
- Refresh RFV scalars after buffer-collapsed re-runs; focus/persist Inputs tab when `formAsTab` is on

## 1.4.1 (2026-04-24)

### Features

- Compositor-thread loading overlay
- Action step type with unified Back/Next pipeline navigation
- Inspector improvements: FilterDropdown component, unified log filters, improved log and path formatting
- Inspector debug tool available to all users
- Stress test pipeline with onInit and chained steps
- E2e test suite for apps, RFV editors, and pipelines

### Performance

- Buffer RxJS-to-Vue reactive updates during global lock
- Suppress Vue rerenders during global lock
- Webpack filesystem cache and babel-loader caching for faster rebuilds

### UI

- Left-aligned navigation buttons in RFV and PipelineView bottom bar

## 1.4.0 (2026-03-19)

- 1.27 platform support
- More export improvements

## 1.3.22 (2026-02-19)

- Workflows export customization options

## 1.3.21 (2025-11-24)

- Show workflow help button in ribbon panel as well

## 1.3.20 (2025-11-20)

- Workflows new buttons-based navigation

## 1.3.19 (2025-10-30)

- Fitting with formulas support

## 1.3.18 (2025-09-11)

- Hisotry and export fixes for workflows

## 1.3.17 (2025-09-09)

- Save focused RFV tab in session storage

## 1.3.16 (2025-09-08)

- Fix funcCall reactive updates edge case

## 1.3.15 (2025-09-05)

- Fix optimizer returned funcCall redraw

## 1.3.14 (2025-09-04)

- Small bug fixes

## 1.3.13 (2025-08-15)

- Fix DataFrame input viewers updates

## 1.3.12 (2025-08-15)

- Optimization view support primary params filtering
- Optimization view support support per script defaults

## 1.3.11 (2025-08-08)

- Optimization tool integration fixes

## 1.3.10 (2025-07-28)

- API updated for 1.26.0

## 1.3.9 (2025-07-22)

- Template links syntax
- Allow multiple meta handlers targiting the same io

## 1.3.8 (2025-06-30)

- Fix adding steps from pipeline view

## 1.3.7 (2025-06-09)

- Add StartWorkflow return

## 1.3.6 (2025-06-03)

- Allow all users to use

## 1.3.5 (2025-06-03)

- Help tweaks

## 1.3.4 (2025-05-28)

- Vue and vueuse versions Bump
- Help tweaks

## 1.3.3 (2025-05-21)

- Fix step reordering ui sync

## 1.3.2 (2025-05-12)

- Compute-utils update

## 1.3.1 (2025-05-02)

- Small changes/fixes

## 1.3.0 (2025-04-24)

- Update Compute2 API

## 1.2.2 (2025-04-22)

- Compute-utils update

## 1.2.1 (2025-04-21)

- Fix npm build

## 1.2.0 (2025-03-30)

- API updated for 1.25.0

## 1.1.0 (2025-02-19)

- API updated for 1.24.0
- RFV/CFV improvements
- Routing improvements

## 1.0.0 (2024-12-25)

- Initial release
