# Compute2 changelog

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
