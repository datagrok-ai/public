# Compute2 stack: platform-version and options gating

## Summary

- No code in Compute2, compute-utils, compute-api, webcomponents or webcomponents-vue branches on the server version; the only version gate the stack reaches is `getDefaultFloatFormat()` in `@datagrok-libraries/utils`, which compares `grok.shell.build.client.version` with `semver`.
- Platform capabilities are gated by feature detection: `typeof api.grok_* === 'function'` on `window` globals and `typeof (obj as any).method === 'function'` on js-api objects, each with a fallback to the older path.
- Behaviour toggles come from two package properties (`roleOnlyModelFilter`, `sharingMethod`), `PipelineConfiguration` flags, `Driver`/`StateTree` construction options (`mockMode`, `defaultValidators`, `batchLinks`) and URL search params; no env or build-time flags gate runtime behaviour.
- All five packages depend on `datagrok-api` as `workspace:^` (js-api is 1.27.11); platform-served libraries are pinned once in `build-config/platform-deps.json` and enter the pnpm catalog via `.pnpmfile.cjs`.

## Platform version gating

| Where | API | What is gated |
|---|---|---|
| `libraries/utils/src/format-version-utils.ts:16-18` | `semver.coerce(grok.shell.build.client.version)`, `semver.gte(v, '1.27.7')` | Default float format: `G7` on js-api >= 1.27.7, `#0.###` otherwise; cached per session |
| `libraries/webcomponents-vue/src/InputForm/utils.ts:6,107` | `getDefaultFloatFormat()` | Consumes the gate above for input-form number formatting (only consumer in scope) |

Direct reads of the server or js-api version inside Compute2, compute-utils, compute-api, webcomponents, webcomponents-vue: none found.

js-api surface for version and build info (`js-api/src`):

| API | Location | Notes |
|---|---|---|
| `grok.shell.build.client: ComponentBuildInfo` | `shell.ts:23-25,46` | `{branch, commit, date, version}` via `grok_Shell_GetClientBuildInfo()` |
| `grok.shell.build.server: ComponentBuildInfo` | `shell.ts:25` | Field is initialised as an empty `new ComponentBuildInfo()`; no TS code assigns it, no consumer found in `packages/`, `libraries/` or `help/` |
| `ComponentBuildInfo` | `dapi.ts:82-87` | |
| `grok.dapi.admin.getServiceInfos()` | `dapi.ts:501`, `ServiceInfo` at `dapi.ts:605-613` | Service status list (name, enabled, status, type); used by `libraries/tutorials/src/tutorial.ts:169` and `packages/UsageAnalysis/src/widgets/usage-widget.ts:40`, not for version gating |
| `DG.Package.version` | `entities/misc.ts:94` | Package version, not platform version |

Prior art outside the stack: `packages/UsageAnalysis/src/release/data.ts:174` and `src/test-track/app.ts:79` read `grok.shell.build.client.version` for display only. No `compareVersions`/`versionCompare`/`semver` usage exists in `js-api/src`.

## Feature detection

| Where | Check | What is gated |
|---|---|---|
| `packages/Compute2/src/utils.ts:405-414` | `typeof api.grok_View_Get_IsPinned/grok_View_Set_IsPinned/grok_View_Pin === 'function'` on `window`, then `typeof view.pin === 'function'` | `pinView()`: js-api 1.27.5 replaced `grok_View_Pin` with `grok_View_Set_IsPinned`; probes which interop the platform provides, falls back to `view.pin()` |
| `packages/Compute2/src/project-export.ts:9-21` | `typeof (DG.Project as any).showSaveDialog === 'function' \|\| typeof window.grok_Project_OpenSaveDialog === 'function'` | `canSaveProject()` / `showProjectSaveDialog()`: prefers the js-api method, falls back to the raw Dart binding; the export action is hidden when neither exists |
| `libraries/compute-utils/shared-utils/utils.ts:197-201` | `typeof (grok.dapi.groups as any).currentUserGroups === 'function'` | `getCurrentUserGroups()`: falls back to raw `fetch('/api/groups/all_parents')` on older js-api |
| `libraries/compute-utils/shared-utils/utils.ts:205-210` | `typeof (grok.dapi.groups as any).requestMembership === 'function'` | `requestGroupMembership()`: falls back to raw `fetch` POST on older js-api |
| `libraries/compute-utils/reactive-tree-driver/src/runtime/rule-sources.ts` | `typeof call.evalParamValidators === 'function'` on the step's `DG.FuncCall` | Annotation-derived `validators` checks (the script's own validators run by the platform): resolve to no verdicts with one console warning on clients before 1.28.0; `check` links and rules that name their validators call them directly and are not gated |
| `libraries/compute-utils/function-views/src/run-comparison-view.ts:62` | `(window as any).grok_TableView(...)` | Unconditional raw Dart binding, no detection |

Reference pattern in `libraries/u2/src/dg/funcs`:

| Where | Check | What is gated |
|---|---|---|
| `param-rules.ts:9-13,39-60` | `typeof api.grok_ScriptSync === 'function'` (`api = globalThis`) | Rule expressions run through the raw `grok_ScriptSync` global; when absent, the rule keeps the previous state / skips validation and records a warning |
| `param-tables.ts:17-28` | `typeof api.grok_Property_Get !== 'function'` | Reading a Dart property by name returns `null` on clients without the binding |
| `param-tables.ts:48,81` | `'func' in v`, `typeof v.names === 'function'` | Shape checks on values received from the platform |

Not capability checks, listed for completeness: `globalThis.initialURLHandled` (`packages/Compute2/src/apps/RFVApp.tsx:133-139`, `components/TreeWizard/TreeWizard.tsx:339-347`, `libraries/compute-utils/function-views/src/custom-function-view.ts:219-223`) is a process-wide once-flag for start-URL handling; `libraries/compute-utils/old-views/function-view.ts:140` uses `grok.shell.getVar('isLoaded')` for the same purpose; `navigator.hardwareConcurrency` checks in `function-views/src/fitting/worker/{pool.ts:254,executor.ts:58}` size the worker pool.

`'x' in grok/DG` checks, optional-chaining probes on platform APIs, and checks of `grok_*` globals other than the ones above: none found in scope.

## Package options and config flags

Package settings (`_package.settings`, backed by `grok_Package_Get_Settings_Sync`, `js-api/src/entities/misc.ts:186`):

| Property (`packages/Compute2/package.json` `properties`) | Read at | What is gated |
|---|---|---|
| `roleOnlyModelFilter` (bool, default `false`) | `packages/Compute2/src/package.ts:103` | `modelCatalogOptions.roleOnlyFilter`: faster role-only model filter, drops the legacy `#model` tag branch |
| `sharingMethod` (`none` \| `workspaces`, default `none`) | `packages/Compute2/src/sharing/sharing.ts:23-35` | `getShareAction()`: share-to-workspace action shown only for `workspaces` and when `isWorkspaceSharingAvailable()` |

Other settings sources: `libraries/compute-utils/function-views/src/custom-function-view.ts:208-212` and `old-views/src/function-view.ts:268-272` read `REPORT_BUG_URL` / `REQUEST_FEATURE_URL` from the host package's deprecated `getProperties()`; `old-views/src/pipeline-view.ts:605-610` persists help-panel state in `grok.userSettings`. `grok.shell.settings`, `process.env`, `import.meta.env`, `__DEV__`, `NODE_ENV` and `mode === 'development'` in runtime source: none found. The only build-time switch is `packages/Compute2/rspack.config.js:33-38` (`env.enable_vue_dev_tools` -> `mode: 'development'`), which does not reach application code.

`PipelineConfiguration` flags (`libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration.ts`):

| Flag | Declared | Consumed | What is gated |
|---|---|---|---|
| `enableHistory` (FuncCall step) | `:263` | `runtime/StateTreeNodes.ts:168`; `packages/Compute2/src/components/TreeWizard/TreeWizard.tsx:485-488` | Per-step opt-in save-to-history and history panel |
| `disableHistory` | `:287` | `StateTreeNodes.ts:370`; `Compute2/src/utils.ts:171`; `components/PipelineView/PipelineView.tsx:85,152,250` | Hides history controls, blocks `couldBeSaved` |
| `disableDefaultExport` | `:288` | `TreeWizard.tsx:398,879` | Suppresses the built-in Excel export |
| `forceNavigate` | `:285` | `StateTreeNodes.ts:372`; `Compute2/src/utils.ts:113` | Includes otherwise-skipped pipeline nodes in navigation |
| `disableUIControlls/Adding/Removing/Dragging` (`NestedItemContext`) | `:293-296` | `StateTreeNodes.ts:464-468`; `Compute2/src/utils.ts:169`; `TreeNode.tsx:295-301`; `TreeWizard.tsx:653,687`; `package.ts:653-702` | Add/remove/drag controls on dynamic steps |
| `isActionStep` | `:318` | `StateTreeNodes.ts:397`; `PipelineView.tsx:86-97,152,169`; `TreeWizard.tsx:958` | Placeholder step that only shows actions |
| `runOnInit` (link) | `:93,174` | `runtime/LinksState.ts:100,391` | Links run at tree init |

Driver / StateTree construction options:

| Option | Set at | Consumed | What is gated |
|---|---|---|---|
| `mockMode` | `Driver.ts:45` (ctor, default `false`; `Compute2/src/composables/use-reactive-tree-driver.ts:33` uses the default); `mockMode: true` only in `packages/LibTests/src/test-utils.ts` and tests | `runtime/StateTree.ts:107,122,295,583`; `runtime/StateTreeFactory.ts:179-190` | Uses `FuncCallMockAdapter` instead of real FuncCalls; save throws; `runStep` accepts mock results |
| `defaultValidators` | `Driver.ts:299,330,338` (always `true`); default `false` in `StateTreeFactory.ts:22,46,111` | `runtime/LinksState.ts:41,106,122`; `runtime/links-dependencies.ts:89-133` | Auto-created per-IO validators and initial meta run |
| `batchLinks` | `Driver.ts:300,331,339` (always `true`) | `runtime/Link.ts:200` | Batched link execution for batchable links |

`localValidation`: a prop of `RichFunctionView` (`packages/Compute2/src/components/RFV/RichFunctionView.tsx:260,330,609`), set to `true` only by the standalone `RFVApp` (`src/apps/RFVApp.tsx:291`). When true the run gate is the Dart form's own validity; otherwise it is the driver's `callState.isRunnable`. It is the one switch that separates the standalone RFV from a workflow step.

URL params: `packages/Compute2/src/url-inputs.ts:48-126` (`parseUrlInputs`, `buildInputsUrl`) map search params to FuncCall inputs; `apps/RFVApp.tsx:75-149` and `components/TreeWizard/TreeWizard.tsx:119-387` use `useUrlSearchParams` (`id`, `currentStep`) to restore runs and step position.

## Version pinning and dependencies

| Package | `datagrok-api` | Other pinning fields |
|---|---|---|
| `packages/Compute2/package.json` | `workspace:^` | `properties` (above), `sources: ["common/vue.js"]`, `meta: {url: "/Modelhub", dartium: false}`, `grok.testDependencies`; no `peerDependencies`, `engines` or platform-version field |
| `libraries/compute-utils/package.json` | `workspace:^` | none |
| `libraries/compute-api/package.json` | `workspace:^` | none |
| `libraries/webcomponents/package.json` | `workspace:^` | none |
| `libraries/webcomponents-vue/package.json` | `workspace:^` | none |

`js-api/package.json` version: 1.27.11. `workspace:^` resolves to the in-repo js-api at build time and is rewritten to `^1.27.11` on publish; there is no declared minimum platform (server) version anywhere in the stack.

Platform-served libraries: `build-config/platform-deps.json` pins rxjs 6.6.7, dayjs 1.11.23, cash-dom 8.1.5, wu 2.1.0, exceljs 4.4.0, html2canvas 1.4.1, openchemlib 8.21.0, vue 3.5.42, ngl 2.5.0, codemirror 5.65.21. `.pnpmfile.cjs` (workspace root) merges these into the default pnpm catalog so `"vue": "catalog:"` resolves to the platform version; every entry is a default bundler external; `grok check` warns on packages that pin their own version (`build-config/README.md:42-50`). Root `pnpm-workspace.yaml` `catalog:` adds only typescript, `@types/node`, `@types/wu`. Opt-out per package: `externals: {codemirror: false}` in `rspack.config.js`.

## Reusable helpers

| Helper | Location | Fit for a version gate |
|---|---|---|
| `getDefaultFloatFormat()` | `libraries/utils/src/format-version-utils.ts:14-21` | Working example of the pattern (`semver.coerce` + `semver.gte` on `grok.shell.build.client.version`, cached); `@datagrok-libraries/utils` already depends on `semver ^7.7.4`, so a generic `isApiAtLeast(min)` could sit beside it |
| `NodePackagesDataSource.compareVersions(a, b)` | `tools/bin/utils/node-dapi.ts:719-727` | Numeric dot-group compare, but lives in the Node-only `datagrok-tools` CLI; not importable from browser packages |
| `getCurrentUserGroups()` / `requestGroupMembership()` | `libraries/compute-utils/shared-utils/utils.ts:197-210` | Template for method-presence detection with an older-API fallback |
| `pinView()` | `packages/Compute2/src/utils.ts:400-416` | Template for probing raw `grok_*` interop globals |

Generic `compareVersions`/`versionCompare` helper in `@datagrok-libraries/utils` or `js-api`: none found.

## Recommendation for gating a new js-api capability

Follow the pattern the stack already uses: detect the capability at the call site (`typeof (obj as any).method === 'function'` for js-api methods, `typeof (window as any).grok_X === 'function'` for raw interop) and keep the old path as the fallback, as in `project-export.ts:9-21` and `shared-utils/utils.ts:197-210`. Feature detection is preferred because `workspace:^` ties the bundle to whatever js-api it was built against, while the running client may be older or newer; a version compare is only justified when the capability exists on both sides but behaves differently (the `G7` case), and then it should reuse `semver` on `grok.shell.build.client.version` as `format-version-utils.ts` does rather than add a new comparison helper. Do not gate on `grok.shell.build.server` (never populated in TS, no consumers) or on package.json fields (no minimum-platform field exists).
