# Compute2 stack: gating on the js-api version and on package options

Scope: `packages/Compute2`, `libraries/compute-utils`, `libraries/webcomponents`, `libraries/webcomponents-vue`.

## Gating on the js-api version

| Where | API | What is gated |
|---|---|---|
| `libraries/utils/src/format-version-utils.ts:16-18` | `semver.coerce(grok.shell.build.client.version)`, `semver.gte(v, '1.27.7')` | Default float format: `G7` on js-api >= 1.27.7, `#0.###` otherwise; cached per session |
| `libraries/utils/src/format-version-utils.ts` (`isClientAtLeast`) | same `semver` compare, parameterised | Shared helper for the rule-expressions gate below |
| `libraries/webcomponents-vue/src/InputForm/utils.ts:6,107` | `getDefaultFloatFormat()` | Consumes the gate above for input-form number formatting |

No other code in scope reads `grok.shell.build.client.version` or `grok.shell.build.server` to branch on it.

## Gating on js-api methods (feature detection)

| Where | Check | What is gated |
|---|---|---|
| `packages/Compute2/src/utils.ts:405-414` | `typeof window.grok_View_Get_IsPinned / grok_View_Set_IsPinned / grok_View_Pin === 'function'`, then `typeof view.pin === 'function'` | `pinView()`: js-api 1.27.5 replaced `grok_View_Pin` with `grok_View_Set_IsPinned`; probes which interop the platform provides, falls back to `view.pin()` |
| `packages/Compute2/src/project-export.ts:9-21` | `typeof DG.Project.showSaveDialog === 'function' \|\| typeof window.grok_Project_OpenSaveDialog === 'function'` | `canSaveProject()` / `showProjectSaveDialog()`: prefers the js-api method, falls back to the raw Dart binding; the export action is hidden when neither exists |
| `libraries/compute-utils/shared-utils/utils.ts:197-201` | `typeof grok.dapi.groups.currentUserGroups === 'function'` | `getCurrentUserGroups()`: falls back to raw `fetch('/api/groups/all_parents')` on older js-api |
| `libraries/compute-utils/shared-utils/utils.ts:205-210` | `typeof grok.dapi.groups.requestMembership === 'function'` | `requestGroupMembership()`: falls back to a raw `fetch` POST on older js-api |
| `libraries/compute-utils/reactive-tree-driver/src/runtime/rule-expressions.ts` (`hasScriptSupport`) | `typeof grok.functions.scriptSync === 'function' && isClientAtLeast('1.28.0')` | `script`/`scriptVerdict` rule operations behind `visible` and GrokScript `validator` checks: the method exists on older clients but ignores the variables map, so the version decides; unsupported clients skip these checks with one console warning |
| `libraries/compute-utils/reactive-tree-driver/src/runtime/rule-sources.ts` | `typeof call.evalParamValidators === 'function'` on the step's `DG.FuncCall` | Annotation-derived `validators` checks (the script's own validators run by the platform): resolve to no verdicts with one console warning on clients before 1.28.0; `check` links and rules that name their validators call them directly and are not gated |

## Gating on package options

Read through `_package.settings` (`grok_Package_Get_Settings_Sync`, `js-api/src/entities/misc.ts:186`), declared as `properties` in `packages/Compute2/package.json`:

| Property | Read at | What is gated |
|---|---|---|
| `roleOnlyModelFilter` (bool, default `false`) | `packages/Compute2/src/package.ts:103` | `modelCatalogOptions.roleOnlyFilter`: faster role-only model filter, drops the legacy `#model` tag branch |
| `sharingMethod` (`none` \| `workspaces`, default `none`) | `packages/Compute2/src/sharing/sharing.ts:23-35` | `getShareAction()`: share-to-workspace action shown only for `workspaces` and when `isWorkspaceSharingAvailable()` |

Host-package properties read through the deprecated `getProperties()`: `REPORT_BUG_URL` / `REQUEST_FEATURE_URL` in `libraries/compute-utils/function-views/src/custom-function-view.ts:208-212` and `old-views/src/function-view.ts:268-272`.
