# Compute2 stack: platform requirement and package options

Scope: `packages/Compute2`, `libraries/compute-utils`, `libraries/webcomponents`, `libraries/webcomponents-vue`.

## Platform requirement

The stack requires platform 1.28.0 or later and doesn't branch on the platform version. The docs describe no way to
block installing a package on an older platform: the `datagrok-api` constraint in `package.json` states the
requirement, and nothing on the client enforces it. A package on an older platform fails on the first missing API.

## Fallbacks kept for Compute v1

`libraries/compute-utils/shared-utils/utils.ts` is shared with Compute v1's Model Hub, so these stay:

| Where | Check | What is gated |
|---|---|---|
| `getCurrentUserGroups()` | `typeof grok.dapi.groups.currentUserGroups === 'function'` | Falls back to a raw `fetch('/api/groups/all_parents')` on older js-api |
| `requestGroupMembership()` | `typeof grok.dapi.groups.requestMembership === 'function'` | Falls back to a raw `fetch` POST on older js-api |

## Gating on package options

Read through `_package.settings` (`grok_Package_Get_Settings_Sync`, `js-api/src/entities/misc.ts:186`), declared as `properties` in `packages/Compute2/package.json`:

| Property | Read at | What is gated |
|---|---|---|
| `roleOnlyModelFilter` (bool, default `false`) | `packages/Compute2/src/package.ts:103` | `modelCatalogOptions.roleOnlyFilter`: faster role-only model filter, drops the legacy `#model` tag branch |

Host-package properties read through the deprecated `getProperties()`: `REPORT_BUG_URL` / `REQUEST_FEATURE_URL` in `libraries/compute-utils/function-views/src/custom-function-view.ts:208-212` and `old-views/src/function-view.ts:268-272`.
