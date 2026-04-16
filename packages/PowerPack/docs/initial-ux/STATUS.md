# Initial UX — Workspace Widget Status

## Date: 2026-04-15

## What was built

### New files
- `src/spotlight/group-favorites.ts` — `pin()`, `unpin()`, `getMyGroupFavorites()`, `getAdminGroups()`. Storage: per-group projects named `Home: <GroupName>` using `DG.Project.addLink()` + `grok.dapi.permissions.grant()`.
- `src/spotlight/workspace-tab.ts` — `WorkspaceTab` class: card-based layout with greeting, per-group pinned entity cards, "Continue where you left off" recent items, new-user fallback.
- CSS appended to `css/power-pack.css` — card tiles with rounded borders, blue hover glow, flex-wrap layout (classes prefixed `pp-workspace-`).

### Modified files
- `src/spotlight/spotlight-widget.ts` — Added `'Workspace'` as the first tab in `TabControl`.
- `src/package.ts` — Added `onContextMenu` handler: right-click entity → "Add to group favorites" → submenu of groups user is admin of.

## Current state

- TypeScript compiles cleanly (`tsc --noEmit` passes).
- Webpack builds successfully (all 3 new files bundled).
- NOT yet published or visually tested against a running Datagrok instance.

## User modifications after initial implementation

1. **package.ts context menu handler** — user changed entity extraction to `DG.toJs(args?.args?.item?.value)` (tree view node value) instead of `args.args.context`. Added `getEntity()` helper, debug statements (`grok.shell.info`, `console.log`, `debugger`). The debug statements should be removed before shipping.
2. **workspace-tab.ts** — user commented out `this.spotlight.initRecentData()` from the `Promise.all` call (line 114), likely during debugging. The recent data section won't populate until this is re-enabled or data is loaded another way.

## What remains

- Remove debug statements from `package.ts` (`grok.shell.info`, `console.log`, `debugger`).
- Re-enable `initRecentData()` in workspace-tab.ts or confirm alternative data loading path.
- Publish to a Datagrok instance and visually verify:
  - Workspace tab renders as default first tab
  - Cards display correctly for pinned entities
  - "Add to group favorites" context menu works end-to-end
  - Recent items section populates
  - New user (<5 days) sees onboarding content
  - Empty state shows hint text
- Consider adding "Remove from group favorites" (unpin) — context menu on cards or a manage dialog.
