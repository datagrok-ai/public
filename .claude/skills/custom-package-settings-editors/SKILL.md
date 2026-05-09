---
name: custom-package-settings-editors
description: Replace Datagrok's auto-generated package settings UI with a custom widget that edits the package's declared properties
---

# custom-package-settings-editors

## When to use

Your package declares `properties` in `package.json` and Datagrok
auto-renders their default editors in the "Settings" pane — but the
defaults don't fit (you need a `userGroups` picker, folder-existence
check, grouping by category, async lookups, etc.). Triggers: "custom
UI for package settings", "settings editor widget", "replace the
auto-generated form", "render userGroups picker inside Settings". For
declaring the properties themselves, see `develop.md → Package
settings` first — this skill assumes they already exist.

## Prerequisites

- A package scaffold (`grok create <Name>`); commands run from the package root.
- Properties declared in `package.json` under `"properties": [...]`
  (`help/develop/develop.md:542-559`); the editor complements them, it
  does NOT replace them (`DG-FACT-121`).
- `datagrok-api` available (`import * as grok / ui / DG`); the article's
  snippet uses `ui`/`DG` without imports.
- Familiarity with `DG.Widget` — same base class as dashboard and
  tooltip widgets (`DG-FACT-119`).

## Steps

1. **Pick a registration form (decorator vs. header comments).**
   The platform discovers settings editors via the function role
   `packageSettingsEditor` (`DG-FACT-118`). Two surfaces in production:
   - **Generic `func` decorator (canonical, used by HitTriage):**
     `@grok.decorators.func({name: '...', tags: ['packageSettingsEditor'],
     meta: {role: 'packageSettingsEditor'}})` on a `static` method of
     `PackageFunctions`. Prefer this — the specialized
     `@grok.decorators.packageSettingsEditor()` exists but no shipping
     package uses it (`DG-FACT-120`).
   - **Header comments (article form):** `//tags: packageSettingsEditor`,
     `//meta.role: packageSettingsEditor`, `//output: widget result`.
   Either way `src/package.g.ts` MUST carry both
   `//tags: packageSettingsEditor` AND `//meta.role: packageSettingsEditor`.

2. **Pick a signature (with or without `propList`).**
   Article shows `packageSettingsEditor(): DG.Widget` (zero args) —
   PowerPack matches (`packages/PowerPack/src/settings-editor.ts:5-9`).
   Every actively-shipping package uses
   `(propList: DG.Property[]): Promise<DG.Widget>` instead — the only
   form that exposes the `package.json` property metadata to the editor
   (`DG-FACT-DRIFT-047`). Default to `propList` unless you truly don't
   need it.
   ```typescript
   // src/packageSettingsEditor.ts
   import * as grok from 'datagrok-api/grok';
   import * as ui   from 'datagrok-api/ui';
   import * as DG   from 'datagrok-api/dg';

   export async function myPackageSettingsEditorWidget(
     propList: DG.Property[],
   ): Promise<DG.Widget> {
     const inputs: DG.InputBase[] = [];
     for (const prop of propList) {
       const cur = prop.get(null) as string | null;       // null = read staged value
       const input = ui.input.string(prop.name, {value: cur ?? ''});
       input.onChanged.subscribe(() => prop.set(null, input.value));
       inputs.push(input);
     }
     return DG.Widget.fromRoot(ui.form(inputs));
   }
   ```
   Expected: build succeeds; widget's root is a `<div>` containing the
   inputs. The `null` scope passed to `prop.get`/`prop.set` is the
   intentional convention — the platform's setter override stages the
   value and persists it on the user's "save" click
   (`packages/HitTriage/src/packageSettingsEditor.ts:71-82`,
   `DG-FACT-121`).

3. **Wire the function on `PackageFunctions`.**
   Match the canonical HitTriage shape — generic `func` decorator,
   `static`, single `propList: DG.Property[]` param, returns
   `Promise<DG.Widget>` (`packages/HitTriage/src/package.ts:350-358`).
   ```typescript
   // src/package.ts
   import * as grok from 'datagrok-api/grok';
   import * as DG   from 'datagrok-api/dg';
   import {myPackageSettingsEditorWidget} from './packageSettingsEditor';
   export const _package = new DG.Package();

   export class PackageFunctions {
     @grok.decorators.func({
       name: 'My package settings editor',
       tags: ['packageSettingsEditor'],
       meta: {role: 'packageSettingsEditor'},
     })
     static async settingsEditor(
       @grok.decorators.param({name: 'propList', type: 'object'})
       propList: DG.Property[],
     ): Promise<DG.Widget> {
       return myPackageSettingsEditorWidget(propList);
     }
   }
   ```
   Expected: `npm run build` regenerates `src/package.g.ts` with a
   wrapper that includes `//tags: packageSettingsEditor`,
   `//input: object propList`, `//output: widget result`,
   `//meta.role: packageSettingsEditor` (compare
   `packages/HitTriage/src/package.g.ts:127-134`).

4. **(Optional) Filter `propList` by category to render a subset.**
   Filter inside the widget to render only one `package.json`
   `"category"` group (HitTriage renders only `DefaultUserGroupsSharing`
   and appends the `Storage` group separately —
   `packages/HitTriage/src/packageSettingsEditor.ts:67,95`). Properties
   NOT rendered by your widget keep using the auto-generated UI in
   their own group (`DG-FACT-121`).

5. **Build, publish, and let the platform auto-register.**
   ```bash
   npm install
   npm run build                  # runs `grok api && grok check --soft && webpack`
   grok publish <host>            # add --release once stable
   ```
   No explicit `register(...)` call is needed — the function role
   `packageSettingsEditor` IS the registration (`DG-FACT-118`).

## Common failure modes

- **Settings pane shows the auto-generated form, not your widget.**
  Header is missing `//meta.role: packageSettingsEditor` OR the decorator
  is missing `meta: {role: 'packageSettingsEditor'}`. The `tags: [...]`
  alone is NOT enough — both `tags` and `meta.role` are written by the
  HitTriage canonical (`packages/HitTriage/src/package.ts:352-353`,
  `DG-FACT-118`). Inspect `src/package.g.ts`: it MUST contain
  `//tags: packageSettingsEditor` AND `//meta.role: packageSettingsEditor`.
- **Widget renders empty / `propList` is `undefined`.** You declared the
  function with no parameters (article form), but tried to use a
  `propList` argument. Either drop the parameter and rely on
  `_package.getProperties()` for reads only, or declare the param with
  `@grok.decorators.param({name: 'propList', type: 'object'})` so the
  generated `package.g.ts` includes `//input: object propList`
  (`DG-FACT-DRIFT-047`).
- **Edits don't persist after "save".** You called `prop.set(value)` or
  `prop.set('default', value)`. Use `prop.set(null, value)` — the `null`
  scope routes through the core's setter override that stages the value
  and writes it on save
  (`packages/HitTriage/src/packageSettingsEditor.ts:82`, `DG-FACT-121`).
- **Read-back returns stale data.** Same root cause — use
  `prop.get(null)`, not `prop.get()` or `_package.settings[name]`; the
  `null` scope routes through the getter override
  (`packages/HitTriage/src/packageSettingsEditor.ts:71`, `DG-FACT-121`).
- **Article snippet won't compile.** The article uses `ui`/`DG` without
  imports and the no-arg signature — works, but blocks per-property
  logic. Add `import * as DG from 'datagrok-api/dg';` (plus `ui`,
  `grok`) and prefer the `propList` form.
- **Inputs visually clipped.** The Settings pane host applies
  `overflow: hidden`. Workaround: pass `ui.form(inputs, {style:
  {overflow: 'visible'}})` and defer
  `parentElement.style.overflow = 'visible'` after mount
  (`packages/HitTriage/src/packageSettingsEditor.ts:92,131`).

## Verification

- `npm run build` exits `0`; `src/package.g.ts` contains an exported
  wrapper with `//tags: packageSettingsEditor`,
  `//meta.role: packageSettingsEditor`, and (if you used the `propList`
  form) `//input: object propList` + `//output: widget result`.
- `grok publish <host>` exits `0`.
- In Datagrok, open **Browse → Apps**, click your package, expand the
  **Settings** pane on the right context panel — your widget renders
  in place of the auto-generated form.
- Edit a value, click **Save**, reload the page, re-open Settings —
  the new value persists.
- (Multi-category) Properties whose `category` your widget filters out
  still render via the auto-generated form in their own group.

## See also

- Source articles:
  - `help/develop/how-to/packages/custom-package-settings-editors.md`
  - `help/develop/develop.md:542-569` (declaring properties + the
    "special editor function" reference).
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-118` … `DG-FACT-121` and drift `DG-FACT-DRIFT-047`.
- Reference packages:
  - `packages/HitTriage/src/packageSettingsEditor.ts:66-133` +
    `packages/HitTriage/src/package.ts:350-358` — canonical decorator,
    `propList`-driven, async, `userGroups` + folder-path with validation.
  - `packages/HitTriage/src/package.g.ts:127-134` — codegen output.
  - `packages/PowerPack/src/settings-editor.ts:5-41` — minimal no-arg
    `extends DG.Widget` form (article shape; export header commented out).
- Related skills:
  - `column-tooltip` (sibling — `DG.Widget` rendered in the column
    tooltip pane, same base class, different host).
