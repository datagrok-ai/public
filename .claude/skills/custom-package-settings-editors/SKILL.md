---
name: custom-package-settings-editors
version: 0.1.0
description: |
  Replace Datagrok's auto-generated Settings pane for a package with a
  custom-built widget — sliders, dropdowns, user-group pickers, folder
  pickers, validation, grouping — so plugin authors can present
  user-tunable preferences (compute modes, thresholds, default sharing
  groups, storage folders) in a richer form than the default
  one-property-per-input list. Produces a `static` method on
  `PackageFunctions` decorated with `@grok.decorators.func({meta:
  {role: 'packageSettingsEditor'}, ...})` returning a `DG.Widget`.
  Use when asked to "build a custom UI for plugin preferences",
  "render my own form for package properties", or "replace the
  auto-generated Settings pane with sliders and dropdowns".
triggers:
  - custom ui for plugin preferences
  - own form for package properties
  - replace the auto-generated settings pane
  - richer widgets in the settings panel
  - group preferences with validation
  - control how plugin preferences look
allowed-tools:
  - Read
  - Write
  - Edit
  - Bash
harness-authored: true
---

# custom-package-settings-editors

## When to use

Your package declares properties in `package.json` and you want the
right-hand *Settings* pane (under the package's context panel) to render
custom inputs — user-group pickers, folder pickers, validated paths,
grouped sections, sliders, color pickers — instead of the default
one-input-per-property list. Real-world phrasings: "default sharing
groups for new campaigns", "a folder browser for the storage root",
"color palette + threshold sliders in one panel".

## Prerequisites

- A package scaffold (`grok create <Name>`); commands run from package root.
- `datagrok-api` imported — `import * as grok`, `import * as DG`,
  `import * as ui` (article omits these —
  `custom-package-settings-editors.md:13-25`).
- One or more properties declared in `package.json` under
  `"properties": [...]` (`develop.md:542-559`). The editor *complements*
  the declarations — the platform reads the metadata, passes them as
  `DG.Property[]`, and your widget renders the UI (`DG-FACT-121`).
- Familiarity with `DG.Widget` — subclass it or
  `DG.Widget.fromRoot(htmlElement)` (`DG-FACT-119`).

## Steps

1. **Declare the properties in `package.json` that your editor will
   render.** The platform generates the default UI from this block and
   passes the same `DG.Property[]` into your custom editor; without
   entries here, there is nothing to bind to (`DG-FACT-121`). Each entry
   needs at minimum `name`, `propertyType`, `defaultValue`
   (`develop.md:549-559`); use `category` to group properties into
   sections and `inputType` to pick a non-default input widget
   (HitTriage's `userGroups` + grouped category at
   `packages/HitTriage/package.json:67-96`).

   ```json
   "properties": [
     {"name": "view", "propertyType": "string", "inputType": "userGroups",
      "defaultValue": "", "category": "DefaultUserGroupsSharing",
      "friendlyName": "View", "nullable": true}
   ]
   ```

   Expected: `grok check --soft` exits 0 — the properties block parses.

2. **Implement an editor returning `DG.Widget` and register it with
   role `packageSettingsEditor`.** The platform discovers the editor by
   exact-match role — `packageSettingsEditor` (camelCase,
   `js-api/src/const.ts:392-394`, `DG-FACT-118`). The article shows the
   legacy header-comment form (`//meta.role: packageSettingsEditor` —
   `custom-package-settings-editors.md:13-25`); that block is what the
   codegen *emits* into `package.g.ts`, not what you author. The
   authoring surface is `@grok.decorators.func` on a `static` inside
   `PackageFunctions`, with both `meta.role` and `tags` set
   (`DG-FACT-118`, `DG-FACT-120`). The function takes one input —
   `propList: DG.Property[]` — that the platform fills from the
   `package.json` properties. Sync or async return both work.

   ```typescript
   // src/package.ts
   import * as grok from 'datagrok-api/grok';
   import * as DG from 'datagrok-api/dg';
   import * as ui from 'datagrok-api/ui';

   export const _package = new DG.Package();

   export class PackageFunctions {
     @grok.decorators.func({
       name: 'My package settings editor',
       meta: {role: 'packageSettingsEditor'},
       tags: ['packageSettingsEditor'],
     })
     static async settingsEditor(
       @grok.decorators.param({name: 'propList', type: 'object'})
         propList: DG.Property[],
     ): Promise<DG.Widget> {
       const sharing = propList.filter((p) => p.category === 'DefaultUserGroupsSharing');
       const inputs = sharing.map((prop) => {
         const input = ui.input.string(prop.name, {value: prop.get(null) ?? ''});
         input.onInput.subscribe(() => prop.set(null, input.value));
         return input;
       });
       return DG.Widget.fromRoot(ui.form(inputs));
     }
   }
   ```

   Read the current value with `prop.get(null)`, stage changes with
   `prop.set(null, value)`; the platform persists on the user's "Save"
   click (`DG-FACT-121`; HitTriage at
   `packages/HitTriage/src/packageSettingsEditor.ts:66-130`, `:71-82`).
   The `null` `scope` is required — the core overrides the
   getter/setter to route through a per-pane buffer. Filtering
   `propList` by `.category` renders only the properties you styled; the
   rest keep the default rendering. Subclassing `DG.Widget` (article
   form at `custom-package-settings-editors.md:20-24`; PowerPack at
   `packages/PowerPack/src/settings-editor.ts:11-40`) is equivalent.

3. **Build, and verify the auto-emitted wrapper in `package.g.ts`.**
   `npm run build` runs `FuncGeneratorPlugin`, which scans
   `PackageFunctions` and regenerates `src/package.g.ts`. Do NOT
   hand-edit `package.g.ts`. The role IS the registration — no
   `grok.functions.register(...)` call needed.

   ```bash
   npm install && npm run build && grok check --soft
   grok publish <host>   # add --release once stable
   ```

   Expected: `npm run build` exits 0; the wrapper carries both
   `//tags: packageSettingsEditor` and `//meta.role: packageSettingsEditor`
   plus `//input: object propList` and `//output: widget result`
   (compare HitTriage at `packages/HitTriage/src/package.g.ts:127-134`).

## Common failure modes

- **Custom editor never replaces the default UI.** Inspect
  `src/package.g.ts` — the wrapper MUST carry
  `//meta.role: packageSettingsEditor` (camelCase, exact-match);
  `PackageSettingsEditor` or `settings-editor` won't register
  (`DG-FACT-118`).
- **`propList` is empty inside the editor.** The `package.json`
  `properties` block is missing or malformed — without declared
  properties, there is nothing for the platform to pass in
  (`DG-FACT-121`). Re-check Step 1.
- **Changes typed in the editor don't persist after "Save".** The body
  called `prop.set(value)` instead of `prop.set(null, value)` — the core
  overrides the setter and keys the buffer on `null` scope
  (`packages/HitTriage/src/packageSettingsEditor.ts:79-87`).
- **`@grok.decorators.packageSettingsEditor(...)` compiles but the
  editor doesn't fire.** That specialized decorator
  (`js-api/src/decorators/functions.ts:310-316`) has an empty body —
  a codegen marker current production templates don't use. Switch to
  `@grok.decorators.func({meta: {role: 'packageSettingsEditor'}, ...})`
  (`DG-FACT-120`).
- **Article snippet won't compile.** It omits
  `import * as ui`/`import * as DG` and uses the legacy `//meta.role:`
  header — the `package.g.ts` shape, not the authoring surface
  (`custom-package-settings-editors.md:13-25`). Add imports and switch
  to the decorator form.

## Verification

- `npm run build`, `grok check --soft`, and `grok publish <host>` all
  exit 0.
- Regenerated `src/package.g.ts` carries both
  `//tags: packageSettingsEditor` and `//meta.role: packageSettingsEditor`
  (compare `packages/HitTriage/src/package.g.ts:127-134`).
- In Datagrok, open the package's context panel and expand *Settings* —
  your custom widget renders. Tweaking a value and clicking "Save"
  persists; re-opening the pane shows the saved value.

## See also

- Source: `help/develop/how-to/packages/custom-package-settings-editors.md`
  (mirror: `docs/_internal/articles-mirror/how-to/packages/custom-package-settings-editors.md`);
  background: `help/develop/develop.md:542-569`,
  `help/develop/function-roles.md:136-142`.
- Knowledge: `DG-FACT-118` (role token), `DG-FACT-119` (`DG.Widget`
  return), `DG-FACT-120` (decorator surface), `DG-FACT-121`
  (`propList` flow + persistence).
- Reference packages:
  `packages/HitTriage/src/package.ts:350-358` +
  `packages/HitTriage/src/packageSettingsEditor.ts:66-133` +
  `packages/HitTriage/package.json:67-96` — full production form;
  `packages/PowerPack/src/settings-editor.ts:11-40` — minimal subclass.
- Related skills: `custom-cell-renderers`, `column-tooltip` (same
  decorator + role-token + `package.g.ts` registration discipline).
