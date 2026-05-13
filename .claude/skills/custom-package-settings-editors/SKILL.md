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
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

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

- One or more properties declared in `package.json` under
  `"properties": [...]` (`develop.md:542-559`). The editor *complements*
  the declarations — the platform reads the metadata, passes them as
  `DG.Property[]`, and your widget renders the UI (`DG-FACT-121`).
- Familiarity with `DG.Widget` — subclass it or
  `DG.Widget.fromRoot(htmlElement)` (`DG-FACT-119`).

## Steps

1. **Declare the properties in `package.json` that your editor will render.** Without entries here, there is nothing to bind to — the editor *complements* the declarations (see DG-FACT-121). Use `category` for grouping and `inputType` for non-default input widgets.

   ```json
   "properties": [
     {"name": "view", "propertyType": "string", "inputType": "userGroups",
      "defaultValue": "", "category": "DefaultUserGroupsSharing",
      "friendlyName": "View", "nullable": true}
   ]
   ```

2. **Implement an editor returning `DG.Widget` and register it with role `packageSettingsEditor`.** Use `@grok.decorators.func` on a `static` inside `PackageFunctions`, with both `meta.role` and `tags` set (see DG-FACT-118, DG-FACT-120). The function takes one input — `propList: DG.Property[]` — that the platform fills from `package.json` properties.

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

   Read with `prop.get(null)`, stage changes with `prop.set(null, value)` — the `null` scope is required (the core routes through a per-pane buffer); the platform persists on "Save" (see DG-FACT-121). Subclassing `DG.Widget` is equivalent (see DG-FACT-119).

3. **Build — the role IS the registration.** `npm run build` runs `FuncGeneratorPlugin` and regenerates `src/package.g.ts`. Do NOT hand-edit.

   ```bash
   npm install && npm run build && grok check --soft
   grok publish <host>   # add --release once stable
   ```

## Common failure modes

- **Custom editor never replaces the default UI.** Wrapper missing `//meta.role: packageSettingsEditor` (camelCase, exact-match — see DG-FACT-118).
- **`propList` is empty inside the editor.** `package.json` `properties` block is missing or malformed — re-check Step 1 (see DG-FACT-121).
- **Changes typed in the editor don't persist.** Body called `prop.set(value)` — must be `prop.set(null, value)` (see DG-FACT-121).
- **`@grok.decorators.packageSettingsEditor(...)` compiles but doesn't fire.** That specialized decorator has an empty body. Switch to the generic `func` decorator form (see DG-FACT-120).
- **Article snippet won't compile.** Add `import * as ui`/`import * as DG` and use the decorator form, not the legacy `//meta.role:` header.

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
