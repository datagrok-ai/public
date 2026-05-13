---
name: extensions
description: Pick the right Datagrok extension point for a plugin behavior, scaffold the function with datagrok-tools, and declare the matching function-role tag
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# extensions

## When to use

You want a Datagrok plugin to add platform-level behavior — an app, a
viewer, a column action, a file format, a context-panel widget, a column
type detector — and need to know which extension point fits, which
function-role tag activates it, and where the dedicated how-to lives.
Triggers: "extend Datagrok with X", "register a file type", "add
something to the column context panel", "what role tag do I use for Y".

## Prerequisites

- A package scaffold (`grok create <Name>`); `datagrok-tools` (`grok`)
 on PATH.
- `~/.grok/config.yaml` populated for at least one host
 (`grok config show <alias>`).
- The package builds locally (`webpack`) and publishes
 (`grok publish <alias>`).

## Steps

1. **Identify the extension type from the intended behavior.**
 The overview article enumerates eleven extensible surfaces; the
 `function-roles.md` reference enumerates the role tags that activate
 each. Map your intent against this table — pick the narrowest row
 that fits (e.g. "render molecules in a grid cell" → `cellRenderer`,
 not the broader `viewer`).
 ```bash
 ls help/develop/how-to/ # the per-extension how-to dirs
 ```
 Expected: a single chosen row from the table below.

 | Behavior | Role tag | How-to article | Reference package | Knowledge |
 |---|---|---|---|---|
 | Standalone app / dashboard | `app` | `how-to/apps/build-an-app.md` | `packages/HitTriage` | `DG-FACT-001` |
 | Context-panel info card | `panel` | `how-to/ui/add-info-panel.md` | `packages/PowerPack` | (function-roles.md §Info-panels) |
 | Custom viewer (subclass `JsViewer`) | (output `viewer` in package fn) | `how-to/viewers/develop-custom-viewer.md` | `packages/Charts` | `DG-FACT-186..197` |
 | File viewer | `fileViewer` + `meta.fileViewer: <ext>` | `how-to/files/create-custom-file-viewers.md` | `packages/EpsViewer`, `packages/Plates` | `DG-FACT-073..078` |
 | File handler (import) | `fileHandler` (decorator-config `ext` REQUIRED) | `how-to/files/file-handlers.md` | `packages/nmrium` | `DG-FACT-083..086` |
 | File exporter | `fileExporter` (description REQUIRED) | `how-to/files/file-exporters.md` | `packages/Chem`, `packages/Bio` | `DG-FACT-079..082` |
 | Folder content preview | `folderViewer` | `how-to/files/folder-content-preview.md` | `packages/ClinicalCase`, `packages/Plates` | `DG-FACT-088..094` |
 | Cell renderer | `cellRenderer` (+ `meta.cellType`) | `how-to/grid/custom-cell-renderers.md` | `packages/PowerGrid`, `packages/Chem` | `DG-FACT-099..105` |
 | Column tooltip | generic `func` + `meta.role: tooltip` | `how-to/grid/column-tooltip.md` | `packages/Bio`, `packages/Chem` | `DG-FACT-095..098` |
 | Semantic-type detector | `semTypeDetector` (in `detectors.js`) | `how-to/functions/define-semantic-type-detectors.md` | `packages/Chem` | (function-roles.md §Semantic-type-detectors) |
 | Custom filter (subclass `DG.Filter`) | n/a — registered via viewer factory | `how-to/viewers/custom-filters.md` | `packages/Widgets` | (no fact yet) |
 | Custom script handler | `scriptHandler` (+ four required `meta.scriptHandler.*`) | `how-to/scripts/custom-script-handlers.md` | `packages/Pyodide` | `DG-FACT-173` |
 | Settings editor | `packageSettingsEditor` | `how-to/packages/custom-package-settings-editors.md` | `packages/HitTriage`, `packages/PowerPack` | `DG-FACT-118..119` |
 | Package init hook | `init` | function-roles.md §Package-initialization | `packages/Chem` | — |
 | Autostart hook | `autostart` (optional `meta.autostartImmediate: true`) | function-roles.md §Autostart | `packages/Chembl`, `packages/PowerPack` | — |
 | Custom views / accordion sections | n/a — `grok.shell.dockManager.dock(...)`, `grok.shell.addView(...)` | `how-to/views/custom-views.md`, `advanced/ui.md#accordions` | `packages/PowerPack` | `DG-FACT-213..217` |

2. **Scaffold from a `grok add` template when one exists.**
 Five extension types ship a `datagrok-tools` template
 (`function-roles.md` §applications, §info-panels,
 §package-initialization, §semantic-type-detectors). For everything
 else, copy the function header from the dedicated how-to.
 ```bash
 cd <package-name>
 grok add app <Name> # role: app
 grok add function init <Pkg>Init # role: init
 grok add function panel <Name> # role: panel
 grok add detector <SemTypeName> # role: semTypeDetector → detectors.js
 grok add viewer <ViewerName> # JsViewer subclass under src/
 ```
 Expected: a stub file under `src/` (or `detectors.js` for the
 detector); the next `grok link && webpack` regenerates `package.g.ts`
 with the correct annotations preserved.

3. **Declare the function role explicitly in the source.**
 One role per function; do not combine `init, autostart`. Decorator
 form is canonical; codegen emits header-form into `package.g.ts`.
 ```typescript
 // header form
 //meta.role: panel
 //condition: semType == "Molecule"
 //input: string smiles {semType: Molecule}
 //output: widget result
 export function showMolecule(smiles: string) { /* … */ }

 // decorator form
 @grok.decorators.panel({condition: 'semType == "Molecule"'})
 static showMolecule(
 @grok.decorators.param({options: {semType: 'Molecule'}}) smiles: string,
 ): DG.Widget { /* … */ }
 ```

4. **(Detector / fast-path autostart only) place the function in `detectors.js`.**
 `semTypeDetector` lives in `detectors.js` so the platform can run
 detectors without loading the main bundle. Do NOT re-export it from
 `src/package.ts` (triggers a full-bundle load on first cell read).
 All other roles belong in `src/package.ts`.

5. **Build and publish.**
 ```bash
 webpack && grok publish <alias>
 ```
 Expected: webpack exits `0`; `grok publish` exits `0`; the platform
 loads the new function on next refresh.

6. **Verify the function was registered with the right role.**
 In the platform JS console (or any package code path):
 ```javascript
 const found = await DG.Func.find({meta: {role: '<role-token>'}})
.filter((f) => f.package?.name === '<PackageName>');
 console.log(found.map((f) => f.name));
 ```
 Expected: a non-empty array containing your function name. The UI
 also exposes role-filtered search at
 `<host>/functions?q=%23<role>`. To rule out a registration bug
 versus a runtime bug, reload with
 `<host>?initPackageFunctions=false` — your function must NOT fire
 (function-roles.md tip).

## Common failure modes

- **Function not discovered after publish.** `//meta.role:` missing or
  mistyped. `grok link` to regenerate `package.g.ts`, verify the role,
  re-publish.
- **Compound `init, autostart` annotation.** Declare ONE role per
  function; split if both behaviors are needed.
- **Role-token casing wrong.** Tokens are camelCase: `cellRenderer`,
  `semTypeDetector`, `fileExporter`, `fileViewer`, `fileHandler`,
  `folderViewer`, `scriptHandler`, `packageSettingsEditor`. Copy
  verbatim from step 1's table.
- **Detector defined in `src/package.ts` instead of `detectors.js`.**
  Detector code never runs because the main bundle is not loaded until
  invocation.
- **Required companion annotation missing.** `fileViewer` needs
  `meta.fileViewer: <ext>` (`DG-FACT-074`); `fileExporter` needs
  `description` (`DG-FACT-080`); `cellRenderer` needs `meta.cellType`
  (`DG-FACT-100/103`); `scriptHandler` needs four
  `meta.scriptHandler.*` keys (`DG-FACT-173`). Cross-check against the
  per-extension how-to.

## See also

- Source: `help/develop/packages/extensions.md`,
  `help/develop/function-roles.md`, per-extension how-tos in step 1.
- Knowledge: `DG-FACT-001`, `073-078`, `079-082`, `083-086`, `088-094`,
  `095-098`, `099-105`, `118-119`, `173`, `186-197`, `213-217`.
- Related skills: `build-an-app`, `routing`, `custom-cell-renderers`,
  `create-custom-file-viewers`, `file-exporters`, `file-handlers`,
  `folder-content-preview`, `column-tooltip`, `custom-script-handlers`,
  `custom-package-settings-editors`, `custom-views`,
  `manipulate-viewers`, `register-identifiers`, `data-enrichments`,
  `publish-packages`.
