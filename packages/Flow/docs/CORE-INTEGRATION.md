# Flow ↔ Core Integration Map

Verified against source 2026-07-13. Complements [FIRST-CLASS-INSIGHTS.md](FIRST-CLASS-INSIGHTS.md)
(which describes the *pre*-implementation platform mechanics; the gaps listed there in §1.6-1.7 are
now filled and committed — core `d421759058`, `5101c2c1d8`, `e0e2ad3904`, `4bb6e64217`; public
`77cea89d1e`+). Paths relative to the monorepo root.

## The five seams

### 1. Script language (execution)

- [src/package.ts](../src/package.ts) `flowScriptHandler` — annotation-registered (the
  `templateScript` meta contains `//` lines that don't survive the decorator generator):
  `meta.role: scriptHandler`, `scriptHandler.language: flow`, `.extensions: flow`,
  `.commentStart: //`, `.codeEditorMode: javascript`, `.editorFunction: Flow:flowScriptEditor`.
- Core registration at login func-sync: `core/client/xamgle/lib/src/features/functions/func_sync.dart:114-116`
  → `ScriptHandler.register(new PackageScriptHandler.create(f))`
  (`core/shared/grok_shared/lib/src/scripting/script_handler.dart:146-216`; language const
  `ScriptLanguage.Flow` at :14, `editorFunc` read at :161).
- Execution: `Script.runImpl` → `ScriptHandler.forLanguage('flow').run(call)` → package func with
  `call.toJs()` → [FlowEntityHandler.run](../src/entity/flow-entity-handler.ts): parse `.flow` body →
  compile to clean JS on a detached editor → run → copy outputs back via `setParamValue` (the JS
  FuncCall wraps the same Dart call — the Pyodide contract).
- **Client-only**: package handlers never register server-side → cron/`grok s` runs of flow scripts
  throw "Handler for flow is not registered". Registration is login-time only (mid-session publish
  needs a reload).

### 2. Editor seam (`scriptHandler.editorFunction`)

- Option const: `core/shared/ddt/lib/src/func/func_roles.dart:150-155` (`ScriptHandlerEditorFunc`).
- Consumed by `openCustomScriptEditor` — `core/client/xamgle/lib/src/views/script_view.dart:815-828`:
  handler by language → `Funcs.byName(handler.editorFunc)` → `apply(values:[script])` → sets
  `view.entity` + `/script/<id>` path; any miss falls back to the CodeMirror `ScriptView`.
  Honored at double-click, Edit command, `/script/<id>` URL, and template creation.
- Package side: `flowScriptEditor` → `FuncFlowView.forScript(script)` (bind + load body +
  `updatePath`); Save routes on `boundScript?.id` → `grok.dapi.scripts.save`.

### 3. Entity meta + Browse

- `FlowScriptMeta` — `core/client/xamgle/lib/src/meta/flow_script_meta.dart`
  (`isApplicable: x is Script && x.language == ScriptLanguage.Flow`), registered **before**
  `ScriptMeta` in `entity_renderers.dart` (first-wins). Delegates content to
  `Flow:flowScriptWidget` / `flowScriptPreview`, with install balloons when the package is absent.
  'Open Editor' is the default command. `ScriptMeta` excludes flow from Debug (script_meta.dart:79).
- **Deliberately NOT a JS `DG.ObjectHandler`** — it would be inserted first and shadow the Dart
  meta, losing Sharing/Activity/Chats panes and default-command wiring (plan decision D4).
- Browse: no top-level Flows node; a **Flows group under My stuff** via `ProjectMeta.FLOWS`
  (`project_meta.dart`, `#flow` tag smart-filter source; Scripts group excludes flows), plus a
  `FlowsView` gallery (`DataSourceCardView<Script>`, permanentFilter `language="flow"`).

### 4. Creation scripts / data-sync

- Tags: `Tags.CreationScript = '.script'`, `.DataSync`, `.VariableName` —
  `core/shared/ddt/lib/src/utils/metadata.dart:222+`.
- Recording: `DataHistory.logDataFrameCreationScript`
  (`core/client/xamgle/lib/src/features/functions/data_history.dart`); replay on project open:
  `ProjectDataSyncHandler` in `core/client/xamgle/lib/src/shell/shell_project.dart` (lines ordered
  by the `//{"timestamp"}` comment); single-table: `TableInfo.execDataSync`.
- Core → Flow delegation: `Funcs.byRole('creationScriptEditor')`
  (`core/client/xamgle/lib/src/shell/creation_script_editor.dart:17-34`) →
  [openCreationScriptFlowDialog](../src/package.ts) (script + tableIds). Save splits the graph into
  one script per table (`emitCreationScriptsForTables`; owner matched by the `.VariableName` tag)
  → `TableInfo.saveCreationScript` (`core/shared/grok_shared/lib/src/table_info.dart:142-149`,
  JS at `public/js-api/src/entities/table-info.ts:38`).
- **Multi-output binding (2026-07-08, core `5b647dbcdd`)**: trailing JSON comment carries
  `{"timestamp", "call": <guid>, "output": <paramName>, "outputIndex"?}`; single decoder
  `FuncCall.decodeComment`/`commentMeta` (`core/shared/ddt/lib/src/func/func_call.dart:1005-1021`).
  SetVar binds the named output; project reopen dedups by call GUID; toolbox refresh updates all
  siblings. Plan: `core/docs/plans/multi-output-creation-scripts.md`.
  **Flow's emitter/importer must preserve these keys when re-serializing.**

### 5. Files & deployment

- `fileViewer: ffjson` / `fileViewer: flow` open share files in the editor.
- `scripts/*.flow` auto-register as flow Script entities on publish. Server fix required
  (`deploy_demo_projects.dart`): package handlers don't exist server-side, so unknown extensions
  were parsed with `#` and mis-stamped `javascript`; fixed by sniffing `//` from content.
- `detectCommentStart` (`scripting.dart:394-405`) only knows registered handlers — not usable for
  package-declared languages on the server.

## JS-API surface Flow leans on

`grok.dapi.scripts` (dapi.ts:142-145), `grok.dapi.spaces.id(...).addEntity`, `grok.functions.parse`
(importer), `Func.prepare(params).toString()` → Dart `toConsole`/`valueToString` (creation-script
emitter's single source of truth — needs a live backend), `TableInfo.saveCreationScript`,
`DG.Script.create(js).prepare().call()` (execution), `d4-before-run-action` global event (the
func-editor dialog hijack — see CLAUDE.md Autorun section for the race).

## Constraints

- Body ceiling: `scripts.script varchar(262144)`; bodies ride the startup func-sync payload — keep
  ffjson compact.
- `.flow` body single-writer invariant: only `flowScriptText` writes header+JSON together.
- The committed `scripts/*.flow` demos are generated from `files/*.ffjson` and can drift — the
  "Flow: bundled flow scripts" test category locks them.
