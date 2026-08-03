# Flow as a First-Class Entity — Codebase Insights

Companion to [FIRST-CLASS-PLAN.md](FIRST-CLASS-PLAN.md). This file captures *how the platform works*
in every area the plan touches, with file:line anchors, so a fresh session can implement the plan
without re-doing the research. Paths are relative to the monorepo root (`core/...`, `public/...`).

---

## 1. TL;DR — the load-bearing facts

1. **`Script` is a `Func` is an `Entity`.** `class Script extends Func with GrokJsObject, $Script,
   TagsMixin, AuthorMixin, PackageEntityMixin, EmbeddableEntity`
   ([scripting.dart:804](../../../../core/shared/grok_shared/lib/src/scripting/scripting.dart)).
   Everything a Func gets — params, FuncCall lifecycle, run history, caching, scheduling, sharing,
   spaces membership, startup sync — a Script (and therefore a Flow-as-Script) gets for free.
2. **`language` is an open string end-to-end.** DB column `scripts.language varchar(256)`, no CHECK
   (`core/server/db/init_db.sql:661`); `ScriptingService.saveScript` does **no** language validation
   (`core/server/datlas/lib/src/services/scripting/scripting_service.dart:63-74`). Language only
   matters at *run time*, when `Script.runImpl()` does `ScriptHandler.forLanguage(language).run(call)`
   (scripting.dart:880-895).
3. **Packages can already register a whole new script language.** `PackageScriptHandler`
   (`core/shared/grok_shared/lib/src/scripting/script_handler.dart:136-205`) wraps a package function
   with role `scriptHandler` + `scriptHandler.*` meta options; registered during func sync at
   `core/client/xamgle/lib/src/features/functions/func_sync.dart:114-118`. It provides: language id,
   file extensions, icon, template, CodeMirror mode, a custom **parser function**, a **vectorization
   function**, and a **run** delegate (`packageHandlerFunc.apply(values: [call.toJs()])`). Reference
   implementation: Pyodide (`public/packages/Pyodide/src/package.ts:349-369`). There is even a
   `create-script-handler` skill. **This is ~60 % of "first-class Flow" already built.**
4. **Core already delegates to the Flow package today.**
   `core/client/xamgle/lib/src/shell/creation_script_editor.dart:17-33`:
   `hasCreationScriptFlowEditor() => Funcs.byRole('creationScriptEditor').isNotEmpty`, then
   `funcs.first.prepare(paramValues: {...}).run()` — with an "Install the Flow package…" balloon when
   absent. This guard/delegate shape is the house pattern.
5. **The entity-handler ("meta") system is dual-sided.** Dart `EntityMeta<T>`
   (`core/client/d4/lib/src/common/entity_renderer.dart:149-308`) drives icons, cards, tooltips,
   context panel, Browse tree nodes, previews, commands. JS `DG.ObjectHandler`
   (`public/js-api/ui.ts:1685-1847`) can be registered from a package and is inserted **first** in
   the Dart registry (`grok_api.dart:972`, `first: true`) — it then wins `EntityMeta.forEntity()`
   resolution for anything its `isApplicable` matches. But JS handlers cannot supply: display name,
   `DbEntityMeta` accordion panes (Sharing/Activity/Chats), Browse-tree child nodes, drag-drop,
   grid properties, or the default double-click command (capability matrix in §5.3).
6. **The main gap in core**: a script's *editor* is hardwired to the CodeMirror `ScriptView` — there
   is no per-language editor hook. `PackageScriptHandler` can set the CodeMirror *mode*, not replace
   the editor view. That single seam (a `scriptHandler.editorFunction` option) is the one genuinely
   new core feature Flow needs.
7. **The main gap in the Flow package**: flows persist as locally downloaded `.ffjson` files. No
   save-to-server, no entity, no sharing (admitted in
   [COMPETITIVE-ANALYSIS.md](COMPETITIVE-ANALYSIS.md) ~line 160). The serializer, compiler, and
   editor are all reusable building blocks (§8).

---

## 2. The scripting stack

### 2.1 Classes

| Class | Location | Role |
|---|---|---|
| `Func` | `core/shared/ddt/lib/src/func/func.dart:30` | base: params, options, `run()` lifecycle (:340-420), `prepare()` (:214), `apply()` (:240), static `Func.eval` (:625) |
| `Script` | `core/shared/grok_shared/lib/src/scripting/scripting.dart:804` | adds `language`, `script` (body), `sample`; `runImpl` dispatches to handler (:880-895) |
| `FuncCall` | `core/shared/ddt/lib/src/func/func_call.dart:7` | one execution; `setInput`, typed results, serialization |
| `ScriptHandler` | `core/shared/grok_shared/lib/src/scripting/script_handler.dart:16` | per-language strategy: `run`, `parseCode`, template, icon, extensions, `codeEditorMode`, capabilities |
| `PackageScriptHandler` | script_handler.dart:136-205 | adapter wrapping a package func into a handler |
| `ScriptParser` | scripting.dart:24+ | parses `#name:/#input:/#language:` headers; **delegates to `handler.parseCode` when `hasCustomParser`** (scripting.dart:91, 116-119) |

There are **no per-language Script subclasses** in Dart — all per-language behavior is the handler
(strategy pattern) looked up by the `language` string. Per-entity-type behavior is the `EntityMeta`
chain (`ScriptMeta → _FuncMeta → DbEntityMeta → EntityMeta`), dispatched by `isApplicable`
(`ScriptMeta.isApplicable => x is Script`, `core/client/xamgle/lib/src/meta/script_meta.dart:7`).

### 2.2 Handler registry and dispatch

- Registry: `static Map<String, ScriptHandler> _handlers`, `register()/forLanguage()/forExtension()`
  (script_handler.dart:18, 76-85). `forLanguage` **throws** `'Handler for X is not registered'`.
- Client registration (`core/client/xamgle/lib/src/features/scripting_plugin.dart:13-18`): server
  handlers from `startupData.scriptHandlers` + `JsScriptHandler` + `GrokScriptHandler`; package
  handlers added during func sync (func_sync.dart:114-118). Node/JS host mirror:
  `core/shared/grok_shared/web/grok_shared.dart:296-298`.
- Server registers only built-ins (jupyter_script_handler.dart:219-234). **Package handlers exist
  client-side only** → a flow script executed server-side (cron, `grok s`) would throw at
  `getCapabilities()`. This is an accepted limitation initially (see plan §9).
- `PackageScriptHandler.create(f)` reads options (`FuncOptions.*` in
  `core/shared/ddt/lib/src/func/func_roles.dart:119-148`): `scriptHandler.language` (required),
  `.extensions` (required), `.commentStart`, `.templateScript`, `.codeEditorMode`,
  `.parserFunction`, `.vectorizationFunction`, plus `icon` resolved against the package root.
  `requiresServer=false` forced. **Note:** `FuncOptions.ScriptHandlerName`
  (`'scriptHandler.friendlyName'`, func_roles.dart:129) exists but is consumed nowhere —
  `friendlyName` falls back to `capitalize(language)` (script_handler.dart:50-52). The run closure
  requires the handler func's **first** param to be a `funccall` (by convention the only one;
  extras are not rejected) and calls it with `call.toJs()` (script_handler.dart:180-185) — the JS
  `DG.FuncCall` wraps the *same* Dart call, so outputs set via `scriptCall.setParamValue(...)`
  propagate back (production-proven by Pyodide's `setOutputs`,
  `public/packages/Pyodide/src/package.ts:276-304`). The parser closure resolves the parser func by
  `Funcs.byName` with an explicit throw when missing (:154-163).
- `Script.prepareAsync()` re-parses the code through the custom parser when `hasCustomParser`
  (scripting.dart:932-942) — so a language can derive params from a non-header body format.
  `Script.fromCodeAsync` (:832) exists because package parsers are async JS.

### 2.3 Script UI today

- Editor: `ScriptView extends TableView` (`core/client/xamgle/lib/src/views/script_view.dart:4`),
  `PATH = '/script'`, per-script URL `/script/<id>`; CodeMirror mode from `handler.codeEditorMode`
  (:373); Run = `script.prepareAsync → currentCall.edit()` (:716-721); Save = `dapi.scripting.save`
  (:660). Registered in `View.typeFactories`/`View.urlFactories`
  (`core/client/xamgle/lib/src/views/view.dart:181, 244`).
- **Editor-open choke points** (all funnel into one of two helpers in script_view.dart:757-769):
  - `addScriptView(script)` — gallery double-click (`views/scripts_view.dart:96`), ScriptMeta
    `'Edit...'` command (`features/scripting_plugin.dart:51-57`).
  - `addTemplateScriptView(language)` — `Tools | Scripting | <lang>` commands, generated per
    registered handler in `ScriptMeta.initCommands()` (script_meta.dart:85-94, gated by
    `Permission.CREATE_SCRIPT`).
  - `ScriptView.fromPath` for `/script/<id>` and `/script/<language>` URLs (script_view.dart:80-105).
- Scripts gallery: `ScriptsView extends DataSourceCardView` (`views/scripts_view.dart:4-45`),
  `makeDataSource() => dapi.scripting`, card renderer `renderScriptCard` (:49-98). Browse node:
  `Platform/Functions/Scripts` → preview factory at `browse_panel.dart:200`.
- Icon chain: `ScriptMeta.renderIcon → renderFuncCustomIcon → Script.getImagePath()` → the
  **registered** handler's `iconPath` (scripting.dart:862-867; default
  `/images/entities/<language>.png` only when a handler exists without `meta.icon`,
  script_handler.dart:41) → final fallback `Icons.fa('scroll')` (script_meta.dart:64). With no
  registered handler `getImagePath()` returns null → scroll icon. `Icons.image` renders via CSS
  `background-image`; a 404 yields a silently blank icon, not a fallback (icons.dart:139-147).
- `ScriptMeta.renderAccordion` panes: Details / Script / Run / Activity / Sharing / Chats
  (script_meta.dart:18-33). Browse preview uses inherited `_FuncMeta.renderPreview` → run-form view
  (`meta/function_meta.dart:44-62`), *not* the editor.

### 2.4 Execution flow (condensed)

`call.run()` → `Func.run()` (func.dart:340-420: cache check → `queueInterceptor` → param
resolve/validate → `runImpl`) → `Script.runImpl` → `ScriptHandler.forLanguage(language).run(call)`.
Client-side handlers (grok/javascript/pyodide/package) run in-browser; `ServerScriptHandler` opens a
WebSocket to Datlas, which may re-route to a Jupyter worker via AMQP + grok_pipe (full detail:
`core/docs/SCRIPTING_FLOW.md`). For a `PackageScriptHandler` there is no queue involvement at all.

---

## 3. Entity machinery (persistence, types, clients)

- **Every entity = a row in `entities`** (`core/server/db/init_db.sql:20-30`) + a row in its subtype
  table, upserted automatically by the ORM (`core/server/dinq/lib/src/repository_query.dart:1045-1051`).
  Scripts: `funcs` (base, discriminator `source='script'`) + `scripts` (subtype: `language`,
  `script varchar(262144)`, `sample`) — init_db.sql:655-668, 1494-1519.
- **Type registry**: `initTypes()` maps fixed UUIDs → `$X.TYPE` strings
  (`core/shared/grok_shared/lib/src/entity.dart:3-43`; `$Script.TYPE = "Script"` at :24), seeded to
  DB on boot (`core/server/datlas/lib/src/deploy/init_entity_types.dart`). Generated `$Script`
  (grok_shared.g.dart:7461+) carries `TYPE`, `getTableName() => "scripts"`, DB field mappings —
  produced by `@AutoProperties()` + `dart build.dart`.
- **Server**: `ScriptsRepository extends _FuncsRepository<Script>`
  (`core/server/datlas/lib/src/repositories/scripts_repository.dart`), registered under key
  `"Script"` in `DbContext` (db_context.dart:198). Routes: `ScriptingRouter` `/scripts/...`
  (`core/server/datlas/lib/src/routers/scripting.dart`). `_FuncsRepository.save` also caches the
  full func JSON into `entity_jsons` for startup sync (data_actions_repository.dart:345-353).
- **Client**: `dapi.scripting` = `ScriptingClient extends HttpDataSource<Script>`
  (`core/shared/grok_shared/lib/src/http_client/scripting_client.dart:27-62`), self-registered into
  the polymorphic save dispatch via `EntitiesClient.registerClient(Script, ...)`. JS side:
  `grok.dapi.scripts` (`public/js-api/src/dapi.ts:142-145`), entity class `DG.Script`
  (`public/js-api/src/entities/func.ts:126`).
- **Startup func sync**: all viewable non-`'function'` funcs (including script bodies via the
  `entity_jsons` cache) are shipped to every client at login (`STARTUP_FUNCS_QUERY`,
  data_actions_repository.dart:30-65) and registered into the in-memory `Funcs` registry
  (func_sync.dart:7-196). **Implication: flow bodies ride the startup payload — keep them compact.**
- **Cost of a brand-new entity type** (rejected in the plan, kept for reference): ~8 coordinated
  layers — DB migration + init_db, shared model + codegen, `initTypes()` UUID, repository +
  DbContext, service + router, Dart client + registry, JS API + interop, xamgle meta/UI — and it
  forfeits all Func machinery (params, runs, galleries, scheduling). Template exemplar: Notebook
  (`notebook.dart`, `notebooks_repository.dart`, `notebooks_client.dart`).

---

## 4. Core → JS package delegation: the house patterns

### 4.1 Discovery / guard primitives (Dart)

- `Funcs.byName('Pkg:fn')` → `null` when absent — THE existence check
  (`core/shared/ddt/lib/src/func/funcs.dart:134-154`).
- `Funcs.find({packageName, functionName, tags, role, meta})` — Dart twin of JS
  `DG.Func.find` (funcs.dart:108-128; `meta:` matches `func.options`).
- `Funcs.byNamespace('Flow').isNotEmpty` — "is the package installed" (funcs.dart:156-158).
- `Funcs.byRole(role)` (funcs.dart:205) — used by the existing Flow hook.
- `Funcs.handle(condition, handler)` — applies to all current **and future** registrations
  (funcs.dart:311-315) — order-independent wiring (used to turn `view`/`viewer` funcs into
  View/Viewer factories, `core/client/xamgle/lib/src/interop/js_special_functions.dart:3-22`).
- `ClientPackageFunc` lazy-loads the package bundle transparently on first `apply()`
  (`core/client/xamgle/lib/src/interop/js_func.dart:291-315`); `runSync` requires prior init.

### 4.2 Call + result marshalling

Core→JS calls pass live handles, never JSON: `f.apply(values: [...], processed: true)` →
`jsu.callFuncWithDartParameters` → `public/js-api/src/functions.ts:407-413` → `paramsToJs`/`toJs`
(`public/js-api/src/wrappers_impl.ts`). Can pass: scalars, Maps/Lists, DataFrames, entities
(anything with `GrokJsObject`/`className`, incl. `Script` → `DG.Script` and `FuncCall` →
`DG.FuncCall`), ProgressIndicator. Can receive: `Types.VIEW` → Dart `View` (add via `addView`),
`Types.WIDGET` → `JsWidgetBase` (`.root`), DataFrame/Column, scalars, Maps; multi-output = JS
returns an object keyed by output-param names (js_func.dart:13-32).

### 4.3 Precedents, ranked by fit for Flow

1. **Notebooks** — first-class entity + package-provided view. `NotebookMeta.getView` →
   `View.byType('Notebook')` whose factory is injected by the package's `tags: view` function
   (js_special_functions.dart:3-11 + `public/packages/Notebooks/src/package.js:524-532`). Guards:
   `checkEnabled: () => Funcs.byNamespace('Notebooks').isNotEmpty ? null : 'Install Notebooks
   plugin first'` (notebook_meta.dart:161); runtime `Funcs.byName(...) == null → throw`
   (notebook.dart:79-87).
2. **Custom function editors** — `func.options['editor']` / `dialogFunc` resolved via
   `getDialogFunction` (`features/functions/func_param_editor.dart:233-238`), swapped in for the
   function view and Run pane (`function_meta.dart:51-62, 118-139`), with graceful
   `'Install X package for better experience'` fallback (function_meta.dart:268-272). Consumed by
   the shell run/edit pipeline at `shell.dart:989-1009` (where `run.aux['forceEditParameters']`
   forces the param dialog — the flag near the user's earlier selection).
3. **PackageScriptHandler** (§2.2) — language-level delegation; Flow's backbone.
4. **ML engines** — multi-role function family keyed by one meta value (`mlname`/`mlrole`),
   each optional (`predictive_modeling_engines.dart:134-159`) — the naming model if Flow ever needs
   several hooks (`meta.flowRole: preview|widget|...`).
5. **fileViewer / file-handler** — extension-keyed view providers with checker funcs
   (`funcs.dart:193-203`, `features/file_editors.dart:209-278`) — how `.ffjson`/`.flow` files in
   shares open today.
6. **UserReportMeta.openReport** — the minimal meta-command delegate:
   `byName → null-guard(Balloon) → prepare(paramValues) → run()` (user_report_meta.dart:89-96).

---

## 5. The entity-handler ("meta") system

### 5.1 Dart side

`EntityMeta<T>` override surface (entity_renderer.dart): `renderIcon` (:77), `renderMarkup` (:24),
`renderTooltip` (:30), `renderCard` (:33), `renderProperties`/`renderAccordion` (:65-68, context
panel), `renderView` (:71, `/e/<id>` route), `openPreview`/`renderPreview` (:59-62, Browse preview),
`addEntityNode` (:36-57, tree node), `getName`/`getCaption`, `getById`/`refresh`, `init()` +
`regCommand(name, run:, check:, isDefault:)` (:224-231 — `isDefault: true` = double-click action,
executed via `ui.bind` at `d4/lib/src/widgets/ui.dart:839-843`), `getDropActions`, `properties`.
`DbEntityMeta` adds server-entity panes (sharing, chats, rename dialogs)
(`core/client/xamgle/lib/src/meta/db_entity_meta.dart`).

Registration: `initEntities()` (`meta/entity_renderers.dart:3-90`) — **order matters**;
`EntityMeta.forEntity` returns the first handler whose `isApplicable` is true
(entity_renderer.dart:208-213). A meta for a Script *subset* (e.g. `language == 'flow'`) must be
registered **before** `ScriptMeta`. Per `meta/CLAUDE.md`: new meta = `part` statement in `meta.dart`
+ entry in `initEntities()`.

Consumption sites (all go through `forEntity`): context panel
(`features/property_panel.dart:79-88`), tooltips (ui.dart:832-834), Browse tree
(`browse_panel.dart:403-417`), Browse preview (`browse_panel_preview.dart:140-147`), `/e/<id>`
routing (routing.dart:466-479), gallery cards (`DataSourceCardView`).

### 5.2 JS side

`DG.ObjectHandler` (`public/js-api/ui.ts:1685-1847`): `type`, `isApplicable`, `getById`, `refresh`,
`renderIcon/Markup/Tooltip/Card/Properties/View`, `async renderPreview → View`, `toMarkup`,
`markupRegexp`, `init`, `registerParamFunc(name, run)` (adds Actions/context-menu commands via
`grok_RegisterParamFunc`). `ObjectHandler.register(meta)` → Dart
`EntityMeta.addMeta(new JsEntityMetaProxy(jsMeta), first: true)` (grok_api.dart:972) — the proxy
forwards the members listed in `js_utils.dart:265-282`. Packages register in their `#init` function
(e.g. Chem `MpoProfileHandler`, `public/packages/Chem/src/package.ts:226-232`).

### 5.3 What a JS ObjectHandler can and cannot influence

| Capability | JS? |
|---|---|
| Icon, markup, tooltip, card, context-panel element, `/e/<id>` view, Browse preview, getById/refresh, helpUrl, actions via `registerParamFunc`, grid cell renderer | **Yes** |
| Display name (`getName`), `DbEntityMeta` accordion panes (Sharing/Activity/Chats), Browse-tree child nodes / lazy loading, drag-drop actions, `properties`/grid columns, default double-click command | **No — Dart only** |

**Design consequence**: registering a JS ObjectHandler for flow scripts would *shadow* the Dart
`ScriptMeta` (first-wins) and lose the platform panes and default-command wiring. The plan therefore
keeps a thin **Dart** `FlowScriptMeta` as the owner and has it *call into* the Flow package for the
flow-specific content (see plan §4, decision D4).

---

## 6. Browse, opening, creation, URLs, icons

- **Tree**: hardcoded skeleton in `BrowsePanelTreeBuilder.refresh()`
  (`browse_panel_tree.dart:9-307`), permission-gated per section; leaves get content via
  `nodePreviewFactories[path]` + `urlNodesMapping[urlKey]` registered in the `BrowsePanel`
  constructor (browse_panel.dart:177-242). Recipe for a new node is in
  `browse_panel/CLAUDE.md`. Scripts node: `funcsNode.addItem('Scripts', icon: Icons.fa('scroll'))`
  (:206) → `ScriptsView`.
- **Preview pipeline**: `BrowsePanelPreview._resolveNodeViews` — path factories → special cases →
  `EntityMeta.forEntity(entity).openPreview(...)` (browse_panel_preview.dart:109-148). Preview tabs
  are temporary; double-click pins (browse_panel.dart:152-164).
- **Galleries**: extend `DataSourceCardView<T>` (`views/data_source_card_view.dart:7`) — supply
  `makeDataSource()`, `meta`, `commands` (`commandsButtonName = 'New'` pattern in
  scripts_view.dart:19-23). Per `views/CLAUDE.md`, new Browse CRUD pages must extend it.
- **Creation**: per-handler `Tools | Scripting | <lang>` commands auto-generated
  (script_meta.dart:85-94); gallery New button = `commands` list; context "New X..." =
  `regCommand` on the parent meta (e.g. `DataConnectionMeta.cmdAddQuery`).
- **URL wiring checklist** (from `views/CLAUDE.md`, for a `/flows` gallery route): permission
  constant + tree leaf + `nodePreviewFactories` + `urlNodesMapping` + `View.typeFactories` +
  `View.urlFactories`. Entity URLs: `Entity.getUrl()` default `/browse/<nqName>`; Func override
  `/func/<nqName>` (run view); ScriptView URL `/script/<id>`.
- **Icons**: `Icons.fa/svg/image` (`d4/lib/src/widgets/icons.dart`); package-relative PNG via
  `options['icon']`/`meta.icon`; language default `/images/entities/<language>.png`.

---

## 7. Spaces (what Flow inherits for free)

- Space = `Project` with storage (`project.dart:57`); membership = `ProjectRelation` rows;
  permissions cascade via `project_relations_all` (`repository_query.dart:148-172`). Full doc:
  `core/docs/SPACES.md`.
- **Containment gate**: `Project.canContain(x)` requires `PackageEntityMixin` (project.dart:328) —
  **Script has it**, so flows-as-scripts are containable, draggable into spaces, and auto-published
  into the user's namespace project on save (repository_query.dart:1052-1060, 1295-1343).
- Space tree categories: hard-coded switch in `ProjectMeta._getProjectStats`
  (project_meta.dart:126-154; `$Script.TYPE → 'Scripts'`) + per-category data sources
  (:420-504). Flows show under **Scripts** unless a language-aware split is added (optional polish).
- Space children rendering uses `EntityMeta.forEntity(e).addEntityNode(...)` and
  `openPreview` — a registered `FlowScriptMeta` is automatically honored inside spaces and
  `NamespaceView` galleries (namespaces_view.dart:122-127, 169).

---

## 8. The Flow package today

### 8.1 Registered functions (complete)

| Function | Annotations | Notes |
|---|---|---|
| `funcflowApp(path?)` | `tags: app`, `meta.role: app`, `output: view` | returns `new FuncFlowView()`; `path` currently ignored |
| `viewFuncFlow(file)` | `meta.role: fileViewer`, `meta.fileViewer: ffjson` | `.ffjson` file double-click → editor |
| `flowFromCreationScript(script)` | `output: view` | builds flow from a creation script |
| `openCreationScriptFlowDialog(script, tableIds, show)` | `meta.role: creationScriptEditor`, `meta.includeInFlow: false` | **the existing core→Flow delegation target** (dialog + per-table save-back) |
| `info`, `testDialog` | — | dev |

No detectors, no file-handler, no save-to-server, no string→view opener, no headless compile/run/
summarize functions callable from Dart.

### 8.2 Persistence format (`.ffjson` v2, `src/serialization/flow-schema.ts:6-75`)

```ts
{ version: '2.0', name, description, author, created, modified,
  nodes: [{id, typeName, label, description?, collapsed?, pos, properties, inputValues}],
  connections: [{id, source, sourceOutput, target, targetInput, waypoints?}],
  annotations?: [...],
  metadata: { settings: {scriptName, scriptDescription, tags} } }
```

`serializeFlow/deserializeFlow/downloadFlow/loadFlowFromFile` in `flow-serializer.ts`. Loader
rejects version ≠ '2.0'; unknown node `typeName`s are skipped; node ids are regenerated on load.
Save/Open today are **local disk only** (funcflow-view.ts:841-865). The only server-side flows are
two templates in the package `files/` share.

### 8.3 Compilation targets

1. **JS script** (`emitScript`, script-emitter.ts:39): standard Datagrok JS script with
   `//input://output:` header derived from Input/Output nodes (incl. qualifiers + captions), body of
   `await grok.functions.call(...)` steps; optional instrumented mode (events + `__ff_stash`
   live-value registry). **Not** reverse-importable.
2. **Creation script** (grok-language cascade): bidirectional (`buildCreationScriptGraph` pure
   importer; `emitCreationScript` — requires live backend since it uses Dart
   `prepare().toString()`), but **lossy** (drops viewers, JS-only nodes, layout, qualifiers).
   Conclusion: **`.ffjson` is the only lossless flow representation.**

### 8.4 Reusable building blocks for the entity integration

`FuncFlowView.loadFromJson(json)` (funcflow-view.ts:1095), `loadFromCreationScript`,
`serializeFlow`, `emitScript` (+ `compileGraph`, `sliceUpTo`, `validateGraph`),
`ExecutionController.runInstrumented/previewNodeData/produceTableForNode`, headless detached-editor
pattern (`makeEditor`, tests/test-utils.ts:18-23 — enables compile/summarize without a visible
view), `summarizeFlow/summarizeNode` (pure — context panel/card material), `FlowSettings` (maps to
entity name/description/tags). Run of the clean compiled script is fully standalone:
`DG.Script.create(script).prepare().call()`.

### 8.5 Known internal quirks that affect the integration

- DG-func node factories register ~100 ms after view construction (`registerAllFunctions` in a
  setTimeout, funcflow-view.ts:105-112); `loadFromJson` defers to compensate. An entity-open path
  should await registration deterministically rather than racing.
- `funcflowApp`/`viewFuncFlow` toggle shell windows in a 200 ms setTimeout — a hack to keep out of
  the view; fine to keep for now.
- The view's Save currently downloads; Save semantics must fork on "entity-bound vs scratch".

---

## 9. Constraints and gotchas

1. **Dart 1.x** in core: `new` constructors, no generic function types in fields/collections
   (analyzer stack-overflow — use `Function`), `Map.removeWhere` missing, `Modal.onCancel` is a
   field. See root `CLAUDE.md` "Dart 1.x gotchas".
2. **Codegen**: adding `FuncOptions` constants (ddt) or `ScriptHandler` fields (grok_shared)
   requires `dart build.dart` in dependency order: ddt → grok_shared → datlas → xamgle → d4.
3. **Body size**: `scripts.script` is `varchar(262144)` — 256 K chars ceiling for header + ffjson.
4. **Startup payload**: script bodies ship to all clients at login via `entity_jsons` (§3). Compact
   serialization matters; a body-stripping optimization would be a separate core task.
5. **Server-side execution**: package script handlers don't exist server-side; scheduled/`grok s`
   runs of flow scripts will fail with "Handler for flow is not registered". Acceptable v1
   limitation; document it.
6. **Meta ordering**: `FlowScriptMeta` must register before `ScriptMeta` in `initEntities()`
   (entity_renderers.dart) or use `addMeta(first: true)`.
7. **JS ObjectHandler first-wins**: if the package ever registers a `DG.ObjectHandler` matching
   flow scripts, it will shadow the Dart meta entirely (grok_api.dart:972) — don't do both.
8. **`ScriptParser` tolerance vs. `ScriptView` fragility**: unknown languages *parse* fine by
   default (`allowInvalidLanguage: true`, scripting.dart:42, 80, 105, 832), and the editor's
   debounced validation disables Run for them (script_view.dart:402-403). **But opening the editor
   for an unregistered language currently breaks**: `ScriptView.init()` calls
   `ScriptHandler.forLanguage(script.language).codeEditorMode` unguarded (script_view.dart:373) and
   `forLanguage` throws; `addScriptView`'s `init().then(...)` has no catchError (:107-118) so the
   view hangs on the loader, and `/script/<id>`'s catchError mislabels it "Script not found"
   (:90-98). A missing-handler guard in `ScriptView.init()` is required work, not existing behavior.
9. **Handler registration is login-time only**: package script handlers register in the startup
   func-sync loop (func_sync.dart:110-117); the runtime entity-sync listeners (:84-97) only call
   `registerFunc` and skip the handler branch, and `ScriptMeta.initCommands()` runs once after
   `initFuncs()` (`xplorer_init.dart:180-181`). A package published mid-session gets its language,
   template command, and editor seam only after a page reload.
10. **`grok api` / publish parsing**: server extracts func metadata from `src/package*.ts`
   annotations by regex (packages_service.dart:156-275); `meta.*` lines become `func.options`
   (scripting.dart:151-153). Never `.npmignore` `src/` or `webpack.config.js`.
