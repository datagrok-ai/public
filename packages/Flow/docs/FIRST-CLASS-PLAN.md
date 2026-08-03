# Flow as a First-Class Entity — Implementation Plan

Companion to [FIRST-CLASS-INSIGHTS.md](FIRST-CLASS-INSIGHTS.md) (read it first — every mechanism
referenced here is anchored there with file:line). Together the two files are self-sufficient to
start implementation in a fresh session.

---

## 1. Motivation

Flow is a visual workflow designer that compiles node graphs into ordinary Datagrok JS scripts.
Functionally it is a *script authoring surface* — but today a flow is a `.ffjson` file downloaded to
the user's disk. It has no identity on the platform: no entity, no sharing, no versioning, no
search, no spaces, no run history, no scheduling, no URL. `docs/COMPETITIVE-ANALYSIS.md` calls this
out as the biggest gap vs. KNIME/Spotfire: *"the opportunity is to bridge Flow into the platform's
existing governance."*

Spotfire's canvas is first-class in Spotfire. Flow must be first-class in Datagrok: users find it in
Browse, create flows from "New", organize them in Spaces, share them like any entity, double-click
to open the visual editor, run them like any function, and see rich previews and context panels —
while the implementation stays in the Flow package (TypeScript), with core (Dart) knowing only a
thin, guarded contract.

## 2. Goal and non-goals

**Goal**: a **Flow entity** = a `Script` with `language: 'flow'`, whose body is the lossless
`.ffjson` document prefixed by a standard annotation header. Core treats it as a script everywhere
(persistence, sharing, spaces, params, FuncCalls, galleries) and delegates flow-specific behavior
(editor view, execution, preview, context-panel content, parsing) to Flow-package functions via
`Funcs.byName`-style lookups, guarded so that a server without the Flow package degrades gracefully.

**Non-goals (v1)**:
- Server-side execution of flows (cron/`grok s`) — package script handlers are client-only today.
- A new DB table / entity type (see decision D1).
- Replacing the `.ffjson` fileViewer or the creation-script dialog — they stay.
- Flow versioning beyond what entities already provide.

## 3. Target UX (acceptance stories)

1. **Create**: `Tools | Scripting | Flow` (auto-appears when the package is present at login — see
   Phase 3.4), a **New**
   button in the Flows gallery, and the Flow app's Start panel → the visual editor opens on a
   template; **Save** asks for a name (and optionally a space) and creates a platform entity.
2. **Find**: Browse → `Platform | Functions | Flows` (own node + gallery with cards); flows also
   show inside Spaces; entity search finds them by name/tags.
3. **Open**: double-click anywhere (tree, card, space, search) → the **visual editor**, never
   CodeMirror. URL `/script/<id>` for a flow also opens the visual editor.
4. **Run**: from the context panel Run pane, `grok.functions.call('user:MyFlow', {...})`, the
   funccall dialog, or the editor — parameters come from the flow's Input nodes; outputs from
   Output nodes; run history recorded like any script.
5. **Preview / context panel**: selecting a flow shows its icon, a plain-language summary
   (`summarizeFlow`), node/link counts, and Details/Sharing/Activity panes; Browse preview shows a
   read-only rendering of the graph.
6. **Organize & share**: drag into a Space, share via the standard dialog — zero new code
   (inherited from Script).
7. **Degrade**: on a client without the Flow package, a flow still lists/opens as a script entity
   (header parses, params visible, icon default), Run and Edit show "Install the Flow package…".

## 4. Design decisions (decision log — review these first)

| # | Decision | Alternatives considered | Rationale |
|---|---|---|---|
| **D1** | **Persist as `Script` row, `language='flow'`** (+ tag `#flow`) | New entity type (new table/repo/router/client/meta — 8 layers) | Zero server changes: `language` is unconstrained end-to-end; CRUD/sharing/search/spaces/params/run-history all key off the Script type. New type forfeits all Func machinery. (Insights §1.2, §3) |
| **D2** | **Body = annotation header + ffjson** (see §5) | Body = ffjson only (needs custom parser for everything); body = emitted JS with ffjson in `options` | Default `ScriptParser` can read name/params even without the package → graceful degradation; ffjson stays the single lossless source of truth; `options` is a worse home (also serialized, less visible) |
| **D3** | **No Dart `FlowScript` subclass** — behavior via `FlowScriptMeta` + `PackageScriptHandler` | `class FlowScript extends Script` | Matches platform architecture: per-language behavior is a handler, per-type UI is a meta (no per-language Script subclasses exist). A subclass needs deserialization dispatch changes in repositories/`EntitiesClient.newInstance` (`#type` is `"Script"`) for near-zero gain. Revisit only if real polymorphism emerges. |
| **D4** | **Dart `FlowScriptMeta` owns dispatch; JS side is one `FlowEntityHandler` class whose methods are exposed as a few role-tagged functions** | Register a `DG.ObjectHandler` via `ObjectHandler.register` | A registered JS handler is inserted `first: true` and *shadows* the Dart meta entirely, losing `DbEntityMeta` panes (Sharing/Activity/Chats), display name, default double-click command, drag-drop (Insights §5.3). Instead the Dart meta stays the owner and calls JS for flow-specific *content*. The user's "one ObjectHandler class" wish is honored organizationally: all JS logic lives in one `FlowEntityHandler` class; the registered functions are 3 thin wrappers. |
| **D5** | **One small generic core seam: `scriptHandler.editorFunction`** — a per-language custom editor view | Flow-specific hooks scattered in ScriptView/ScriptsView/meta; or `editor:` func option (that's the *run dialog* hook, wrong seam) | One option on `ScriptHandler`, honored at the two editor choke points (`addScriptView`, `addTemplateScriptView`) + `/script/<id>` routing. Benefits any future visual/package language, not just Flow. |
| **D6** | **Flows keep entity type `Script`**; Browse gets a dedicated **Flows** node filtered by language; Scripts gallery excludes flows | Distinct EntityType "Flow" | Distinct type ⇒ D1's Option B costs. Filtering by `language='flow'` (`whereSmart`) is cheap and works server-side. Space tree category will read "Scripts" in v1 (optional polish later). |
| **D7** | Execution = handler compiles ffjson → clean JS script on a detached editor → runs it, copying params/outputs through the incoming `FuncCall` | Persist the emitted JS alongside and run that | Single source of truth; emitted JS can drift from the graph. Compile is fast; detached-editor pattern exists (`makeEditor`). |

If any decision is overturned in review, phases below marked with the decision id are affected.

## 5. The `.flow` body format (D2)

Stored in `Script.script` (a.k.a. the body; column `scripts.script`, 256 K chars max):

```
//name: Sales pipeline
//description: Loads, joins and charts quarterly sales
//language: flow
//tags: flow
//input: dataframe orders [Orders table]
//input: int threshold = 10 { min: 0 }
//output: dataframe result
//meta.flowVersion: 2.0
{
  "version": "2.0",
  "name": "Sales pipeline",
  ...full .ffjson v2 document...
}
```

- `commentStart` is `//`. Header lines are **generated** by the Flow serializer from the graph:
  `//input:`/`//output:` derive from Input/Output nodes exactly as `emitScript` already derives
  them (extract the inline header-emission logic — `script-emitter.ts` ~:48, ~:177-209,
  `formatHeaderDefault` :599 — into a shared helper; there is no ready-made `buildHeader`
  function); name/description/tags from `FlowSettings`.
- The JSON payload after the header is the verbatim `serializeFlow` output. Parsers: the default
  `ScriptParser` reads only the header (body opaque) → params correct even without the package;
  the package's `parserFunction` re-derives everything from the JSON (authoritative) and is used by
  `prepareAsync` before every run.
- Invariant: **header and JSON must never disagree** — the only writer is the Flow serializer
  (single `flowScriptText(flow, settings)` helper; saving always regenerates both).
- The `flow` tag enables tag-based filtering alongside the language filter.

## 6. Implementation phases

Ordering is deliberate: each phase is independently shippable and testable.

### Phase 1 — Flow package: language handler + entity persistence (no core changes)

All in `public/packages/Flow`:

1. **`flowScriptText` / `parseFlowBody`** (new `src/serialization/flow-script-format.ts`):
   compose/split header + ffjson; header generation shared with `script-emitter.ts` (extract the
   header builder rather than duplicating).
2. **Handler function** in `src/package.ts`:
   ```ts
   //name: flowScriptHandler
   //input: funccall scriptCall
   //meta.role: scriptHandler
   //meta.scriptHandler.language: flow
   //meta.scriptHandler.extensions: flow
   //meta.scriptHandler.commentStart: //
   //meta.scriptHandler.parserFunction: Flow:parseFlowScript
   //meta.scriptHandler.templateScript: <header + empty-flow json, \n-escaped>
   //meta.icon: files/icons/flow.png
   ```
   Implementation: read `scriptCall.func.script` → `parseFlowBody` → detached `FlowEditor`
   (`makeEditor` pattern from `tests/test-utils.ts`) → `deserializeFlow` → `emitScript` (clean) →
   `DG.Script.create(js)` → `prepare()` with input values copied from `scriptCall` → `call()` →
   set each output via `scriptCall.setParamValue(name, value)` — the JS `DG.FuncCall` wraps the
   *same* Dart call, so values propagate back (production-proven by Pyodide's `setOutputs`,
   `public/packages/Pyodide/src/package.ts:276-304`). Note the handler's `friendlyName` defaults to
   `capitalize('flow') = 'Flow'`; the `scriptHandler.friendlyName` option exists but is wired
   nowhere in core — don't declare it. Await func-node-type registration deterministically before
   deserializing (expose an `ensureFunctionsRegistered(): Promise` in `node-factory.ts` instead of
   racing the 100 ms timer).
3. **Parser function** `parseFlowScript(code): DG.Script` — parse header + JSON, return a
   `DG.Script` with params built from Input/Output nodes (authoritative over the header).
4. **Save-to-server in `FuncFlowView`**: when the view is bound to a script entity, Save →
   regenerate body via `flowScriptText` → `grok.dapi.scripts.save`; when scratch, "Save" dialog
   (name, description, optional space via `grok.dapi.spaces.id(...).addEntity`) then save + bind.
   Keep "Export .ffjson" as a secondary menu item (old Save behavior). Add `Open from platform…`
   (script picker filtered `language = "flow"`).
5. **Editor/view functions** (thin wrappers over one `FlowEntityHandler` class,
   new `src/entity/flow-entity-handler.ts` — D4):
   - `flowScriptEditor(script): view` — `//meta.role: flowScriptEditor` (consumed by the Phase 3
     core seam; until then used internally): `FuncFlowView` bound to the entity,
     `loadFromJson(parseFlowBody(script.script).json)`.
   - `flowScriptWidget(script): widget` — context-panel pane content: `summarizeFlow` text,
     node/link counts, Open/Run buttons.
   - `flowScriptPreview(script): view` — read-only editor (reuse `FuncFlowView` with edit gestures
     disabled, or the editor as-is for v1).
   All marked `//meta.includeInFlow: false` so they stay out of Flow's own toolbox
   (`shouldExcludeFunc` honors it).
6. **Tests** (`src/tests/entity-tests.ts`): format round-trip (`flowScriptText` ↔
   `parseFlowBody`), handler run end-to-end against a live stand (create script entity with
   language flow → `grok.functions.call` it → outputs correct), parser param derivation.
7. `CHANGELOG.md` `## v.next` entries per the repo rule.

**Acceptance**: a flow saved from the editor is a shareable, searchable, space-draggable entity;
`grok.functions.call('u:MyFlow', {...})` executes it with correct params; on a Flow-less client the
entity still lists with correct params and a default icon.

### Phase 2 — Core: the per-language editor seam (D5) + Flow meta

Core changes, all small and generic. Rebuild chain after ddt/grok_shared edits:
`ddt → grok_shared → datlas → xamgle → d4` (`dart build.dart` each).

1. **`FuncOptions.ScriptHandlerEditorFunc = 'scriptHandler.editorFunction'`**
   (`core/shared/ddt/lib/src/func/func_roles.dart`, next to the other `scriptHandler.*` keys).
2. **`ScriptHandler.editorFunc`** field (`core/shared/grok_shared/lib/src/scripting/script_handler.dart`);
   `PackageScriptHandler.create` copies it from options.
3. **Honor it at the choke points** (`core/client/xamgle/lib/src/views/script_view.dart`;
   resolved against source):
   - `addScriptView(script)` (:762-769): if the script's handler has `editorFunc` and
     `Funcs.byName(editorFunc) != null` → `var v = await f.apply(values: [script.toJs()],
     processed: true); addView(v as View);` else current behavior. `addScriptView` is currently
     sync — make it async (callers at scripts_view.dart:96 and scripting_plugin.dart:56 tolerate
     it). This one point covers: gallery double-click, `Edit...` command, the new-script menu
     (:728-732), **and** `addTemplateScriptView` — confirmed to route through `addScriptView`
     (:757-760).
   - **Patched separately** (construct `ScriptView` directly, bypassing `addScriptView`):
     `/script/<id>` URL — urlFactory `views/view.dart:244` → `ScriptView.fromPath` (:80-100);
     `/script/<language>` template URL (:82-84 → `ScriptView.template` → `forScript`); and
     view-state restore — typeFactory `view.dart:181` → `ScriptView.fromParams` (:286-294).
     Apply the same branch after the script is loaded in each.
   - **Missing-handler guard (required, not existing behavior)**: today `ScriptView.init()` calls
     `ScriptHandler.forLanguage(script.language).codeEditorMode` unguarded (:373) and *throws* for
     an unregistered language — `addScriptView`'s `init().then(...)` has no catchError (:107-118,
     view hangs on the loader) and `/script/<id>` mislabels the error "Script not found" (:90-98).
     Add a guard that falls back to a plain-text CodeMirror mode when no handler is registered, so
     a Flow-less client opens the header+JSON body harmlessly (Run stays disabled by the existing
     invalid-language validation, :402-403).
4. **`FlowScriptMeta`** (`core/client/xamgle/lib/src/meta/flow_script_meta.dart`, `part` of
   `meta.dart`, registered in `initEntities()` **before** `ScriptMeta`):
   ```dart
   class FlowScriptMeta extends ScriptMeta {
     @override bool isApplicable(x) => x is Script && x.language == 'flow';
     @override String get typeName => 'Flow';   // display only; type stays $Script.TYPE
   }
   ```
   Overrides (each guarded, `Funcs.byNamespace('Flow')` / `Funcs.byName(...) != null`, with the
   Notebook-style `checkEnabled` message `'Install the Flow package first'`):
   - `init()`: `regCommand('Open Editor', isDefault: true, run: openInFlowEditor)` (double-click →
     visual editor; delegates to `addScriptView`, which the seam redirects), keep inherited
     `Run...`/`Share...`, drop `Debug...`.
   - `renderAccordion`: inherited panes, but replace the 'Script' source pane with a 'Flow' pane
     whose content is `Flow:flowScriptWidget` (info-panel call pattern:
     `f.apply(values: [x], processed: true)` → `Widget.root`); fall back to the source pane when
     the package is absent.
   - `renderPreview`: `Flow:flowScriptPreview` → `resultParamValue as View` (AppMeta pattern,
     keep `parentCall`); fallback `super.renderPreview`.
   - Icon: nothing to do when the package is installed — `Script.getImagePath()` resolves the
     handler's `meta.icon` against the package root. Package absent ⇒ no registered handler ⇒
     `getImagePath()` returns null ⇒ standard scroll icon (`script_meta.dart:64`). The
     `/images/entities/<language>.png` default is only reachable *with* a registered handler, so
     shipping a core image would change nothing — accept the scroll fallback.
5. Optionally add `static const String Flow = 'flow';` to `ScriptLanguage`
   (script_handler.dart:4-14) and use it instead of the literal.

**Acceptance**: double-click / URL / Edit / New all open the visual editor when Flow is installed,
CodeMirror otherwise; context panel shows the Flow pane; Browse preview renders the graph.
`dart analyze` clean; codegen re-run.

### Phase 3 — Core: Browse node, gallery, creation UX

1. **`FlowsView extends DataSourceCardView<Script>`** (`core/client/xamgle/lib/src/views/flows_view.dart`,
   modeled 1:1 on `ScriptsView`): `makeDataSource() => dapi.scripting` with
   **`permanentFilter = 'language = "flow"'`** — do NOT pre-apply `whereSmart` in
   `makeDataSource`: `DataSourceCardView` overwrites `textFilter` on every refresh
   (data_source_card_view.dart:729-732) and instead AND-combines `permanentFilter` into the search
   string (:158-162, :853-858; precedent `PropertySchemasView..permanentFilter` at view.dart:259).
   The server-side smart filter supports language (`getScripts` → `.whereSmart(text)`,
   scripting_service.dart:109; ScriptsView already ships a `'language = "R"'` filter example,
   scripts_view.dart:38-44). Plus `meta => new FlowScriptMeta()`, card renderer showing the flow
   summary line, `commands` = New Flow (opens template via `addTemplateScriptView('flow')`),
   `commandsButtonName = 'New'`.
2. **Browse tree**: `funcsNode.addItem('Flows', icon: ...)` next to Scripts in
   `browse_panel_tree.dart` (~:206); `nodePreviewFactories['Platform/Functions/Flows'] =
   (params) => new FlowsView(...)` + `urlNodesMapping['flows']` in `browse_panel.dart` constructor;
   `View.typeFactories`/`urlFactories` entries for `/flows` (the 5-point checklist in Insights §6).
   Gate: show the node only when flows exist or the package is installed
   (`Funcs.byNamespace('Flow').isNotEmpty` — cheap sync check), to avoid noise on bare servers.
3. **Exclude flows from the Scripts gallery** (`ScriptsView`:
   `permanentFilter = 'language != "flow"'`) — decide in review; alternative is to leave them
   listed in both.
4. `Tools | Scripting | Flow` "New" command appears automatically from `ScriptMeta.initCommands`'s
   per-handler loop; label confirmed to read **"New Flow Script..."**
   (`'New ${handler.friendlyName} Script...'`, friendlyName = capitalized language). Caveat:
   `initCommands` runs **once** after `initFuncs()` (`xplorer_init.dart:180-181`), and the runtime
   func-sync listeners skip the handler-registration branch — a Flow package published
   mid-session gets its language/commands/editor-seam only after a page reload. Acceptable v1;
   a live-sync fix is a Phase-4 candidate.

**Acceptance**: Browse shows Flows with cards; `/flows` URL works; New Flow → editor → Save →
appears in gallery without refresh (`AppEvents.ENTITY_ADDED` fires via normal save path).

### Phase 4 — Polish (each item independent; pick in review)

- **Space category**: split `Flows` out of `Scripts` in `ProjectMeta._getProjectStats` +
  `dsFactory` switch (project_meta.dart:126-154, 420-504) using a language-filtered source.
- **Read-only preview**: proper view-only mode in `FuncFlowView` (disable drag/connect/delete).
- **`.flow` files**: `fileViewer-flow` function so `.flow` exports in shares open in the editor
  (handler `extensions: flow` already claims the extension for script parsing).
- **Run pane**: default `_FuncMeta.renderRunSection` should be fine (params are ordinary); verify
  and only customize if the UX demands it.
- **Migration helper**: "Import .ffjson → entity" command in the Flows gallery.
- **Card/tooltip art**: mini graph thumbnail on cards (canvas render of node boxes — cheap, pure).
- **Startup-payload guard**: measure typical body sizes; if large, consider a core follow-up to
  strip bodies from `entity_jsons` startup sync for scripts above a size threshold (separate
  ticket — touches `cleanupFuncForJson`).
- **Server-side execution**: future ticket — a datlas-side flow handler is only possible once
  server can run package handlers (or Flow pre-compiles JS into a sibling script on save).
- **Live handler sync**: register package script handlers (and refresh `ScriptMeta.initCommands`)
  from the func-sync runtime listeners so a mid-session Flow install works without reload
  (func_sync.dart:84-97 currently skips the handler branch at :110-117).

## 7. Degradation matrix (package absent)

| Surface | Behavior without Flow package |
|---|---|
| Entity list/search/spaces/sharing | Fully functional (plain Script) |
| Params of a flow script | Correct — parsed from the annotation header by core `ScriptParser` |
| Icon | Default scroll icon (`Icons.fa('scroll')`) — a handler icon exists only when the package is installed |
| Double-click / Edit | With the Phase-2 missing-handler guard: plain-text CodeMirror view of header+JSON, Run disabled. (Without the guard, core *today* throws on open — the guard is mandatory Phase-2 work) |
| Run (any path) | `'Handler for flow is not registered'` at run; meta commands show `checkEnabled` message `'Install the Flow package first'` |
| Context panel Flow pane / Browse preview | Fall back to inherited Script panes / run-form preview |
| Browse "Flows" node | Hidden (gated on namespace presence) |

## 8. Testing plan

- **Flow package** (`grok test --host localhost`): format round-trip; parser param derivation
  (incl. qualifiers, defaults, captions); handler execution end-to-end (scalar + dataframe outputs;
  error propagation); save/open entity round-trip preserving layout/labels/annotations.
- **ApiTests + ApiSamples** (rule: new public surface needs both): saving/finding/running a flow
  script via `grok.dapi.scripts` + `grok.functions.call`; sample under
  `ApiSamples/scripts/` showing programmatic flow-entity creation. Clean up created entities in
  `tearDownAll`.
- **Core** (Dart): unit-test the editor-seam branch (handler with/without `editorFunc`, func
  present/absent); `dart analyze` on ddt/grok_shared/xamgle after codegen.
- **E2E/UX pass**: the acceptance stories in §3, on a stand with and without the Flow package
  (use the `user`/`ux-designer` agents per the repo's agent-dev loop).

## 9. Risks and open questions (for review)

1. **256 K body ceiling** — very large flows could hit it. Mitigation: compact JSON (drop nulls,
   round positions); error with a clear message on save. Long-term: widen column or external blob. (user answer: "Not a problem, we will never realistically hit 256K")
2. **Startup sync weight** — every flow body ships to every client at login (Insights §3). Fine at
   small scale; needs the Phase-4 guard if flows proliferate.  (user answer: realistically, each flow is on average 10kb, so its fine, we can address it later)
3. **Header/JSON drift** — mitigated by single-writer invariant (§5) + parser preferring JSON;
   a `grok check`-style validation on save could assert consistency.
4. **Mid-session package publish** — handler/commands/editor-seam only pick up at login
   (Insights §9.9); users who install Flow live must reload. Phase-4 candidate: re-run the
   handler-registration branch and `initCommands` from the func-sync runtime listeners. - user answer: not a realistic problem
5. **Editor race** — entity-open must await function-registry readiness (Phase 1.2); the current
   100 ms timers are not deterministic on slow catalogs. user answer: needs checking using tests
6. **Naming** — "Flow" vs "Workflow" in UI labels; the auto-generated menu item reads
   "New Flow Script..." (Phase 3.4). Decide once. user answer: New Flow is fine
7. **Debug command** — inherited `Debug...` makes no sense for flows (dropped in Phase 2.4); the
   editor's own Run/Debug covers it. user answer: correct, not needed there
8. **Do we exclude flows from the Scripts gallery** (Phase 3.3)? Recommend yes. user answer: yes
9. **JIRA** — core changes need a GROK ticket (`GROK-NNNNN: Scripting: per-language editor seam`
   + `Flow: first-class entity`); commit format per repo rules. user answer: VERY IMPORTANT: do not auto commit anything. no need for jira ticket for this one. do not auto commit anything. I will review and commit myself 
