# Dashboard-First Flow — Making Flow the Primary Dashboard-Creation Tool

**Companion to** [COMPETITIVE-ANALYSIS.md](./COMPETITIVE-ANALYSIS.md) (feature parity vs KNIME/Spotfire)
and [UX-FOR-SCIENTISTS.md](./UX-FOR-SCIENTISTS.md) (the "Dr. Maria" approachability lens). Those two
documents made Flow a better *pipeline editor*. This one argues for a strategic reorientation:
**the pipeline is the recipe — the dashboard is the dish.** A scientist doesn't want a flow; she
wants the thing a flow produces: a living, shareable, refreshable picture of her data. Today Flow
stops one step short of that, and that last step is where all the audience is.

**Date:** 2026-07-10 · grounded in Flow's source, the Datagrok core Dart code
(`GROK-CORE/reddata`), the js-api, the platform help, and the Spotfire/KNIME publishing models.

---

## 1. The thesis

> **Every flow should end in a dashboard, and every dashboard should be able to explain itself as a flow.**

Three observations drive this:

1. **Consumers outnumber authors 10:1.** In every Spotfire/KNIME deployment, a handful of people
   build pipelines and everyone else *consumes* the result — a page of linked charts they can
   filter, with data that refreshes. A pipeline tool whose output is "some tables appear in the
   workspace" has no consumer story, so it caps out at the author population.
2. **Spotfire won scientists with the inversion.** In Spotfire the *analysis file* (pages of
   visuals) is the artifact; the data pipeline ("data canvas") sits **behind** it, reachable when
   you ask "where did this table come from?". KNIME learned the same lesson the other way around —
   the pipeline is primary, but its killer deployment feature is **Data Apps**: a component's
   widgets + views become a browser page for non-KNIME users. Both converged on: *author in the
   graph, deliver a page.*
3. **Datagrok already has the entire dashboard machinery — Flow just doesn't drive it.** A
   dashboard here is not a new thing to build; it's a `Project` with `isDashboard = true`, saved
   views, layouts, data sync, presentation mode, sharing, and iframe embedding — all shipped, all
   reachable from JS (§3). What's missing is *only* the bridge from the Flow canvas to that
   machinery, and a UI that treats the dashboard as the goal instead of the script.

Flow's unique advantage carries over intact: the dashboard produced by a flow is a **glass-box
dashboard**. Unlike a Spotfire DXP someone hand-assembled, every table in it has a visible,
re-runnable, *editable* recipe. "See how this was made" opens the canvas. No other tool in this
space can honestly offer that.

---

## 2. Where Flow stands today (the last-mile audit)

What happens today when someone builds and shares an analysis in Flow:

| Stage | Current behavior | Verdict |
|---|---|---|
| Author adds viewers | `Viewers/*` nodes exist; preview renders them live in the output panel ([value-inspector.ts](../src/execution/value-inspector.ts)) | ✅ good bones |
| Author arranges results | The output preview shows one node's outputs at a time (`ui.splitH` for multi-output); there is **no place to compose** several nodes' outputs into one page | ❌ no canvas→page step |
| Author saves | Flow saves as a `DG.Script` (language `flow`) into a Space; Save button state, save-as, spaces all work ([entity/](../src/entity/)) | ✅ shipped |
| Consumer runs the saved flow | The script handler executes the emitted JS; consumer gets **output tables**. Viewer nodes emit `let v = await table.plot.fromType(...)` — a **detached viewer object that lands nowhere** unless declared an output ([script-emitter.ts:424](../src/compiler/script-emitter.ts#L424)). `Add Table View` opens a bare grid ([script-emitter.ts:301](../src/compiler/script-emitter.ts#L301)) | ❌ the composed experience is lost at the door |
| Consumer interacts | Nothing: no controls, no filters page, no layout, no refresh affordance | ❌ |
| Consumer asks "how was this made?" | Must find the script entity and open the Flow editor | ⚠️ exists but not discoverable from the result |

**The audit in one sentence:** everything the author *sees* while building (live viewers, previews,
arrangement) evaporates at publish time; the consumer receives the raw ingredients, not the dish.

---

## 3. The machinery the platform already gives us

Everything below is shipped, verified in core/js-api/help — this section is the "we don't have to
build the hard part" inventory.

| Capability | Mechanism | Where |
|---|---|---|
| **Dashboard entity** | A `Project` with `isDashboard = true` (`Project.dashboard()` factory); holds relations to `TableInfo`s, `ViewInfo`s, `ViewLayout`s + a `ProjectLayout` (`{views: [...]}`); dedicated `CREATE_DASHBOARD` / `BROWSE_DASHBOARDS` permissions | `grok_shared/lib/src/project.dart`, `privileges.dart`; JS: `DG.Project.create()`, `project.isDashboard`, `dapi.projects.save()` |
| **Compose a page** | `TableView` + `addViewer(type, options)` + `view.dockManager` (fill/left/right/up/down + ratio) — programmatic multi-viewer arrangement | js-api `src/views/view.ts`, `src/docking.d.ts`; Dart `xamgle/lib/src/features/docking.dart` |
| **Persist the arrangement** | `tv.saveLayout()` → `ViewLayout` (`toJson`/`fromJson`, `dapi.layouts`); `project.addTableView(df)` + `saveRelations` capture views into the project | js-api `src/entities/view-layout.ts`, `src/entities/project.ts` |
| **Refreshable data ("data sync")** | Every table produced by a function carries its **generation script**; on project upload, toggling **Data sync** stores the recipe instead of a snapshot — the function **re-executes each time the project opens**; manual refresh + parameter editing under **Toolbox → Source** | help: `datagrok/navigation/basic-tasks` §Dynamic data, `access/databases` §dynamic dashboards; Dart: `TableInfo.metaParams[Tags.DataSync]` |
| **Parameterized dashboards** | Parameterized queries expose their inputs to dashboard *consumers*, who tweak them in place — the platform's existing pattern for interactive dashboards | help: `access/databases.md` §Creating dynamic dashboards |
| **Declared output viewers** | `//output: dataframe result {viewer: 'Scatter plot(x: ..., y: ...)'}` — script outputs can declare how they render; the function results view honors it; RichFunctionView (compute-utils) builds a full input-form + output-tabs page from annotations alone | Dart `scripting.dart`, `function_view.dart`; `libraries/compute-utils/.../rich-function-view.ts` |
| **Consumer mode** | **Presentation mode** flag on project upload: no toolbox, ribbons, or panels — "focus on the data and visualizations you have created" | help: basic-tasks §Save and share |
| **Distribution** | Sharing with users/groups (+email notification), `#dashboard`-tagged functions surface on the platform home, views embed as **interactive iframes** in external sites | help: project sharing / embed; `DG.FUNC_TYPES.DASHBOARD` |
| **Scheduled refresh** | `Func.recurrence` + `nextRunTime` (server-side scheduled runs), `meta.cache` + `meta.cache.invalidateOn: <cron>` | Dart `func_call_cache.dart`, `func_roles.dart` |

**The keystone fact:** a flow already saves as a real `DG.Script`. That means a flow can *be* the
"generation script" of every table in a dashboard project — which makes **data sync, manual
refresh, parameter editing, and scheduling work for Flow dashboards with zero new server
machinery.** The bridge is thinner than it looks.

---

## 4. What the incumbents teach about the dashboard step

- **Spotfire** — the artifact is pages of linked visuals; authoring is *chart-first* ("recommended
  visualizations" as soon as data loads; drag a column onto a page). The pipeline is subordinate
  and discoverable from the artifact. Consumers get the web player: filter, mark, drill — never the
  authoring canvas. Data refreshes on schedule on the server.
- **KNIME** — authoring is *graph-first* (like Flow), and the deployment answer is **components
  with composite views → Data Apps**: widget nodes (dropdown, slider, text) + view nodes compose a
  page; the Hub serves it to browser users who never see a node. The lesson: **input nodes are
  UI controls in disguise.** Flow already has the input-node half of a data app.
- **Both**: publishing is a *single verb* on the primary toolbar, the consumer artifact hides the
  machinery by default, and "refresh" is a property of the artifact, not a re-authoring act.

The synthesis for Flow: keep KNIME's graph-first authoring (it's our glass box), adopt Spotfire's
artifact-first *delivery*, and use Datagrok projects as the delivery vehicle both of them had to
build proprietary servers for.

---

## 5. The gaps (what actually stands between Flow and "primary dashboard tool")

1. **G1 — No composition surface.** Viewers/outputs preview one node at a time; there is no "page"
   the author arranges. (The output panel is a debugger, not a designer.)
2. **G2 — Publish verb missing.** "Save" persists the recipe; nothing creates the consumer
   artifact (project + views + layout + sharing).
3. **G3 — The emitted script drops the experience.** Detached viewers; no view/layout emission;
   `Add Table View` is the only (bare) window into results.
4. **G4 — No consumer surface.** Opening a saved flow lands you in the *editor* (auto-pinned
   preview); there's no read-only, presentation-mode "just the dashboard" route.
5. **G5 — Inputs aren't controls.** Flow's Input nodes only become script parameters; they never
   surface as dashboard-side widgets the consumer can tweak (the KNIME Data-App half we're missing).
6. **G6 — Refresh story invisible.** Data sync/scheduling exist platform-side but nothing in Flow
   sets `Tags.DataSync`, offers "refresh on open," or schedules a re-run.
7. **G7 — UI is still script-oriented at the finish line.** The ribbon's terminal verbs are Run /
   Save / View Script. A scientist's terminal verb is **"Share this dashboard."**

---

## 6. Proposal — the Dashboard-First reorientation

Each item: **what / why / how** (anchored to real APIs), effort (S/M/L/XL), impact.

### D1 — The Dashboard panel: a composition surface inside the editor ⭐ keystone

**What:** Add a **Dashboard tab** next to the output preview (or a full-height right-side mode):
a real `DG.TableView`-backed page where *every* viewer node and *every* output table of the flow
appears as a docked tile the author can drag, resize, and remove. It updates live after each run —
the same `funcflow.exec` events that feed the preview feed the dashboard.

**Why:** This closes G1 and changes what authors *think they're making*. The preview answers
"did this node work?"; the dashboard panel answers "what will my colleagues see?". Every
downstream feature (publish, consumer mode, refresh) needs this artifact to exist first.

**How:** Maintain a `DashboardDraft` model: `tiles: [{nodeId, outputKey, kind: grid|viewer|text,
dock: {style, ratio, target}}]`, persisted in the `.ffjson` (new `dashboard` section in
[flow-schema.ts](../src/serialization/flow-schema.ts)). Render with `TableView.create` +
`addViewer` + `view.dockManager.dock(el, style, ratio)`; capture arrangement via
`tv.saveLayout().toJson()`. Default population: every viewer node → its viewer; every Table
Output → a grid; author curates from there. The instrumented runner already clones/captures every
output (`__ff_stash` registry) — the dashboard reads the same registry the preview does.
**Effort: L · Impact: XL.**

### D2 — "Publish dashboard" as the primary ribbon verb ⭐

**What:** One button — **Publish** — that turns the draft into a platform dashboard: creates
`Project.dashboard()`, adds the output tables (`project.addTableView(df)` / `addLink`), applies the
saved `ViewLayout`, links the flow script entity itself into the project, sets **presentation
mode**, uploads via `dapi.projects.save(project, {saveRelations: true})`, then opens the share
dialog (groups, notify). Re-publishing updates the same project (the bound-entity pattern we
already use for scripts).

**Why:** G2/G7. This is the single verb that makes Flow a *dashboard tool* rather than a script
tool. It also rides every governance feature (permissions, browse-dashboards, project gallery)
for free.

**How:** All APIs exist (§3). The Flow side is a `publishDashboard()` that walks `DashboardDraft`,
plus a small "Publish" dialog (name, description, space, data sync toggle, share-with). Follow the
`saveAsDialog`/space-picker patterns already in [funcflow-view.ts](../src/funcflow-view.ts).
**Effort: M · Impact: XL.**

### D3 — Live data: flow-as-generation-script + refresh ⭐

**What:** Published dashboards default to **Data sync**: tables are stored as their recipe (the
flow), not a snapshot — opening the dashboard re-runs the flow; a **Refresh** button re-runs on
demand; optionally a schedule ("refresh nightly") via `Func.recurrence` / `meta.cache.invalidateOn`.

**Why:** A dashboard that shows last Tuesday's data is a screenshot. "It's always current" is THE
reason query-backed dynamic dashboards exist on the platform; flows must be first-class citizens
of the same mechanism. This is also our answer to Spotfire's scheduled server refreshes.

**How:** The flow is already a `DG.Script`; ensure tables returned from running it carry the
generation-script provenance (they do when produced through a `FuncCall` — verify for the
flow-script handler path, see §8-Q1), then set `Tags.DataSync` on the project's TableInfos at
publish. Scheduling = set `recurrence` on the saved script entity (needs a small UI in the Publish
dialog). **Effort: M (given D2) · Impact: XL.**

### D4 — Input nodes become dashboard controls (the Data-App half)

**What:** Any flow **Input node** (string/number/choice/column/file/datetime) surfaces in the
published dashboard as an editable control in a slim sidebar — change a value → the flow re-runs →
tiles update. Exactly what parameterized query dashboards already do, generalized to any flow.

**Why:** G5. This is KNIME Data Apps and Spotfire document properties in one move, and it falls
out of infrastructure we have: the flow script already declares `//input:` lines with captions and
qualifiers; the platform already renders input forms for functions (`FuncMeta.renderRunSection`,
RichFunctionView proves the full pattern).

**How:** Two implementation candidates — (a) rely on the platform's query-parameter dashboard UI
(**Toolbox → Source**) once data sync points at the flow function (cheapest, weakest layout
control); (b) render our own control strip in the published view from the script's `inputParams`,
calling `grok.functions.call(flowNqName, params)` and re-binding tiles (more work, full control).
Start with (a), graduate to (b). **Effort: M–L · Impact: XL.**

### D5 — Dashboard-grade content nodes

**What:** Close the gap between "viewers" and "a page":
1. **Text/Markdown card** node (title, commentary, methods note) — dashboards need words;
2. **Filters tile** (the platform filter group as a dockable tile — instant interactivity since
   all viewers of a DataFrame are already linked/brushed natively);
3. **Chart-first authoring**: in the Dashboard panel, an "**+ Add chart**" button on any table tile
   (gallery of viewer types, then the existing gear/`ff-viewer-gear` property flow) — so making a
   dashboard doesn't require knowing viewer *nodes* exist;
4. Optional later: image/logo card, KPI/"big number" card (a one-cell aggregation → styled div).

**Why:** These four cover ~90% of what real dashboards contain beyond raw charts (see any Spotfire
page: header text, filter rail, charts). Linked brushing we get free — same-DataFrame viewers in a
TableView already select/filter together, which is *better* than what a casual KNIME data app gets.

**How:** Text card = trivial node + tile kind. Filters tile = `tv.getFiltersGroup()` docked.
Add-chart = UI over `addViewer` + existing viewer-node creation (write the node back into the
graph so the glass box stays truthful — the dashboard must remain a *view of the flow*, never a
divergent copy). **Effort: M · Impact: L.**

### D6 — The consumer route: dashboard first, canvas on request

**What:** Opening a published Flow dashboard opens **the project in presentation mode** — tiles,
controls, refresh; no canvas, no toolbox. One affordance in the corner: **"See how this was
made"** → opens the flow editor on the bound script (our existing preview/auto-pin path). The
entity handler's card ((flow-entity-handler.ts)) gains a "Open dashboard" default action, demoting
"Open editor" to secondary for consumers.

**Why:** G4. This is the Spotfire inversion executed with Datagrok primitives — and it turns our
glass box into a consumer feature: trust-by-inspection is one click away but never in the way.

**How:** Project open + presentation mode is platform behavior; Flow's work is routing (the
`/script/<id>` path already routes to the editor; add `/dashboard` or open-project) and the
entity-handler card. **Effort: S–M · Impact: L (strategically XL).**

### D7 — Auto-draft the dashboard (recommendations)

**What:** After the first successful run, if the author has no dashboard yet, propose one:
every viewer node placed, every output table gridded, a sensible dock split, a title card from the
flow name/description. A toast: *"Your dashboard draft is ready — arrange and publish."*

**Why:** Spotfire's "recommended visualizations" is its most-loved onboarding feature; our
UX-FOR-SCIENTISTS review (F9 "the canvas is never inert") says defaults must produce something
visible. The auto-draft turns *any* working flow into a dashboard in zero clicks — which is
exactly the "easiest thing in the UI" requirement.

**How:** A pure function `draftFromFlow(flow, execRegistry) → DashboardDraft` + the D1 renderer.
Heuristics: ≤2 tables → side-by-side grids + their viewers stacked; >2 → tabs. **Effort: S–M
(given D1) · Impact: L.**

### D8 — Dashboard-outcome templates & guide chapter

**What:** Recast the Start panel templates as *dashboard outcomes*: "Assay QC dashboard from a
plate file", "Compound property dashboard from SMILES", "Two-table join dashboard with filters" —
each ends published-preview-ready, not just "nodes on canvas". Add a guide tutorial: *Build and
share your first dashboard* (reuses the guide engine + `untilSectionExpanded`-style gates).

**Why:** UX-FOR-SCIENTISTS U1/U2 said "never open empty; a win in the first session". Under the
dashboard-first thesis, the first-session win must *be a dashboard*, because that's the artifact
Maria can show her PI the same afternoon.

**How:** Template `.ffjson`s now carry the `dashboard` section (D1's schema). Guide content follows
the existing [guide-content.ts](../src/guide/guide-content.ts) patterns. **Effort: M · Impact: M.**

### D9 — Ribbon & language reorientation (the cheap, loud signal)

**What:** Make the ribbon read as a dashboard tool: **Run · Dashboard · Publish · Share** up
front; *View Script / Creation Script / Export .js / Debug* stay under Advanced (per U8, already
partially done). Rename the output panel tab "Preview" → keep; the new tab is "**Dashboard**".
Status bar gains "Published · synced 2h ago" when bound to a project.

**Why:** G7. Language *is* orientation. This is days of work that changes what the tool claims
to be — the highest signal-per-hour item in this document.

**How:** [funcflow-view.ts](../src/funcflow-view.ts) ribbon arrays + status bar; no engine work.
**Effort: S · Impact: M.**

### D10 — Distribution extras (once D2 exists)

**What:** (a) Register a published dashboard's opener as a `#dashboard`-tagged function so it
appears on the platform home/browse-dashboards surfaces; (b) surface the platform's **Embed…**
(iframe) for any published view; (c) a "Copy share link" one-liner after publish.

**Why:** Meeting consumers where they already look (home screen, wikis, ELNs). Iframe embedding
into an ELN/Confluence page is a real daily pattern for pharma teams.

**How:** All platform features; Flow adds buttons and a tiny function registration at publish
time. **Effort: S · Impact: M.**

---

## 7. Storyboard — Maria ships a dashboard

> Maria opens Flow, picks **"Compound property dashboard"** from the Start panel, and drops her
> SMILES file on the amber *Needs a file* node. **Run.** The nodes tick green — and the
> **Dashboard** tab lights up: a grid of her compounds, a scatter of logP vs MW, a histogram, and
> a title card, already arranged (D7). She drags the histogram under the scatter (D1), retitles
> the card, and adds a **Filters** tile (D5).
>
> She clicks **Publish** (D2): names it "AZ-114 series properties", leaves **Keep data current**
> on (D3), shares with *Chemists* (D10). Her PI gets a link; it opens as a clean page — charts,
> filters, a *Refresh* button, and her flow's `logP cutoff` input as an editable control in the
> sidebar (D4). He drags the cutoff to 4, the page recomputes, he filters to one scaffold — the
> scatter and grid stay in sync (free brushing).
>
> Six weeks later a new chemist asks how logP was computed. He clicks **"See how this was made"**
> (D6) and the recipe unfolds — every step inspectable, nothing to take on faith. He clones the
> flow, swaps the input file, and publishes his own series' dashboard in ten minutes.

Every beat maps to a D-item; none requires new platform machinery.

---

## 8. Open questions to verify before building

- **Q1 — Provenance through the flow-script handler.** Confirm a table produced by running the
  saved flow (via our script handler) carries a generation-script/FuncCall reference the Data-sync
  mechanism accepts (`TableInfo.metaParams[Tags.DataSync]`). If the handler path drops provenance,
  publish must stamp it explicitly.
- **Q2 — Layout fidelity from JS.** `tv.saveLayout()` round-trip fidelity for docked *non-grid*
  viewers created programmatically (`addViewer` + dockManager) — spot-check with 3–4 viewer types.
- **Q3 — Multi-table pages.** One `TableView` binds one primary DataFrame; a dashboard over 2+
  output tables needs either multiple views in the project (`ProjectLayout.views[]` — supported)
  or in-view linked tables. Decide the tile→view mapping rule (proposal: one view per output
  table; tiles group by their source table; project stores all views).
- **Q4 — Re-publish semantics.** Updating a live dashboard someone has open / has bookmarked:
  reuse the project entity (stable URL) and replace relations — verify `dapi.projects.save`
  behaves as an upsert with `saveRelations`.
- **Q5 — Control strip re-runs (D4b).** Debounce + cancellation for consumer-triggered re-runs of
  expensive flows; consider `meta.cache` on the flow function so identical parameter sets are
  served from cache.
- **Q6 — Permissions defaults.** Publishing should request `CREATE_DASHBOARD`; sharing defaults
  (private / space-inherited) need a decision with the platform team.

---

## 9. Sequencing

```
Phase 1 — the artifact exists (the new keystone)
  D1 Dashboard panel (draft model + TableView/dockManager renderer)
  D7 auto-draft after first run
  D9 ribbon/language reorientation                       ← ship early, it's cheap

Phase 2 — the artifact ships
  D2 Publish → project (+ share dialog)
  D3 data sync + refresh (+ Q1 verification first)
  D6 consumer route + entity-handler default action

Phase 3 — the artifact is alive
  D4 input nodes as dashboard controls (a → b)
  D5 text card · filters tile · add-chart gallery
  D10 home-screen registration · embed · share link

Phase 4 — the funnel
  D8 dashboard-outcome templates + guide chapter
```

Dependency spine: **D1 → D2 → {D3, D4, D6}**; D7/D9 ride Phase 1; D5/D8/D10 are parallelizable
fill. Phase 1+2 alone already change the product's category: from "visual script editor" to
"dashboard designer with a glass-box recipe".

---

## 10. The pitch, one paragraph

Datagrok already owns the hard parts of a dashboard platform — linked viewers, layouts, dynamic
data sync, projects, permissions, embedding. KNIME and Spotfire each built half of the ideal:
KNIME the honest graph, Spotfire the consumable page. Flow can be the first tool where **the page
and the graph are the same object seen from two sides** — scientists author by wiring what they
understand, publish with one verb, and every chart a colleague ever looks at can answer the
question no other dashboard can: *"show me exactly how you made this."*
