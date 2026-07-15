---
feature: browse
target_layer: playwright
coverage_type: smoke
priority: p0
realizes_atlas: []
realizes: []
realized_as:
  - apps_matrix.test.ts
related_bugs: [GROK-20254]
---

# Browse — Apps matrix (manual test cases)

Manual test cases for walking through **all applications** in `Browse > Apps`. Each application
is opened by **clicking its Browse tree node (by name)**; applications that expand into a
data subtree are traversed in depth — each sub-node opens its own view, with which a
minimal interaction is performed. The main check is the **absence of any errors** (console,
balloons, Context Panel).

Automation: `apps_matrix.test.ts` (same directory). The case structure maps 1:1 to the `APPS` array and
the traversal logic in the code.

---

## 1. Environment and common preconditions

- **Instance:** `https://dev.datagrok.ai/` under a logged-in user with access to demo/data.
- **Matrix source:** the live `Apps` tree captured from dev (see `browse_tests_context.md`).
- **Open:** the **Browse** side panel, with the **Apps** node expanded.
- Navigation is strictly **by node name** (the stable attribute `name="tree-Apps---…"`).

## 2. Scope

**All top-level nodes under `Apps` are covered, EXCEPT**:
- `Demo` — covered by `demo_apps.test.ts`;
- `Compute > Model Hub` — covered by `modelhub.test.ts`;
- `Alation` — the application is deprecated (remains only on dev), excluded from the matrix.

## 3. How an application opens (platform behavior)

- Clicking an application leaf invokes its `@app` function → a new view opens (view name = node
  name), and the Browse panel may switch to the Toolbox.
- Group applications build their subtree **lazily** (`@appTreeBrowser`, sometimes querying the DB
  or a docker container) — child nodes appear with a delay and **incrementally** (streaming).
- Some applications depend on a **docker container** (MolTrack, Admetica, Retrosynthesis,
  Docking/Boltz via `bio`, Preclinical Case) — if the container is not up, the view does not function.
- Some applications are **external connectors** (Benchling, CDD Vault, Revvity Signals, Chemspace) —
  they depend on a third-party SaaS that is available only on dev.

## 4. Coverage depth ("sample")

To keep the suite within the CI timeout (1 worker) and prevent combinatorial blowup (one
`Preclinical Case` = 216 nodes):

- **Leaf application** (`Flow`, `PopPK`, `Tutorials`) → open, interact, verify.
- **Group application** → open **ALL** top-level child nodes (tools / instances), and
  for each expandable child walk **ONE representative chain** (the first child at each
  level, down to a leaf). This way a single data instance is opened "in depth" without traversing all
  same-type data leaves.
- `Tutorials` — open only (`runOnly`); the automated test for completing tutorials is separate, later.

## 5. Common test case (for each node)

**ID:** `Browse-AppsMatrix-<App>[-<sub-nodes>]`
**Type:** smoke.

**Steps (action → expected reaction):**
1. Expand `Apps` → applications appear.
2. Expand the application (by name) → its child nodes appear (wait for the tree to finish building).
3. Click a node → its view opens.
4. **Poke around the view:** click the view body; if there is a grid — click a cell; hover over the
   main viewer.
5. Open the **Context Panel** (F4) + **Expand all** → all property panels render.

**Final check (for all nodes):**
- ✅ no uncaught JS errors (`pageerror`);
- ✅ no `console.error`;
- ✅ no error balloons;
- ✅ when the Context Panel is expanded, no panel crashes.

## 6. Application matrix (11)

| App | Type | Tags | Notes |
|-----|-----|------|---------|
| Bio | group | heavy | Boltz-1 / Consensus-Pharmacophore / Docking (docker `bio`) |
| Chem | group | heavy | tools + connectors (Benchling, CDD Vault, Revvity, Chemspace) + MolTrack (docker) + Hit Triage/Hit Design |
| Clinical Case | group | heavy | 2 studies (CDISCPILOT01, nida-ctn-0001) × ~19 views + Import study |
| Compute | group | heavy | Diff-Studio, KNIME, Simulations, NCA (Model Hub excluded) |
| Flow | leaf | | |
| Misc | group | | Excalidraw, MetabolicGraph |
| Peptides | group | heavy | Oligo-Toolkit, PeptiHit, PepTriage, PolyTool, Monomer-Libraries |
| Plates | group | | Create / Search plates/wells/analyses / Templates |
| PopPK | leaf | heavy | |
| Preclinical Case | group | heavy | 7 studies × domains (docker `preclinicalcase`) |
| Tutorials | leaf | runOnly | open only |

## 7. Error handling (what is NOT considered a bug)

The suite honestly turns red on real Browse/UI defects, but **tolerates** known non-defects —
scoped to the node on which the error occurred:

1. **Infrastructure (always):** connection drops, `Container … not started`, raw backend HTTP
   responses `4xx/5xx`. Signatures are narrow so as not to mask a real bug that happens to have a
   number in its text.
2. **Docker-dependent apps:** on a docker error, the **live container status** is checked
   (`grok.dapi.docker.dockerContainers`). Container `started` + error = **real bug (fail)**;
   container not `started` (`starting`/`error`/`stopped`/absent) = infrastructure, **tolerated**.
3. **Connectors** (Benchling, CDD Vault, Revvity, Chemspace): depend on an external SaaS available
   only on dev → all their errors are tolerated.
4. **Interaction noise:** `ResizeObserver: parameter 1 is not of type 'Element'` and React
   error-boundary echoes — provoked by our `pokeView` on a mounting React editor; when opened
   calmly in isolation they do not occur → not a navigation defect.
5. **Async-bleed:** a docker-down / connector sub-app spews balloons that "stick" to the neighboring
   node → after leaving such a sub-app a `closeAll` is performed to silence it.

Each toleration is logged with the line `ⓘ <node>: N error(s) tolerated — <reason>`.

## 8. Portability (other stands / CI)

The matrix was captured on dev (the full set). On a minimal stand the application package is not
installed → the application node is absent → the entire app is `test.skip()`-ed; absent sub-nodes are
skipped individually. There is no hard-coded allowlist — the same file is correct on both dev and a
trimmed-down stand. Protection against false-green: if the application is present (a group) but not a
single node opened — this is a subtree loading failure (fail), not "vacuous success".

## 9. Run results on dev (triage)

Full run of all 11 tests on `dev.datagrok.ai` (2026-06-23): **11 passed / 0 failed**.

- 🐞 **GROK-20254** — Hit Triage / Hit Design threw a repeating
  `Cannot read properties of null (reading 'id')` on a timer while the view was open. Found by this
  suite, confirmed in isolation, ticket filed; as of the final run it is **already fixed on dev**
  (the test is green).
- 🔧 **Plates / Search plates** — `RangeSlider.setValues` crashed on `grid.columns.add` — turned out
  to be a **test race** (`closeAll` destroyed the view while its debounced `search().then()` was
  still in flight). Not a platform bug. Fix: do not `closeAll` between sub-nodes + wait for
  `networkidle`+settle.
- 🔧 **ResizeObserver `<Editor>`** (MolTrack Register / Reaction Enumerator) — provoked by
  `pokeView`; clean in isolation → noise, added to the ignore list.
- ⓘ **MolTrack** — `Container not started` (the `moltrack` container = `starting`): infrastructure,
  tolerated by the live docker status.

## 10. Notes

- The application list was captured automatically; when the set changes on the server, regenerate the
  matrix (see `browse_tests_context.md`).
- `heavy` — the application launches server/docker computations and opens more slowly (increased test
  timeout).
- Deep same-type data leaves (datasets of a clinical/preclinical study domain) are **not** traversed
  in full — that is exactly what "sample" means.
