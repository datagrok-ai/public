---
feature: browse
target_layer: playwright
coverage_type: smoke
priority: p0
realizes: []
realized_as:
  - demo_apps.test.ts
related_bugs: [GROK-14320, GROK-18050]
---

# Browse — Apps > Demo (manual test cases)

Manual test cases for viewing **all demo applications** in `Browse > Apps > Demo`.
Each demo is opened by **clicking the Browse tree node (by name)**, after which a minimal
interaction is performed in the view. The main check is the **absence of any errors**.

Automation: `demo_apps.test.ts` (same directory). The cases are structured for a direct port
into Playwright — the demo matrix below matches the `DEMOS` array in the code 1:1.

---

## 1. Environment and common preconditions

- **Instance:** `https://dev.datagrok.ai/` under a logged-in user with permissions for demo data.
- **Tree source:** functions with `meta.demoPath` (format `Category | Subcategory | Demo`),
  captured live from dev — **62 demos** across 6 categories. See section 4 for the list.
- **Open:** the **Browse** side panel (icon on the Sidebar), with the **Apps** node expanded.
- Navigation is done strictly **by node name**, not by index/order.

## 2. How a demo opens (platform behavior)

- Clicking a demo leaf invokes the demo function (`func.apply()`), which opens a new view and
  makes it current. **View name = the last path segment** (e.g. `Scatter Plot`).
- A secondary data view (`Table` / `Grid`) with the source data usually opens alongside it.
- The view header shows the breadcrumbs `Home / Demo / … / <Demo>`.
- Opening the next demo closes the previous demo views (they are marked `temp.demoApp`).
- Some demos are **dashboards** (`isDemoDashboard`): they open as a saved layout.
- **Heavy** demos (Compute, Docking, Admetica, Retrosynthesis, Databases, Map…) launch
  server/docker computations — they take longer to open and may depend on backend availability on dev.

## 3. Common test case (applies to every demo)

**ID:** `Browse-DemoApps-<Category>-<Demo>`
**Type:** smoke (for all); regression — for demos with `ref` (section 4).

**Preconditions:** section 1.

**Steps (action → expected reaction):**

1. In the Browse tree, expand `Apps` → reaction: child nodes appear, including `Demo`.
2. Expand `Apps > Demo` → categories appear (Cheminformatics, Bioinformatics, Data Access,
   Visualization, Compute, Curves).
3. Expand the category (and subcategory if present, e.g. `Visualization > General`) **by name**
   → demo leaves appear.
4. Click the demo leaf **by name** → a view opens with a title = the demo name;
   the header shows the breadcrumbs `Home / Demo / … / <Demo>`.
5. **Interact within the view:** click on the view body; if there is a table — click a grid cell
   (the current row changes); hover over the main viewer → reaction: highlight/tooltip,
   redraw without failures.
6. Open the **Context Panel** (F4) and click **Expand all** → all property panels render.

**Final check (for all demos):**

- ✅ a view with the demo name opened;
- ✅ no uncaught JS errors (`pageerror`);
- ✅ no `console.error` (except noise explicitly added to the ignore list — see `helpers.ts`);
- ✅ no error balloons;
- ✅ when the Context Panel is expanded, no panel fails (a panel error would surface in
  console/pageerror).

**Postconditions / cleanup:** return to `Home` (closes the demo preview); nothing is created on
the server — no server-side cleanup is needed.

## 4. Demo matrix (62)

Columns: **Path** — path under `Apps > Demo`; **View** — the name of the opened view; **Tags** —
`dashboard` (layout), `heavy` (server/docker), `skip` (platform marks `demoSkip`); **Ref** —
ticket/flag.

### Cheminformatics (10)

| Path | View | Tags | Ref |
|------|------|------|-----|
| Cheminformatics / Med Chem | Med Chem | | |
| Cheminformatics / Chemical Space | Chemical Space | skip | GROK-14320 |
| Cheminformatics / Molecule Activity Cliffs | Molecule Activity Cliffs | skip, dashboard | GROK-14320 |
| Cheminformatics / R-Group Analysis | R-Group Analysis | skip, dashboard | GROK-14320 |
| Cheminformatics / Matched Molecular Pairs | Matched Molecular Pairs | | |
| Cheminformatics / Similarity & Diversity Search | Similarity & Diversity Search | | |
| Cheminformatics / Scaffold Tree | Scaffold Tree | | |
| Cheminformatics / Database Queries | Database Queries | heavy | |
| Cheminformatics / Admetica | Admetica | heavy | |
| Cheminformatics / Retrosynthesis | Retrosynthesis | heavy | |

### Bioinformatics (10)

| Path | View | Tags | Ref |
|------|------|------|-----|
| Bioinformatics / Peptide SAR | Peptide SAR | dashboard, heavy | |
| Bioinformatics / Sequence Activity Cliffs | Sequence Activity Cliffs | heavy | |
| Bioinformatics / siRNA | siRNA | skip, heavy | demoSkip |
| Bioinformatics / Antibodies | Antibodies | heavy | |
| Bioinformatics / Sequence Space | Sequence Space | dashboard, heavy | |
| Bioinformatics / Similarity, Diversity | Similarity, Diversity | heavy | |
| Bioinformatics / Atomic Level | Atomic Level | skip, heavy | demoSkip |
| Bioinformatics / Docking | Docking | heavy | |
| Bioinformatics / Docking Conformations | Docking Conformations | heavy | |
| Bioinformatics / Proteins | Proteins | heavy | |

### Data Access (3)

| Path | View | Tags | Ref |
|------|------|------|-----|
| Data Access / Table Linking | Table Linking | | |
| Data Access / Files | Files | | |
| Data Access / Databases | Databases | heavy | |

### Visualization (33)

| Path | View | Tags | Ref |
|------|------|------|-----|
| Visualization / Data Flow and Hierarchy / Network Diagram | Network Diagram | | |
| Visualization / Data Flow and Hierarchy / Tree | Tree | | |
| Visualization / Data Flow and Hierarchy / Tree Map | Tree Map | heavy | |
| Visualization / Data Separation / Trellis Plot | Trellis Plot | | |
| Visualization / Data Separation / Matrix Plot | Matrix Plot | | |
| Visualization / General / Scatter Plot | Scatter Plot | | |
| Visualization / General / Bar Chart | Bar Chart | | |
| Visualization / General / Line Chart | Line Chart | | |
| Visualization / General / Histogram | Histogram | | |
| Visualization / General / Pie Chart | Pie Chart | | |
| Visualization / General / 3D Scatter Plot | 3D Scatter Plot | | |
| Visualization / General / Tile Viewer | Tile Viewer | | |
| Visualization / General / Density Plot | Density Plot | | |
| Visualization / General / Filters | Filters | | |
| Visualization / General / Heatmap | Heatmap | | |
| Visualization / General / Markup | Markup | | |
| Visualization / General / Radar | Radar | | |
| Visualization / General / Sunburst | Sunburst | | |
| Visualization / General / Chord | Chord | | |
| Visualization / General / Sankey | Sankey | | |
| Visualization / General / Surface Plot | Surface Plot | | |
| Visualization / General / Timelines | Timelines | | |
| Visualization / General / Word Cloud | Word Cloud | | |
| Visualization / General / Data Annotations | Data Annotations | | |
| Visualization / Geographical / Map | Map | heavy | |
| Visualization / Input and Edit / Grid | Grid | | |
| Visualization / Input and Edit / Form | Form | | |
| Visualization / Statistical / Box Plot | Box Plot | | |
| Visualization / Statistical / Correlation Plot | Correlation Plot | | |
| Visualization / Statistical / PC Plot | PC Plot | | |
| Visualization / Statistical / Pivot Table | Pivot Table | | |
| Visualization / Statistical / Statistics | Statistics | | |
| Visualization / Time and Date / Calendar | Calendar | | |

### Compute (4)

| Path | View | Tags | Ref |
|------|------|------|-----|
| Compute / Diff Studio | Diff Studio | heavy | |
| Compute / Multivariate Analysis | Multivariate Analysis | heavy | |
| Compute / PK-PD Modeling | PK-PD Modeling | heavy | |
| Compute / Bioreactor | Bioreactor | heavy | |

### Curves (2)

| Path | View | Tags | Ref |
|------|------|------|-----|
| Curves / Curve Fitting | Curve Fitting | | |
| Curves / Assay Curves | Assay Curves | | |

## 5. Results of the first run on dev (triage)

A run of all 62 demos on `dev.datagrok.ai` (2026-06-04): **59 passed, 1 known-fail, 2 fixed
in the test**.

- ✅ **59 demos** open and can be clicked without errors — including all "skip" ones (GROK-14320:
  Chemical Space, Molecule Activity Cliffs, R-Group Analysis; `demoSkip`: siRNA, Atomic Level)
  and all heavy Compute/Docking/Admetica/Retrosynthesis. This means `demoSkip` applies to the
  demo-runner, not to Browse.
- 🐞 **Bioinformatics / Similarity, Diversity** — "random clicks" on the view (the `pokeView` step)
  throw `NullError: method not found: 'gdV' on null`. This is a **regression
  [GROK-18050](https://reddata.atlassian.net/browse/GROK-18050)** (the ticket is in Done status).
  The test is **regular, with no suppression** — it fails for real when the bug triggers (it is
  intermittent, since clicking at fixed coordinates does not always land on the broken element).
  The ticket is linked via `test.info().annotations`. Details — `KNOWN_BUGS.md`.
- 🔧 **Med Chem** — the `console.error` came from a CORS block of an external Wikipedia image in
  the demo description (not a platform defect) → added to the `helpers.ts` ignore list.
- 🔧 **Database Queries** — opens an SQL-query form rather than a data view; marked `opensForm`.
  During the automated run on dev the Chembl demo DB (`db.datagrok.ai:54325`) was unavailable
  (`Connection refused`) — this is **dev infrastructure**, not a Browse/demo defect. Such a
  backend error is tolerated via `envIgnore` (any other error still fails the test). The
  navigation and form opening themselves work.

## 6. Notes and known issues

- **`skip` (`demoSkip`)** — the platform excludes these demos from the automatic demo-runner as
  unstable. In this suite they are **not skipped**: the goal is to catch real errors. Those that
  consistently fail on dev are recorded in `KNOWN_BUGS.md` and marked `test.fail()` (like the
  ModelHub cases), so that a server-side fix immediately surfaces as an "unexpected pass".
- **GROK-14320** — Chemical Space, Molecule Activity Cliffs, R-Group Analysis (Cheminformatics).
- **`heavy`** — may fail not because of the UI, but because of the docker backend being
  unavailable on dev (Admetica, Retrosynthesis, Docking, Boltz, etc.). We treat such failures
  separately from Browse/render errors.
- **Plates** — the sources contain `Plates | Assay Plates`, but it is **absent** from the Demo
  tree on dev (not deployed). Not included in the matrix.
- The list was captured automatically; when the set of demos on the server changes, regenerate the
  matrix and the `DEMOS` array.
