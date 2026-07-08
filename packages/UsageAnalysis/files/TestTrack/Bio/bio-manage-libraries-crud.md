---
feature: bio
target_layer: playwright
coverage_type: regression
priority: p0
realizes: [bio.cp.manage-monomer-libraries]
produced_from: atlas-driven
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
realized_as:
  - bio-manage-libraries-crud-spec.ts
gate_verdicts:
  f:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T04:15:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-02-bio-automate-01
    timestamp: 2026-06-02T12:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-02-bio-automate-01
    timestamp: 2026-06-02T14:05:00Z
    spec_runs:
      - spec: bio-manage-libraries-crud-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 100
        failure_keys: []
  a:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T11:30:00Z
    review_round: 1
    failure_keys: []
---

# Bio — Manage Monomer Libraries: app, tree browser, Monomers view & Match dialog

Covers three sibling **Bio | Manage** surfaces that the shorter
`manage.md` smoke test and the lifecycle scenarios
(`bio-lifecycle-monomer-library.md` /
`bio-lifecycle-monomer-collection.md`) only touch indirectly:

- the standalone **Manage Monomer Libraries** app and its
  tree-browser sidebar (listing each library, click to open its
  per-library manager),
- the **Bio | Manage | Monomers** view (full CRUD for individual
  monomers, with SMILES↔molfile standardization),
- the **Bio | Manage | Match with Monomer Library...** dialog
  (matches molecules in a column against monomers in a chosen
  library, for PEPTIDE / RNA / CHEM polymer types), and
- the underlying `standardiseMonomerLibrary` normalization API
  that library edits run through.

## Setup

Dataset for the Match-with-library scenario: a small-molecules
table loaded from
`System.AppData/Bio/tests/filter_HELM.csv` — its sequence column
provides the HELM-shape input the Match dialog's `Molecules`
column input picks up. (The dataset is required only by
Scenario 3; Scenarios 1 and 2 drive the standalone app + the
Monomers view without a dataset.)

Server-state dependencies:

- The recon account on dev.datagrok.ai surfaces four monomer
  libraries (`HELMCoreLibrary.json`, `NH2.json`,
  `polytool-lib.json`, `sample-lib-Aca-colored.json`) per
  `.claude/skills/grok-browser/references/bio.md#L479-L480`. The
  list varies per-user (FileShare permissions on
  `System:AppData/Bio/monomer-libraries`); assert STRUCTURE
  (≥1 library row, per-row controls present) rather than a fixed
  library count.
- The `Manage Monomer Libraries` app (`bio.manage.libraries-app`)
  is reachable both via the `Peptides` browse path and via direct
  invocation (`grok.functions.call('Bio:libsApp')` per
  `package.ts#L1377` `name: libsApp`); the scenarios below take
  the direct-invocation path to keep the entry deterministic
  across user accounts whose browse-tree state may differ.
- The `Manage Monomers` view (`bio.manage.monomers-view`) is
  driven entirely via its top-menu entry and its UI mirrors the
  `Manage Monomer Libraries` view shape per
  `bio.md#L484-L488` (the recon log records this as out-of-scope
  for selector recon on Bio 2.26.5; Scenario 2 below asserts the
  view-open contract only — see Notes on the bounded-assertion
  posture).

## Scenarios

### Scenario 1 — `Manage Monomer Libraries` app + tree browser

Steps:

1. Invoke the standalone `Manage Monomer Libraries` app entry
   point: `await grok.functions.call('Bio:libsApp')`.
2. Confirm an app-shell view opens with `grok.shell.v.name`
   resolving to the `Manage Monomer Libraries` title (or the
   equivalent app-root view).
3. Drive the `Manage Monomer Libraries` tree browser
   programmatically: `await grok.functions.call(
   'Bio:appTreeBrowser', {treeNode: <root-node-of-Peptides>,
   browseView: <browseView>})`.
4. Verify the tree-browser side panel renders one node per
   library available to the current user (≥1 library node;
   exact set varies per `FileShare` per the recon note above).
   Each library node MUST be clickable and MUST open its
   per-library monomer manager when clicked (the per-library
   monomer manager surface is verified at the structural
   "panel opens" level — assertion bar matches the bounded
   posture documented in Notes).

Expected:

- The `Bio:libsApp` invocation resolves (no thrown error, no
  rejected promise).
- The app view opens with the `Manage Monomer Libraries` title.
- The `Bio:appTreeBrowser` invocation populates the tree
  sidebar with ≥1 library node (number varies with user
  FileShare permissions; the scenario asserts presence and
  click-ability, not a fixed count).
- Clicking a library node opens its per-library manager (a
  child panel inside the app shell; the manager renders without
  raising an error balloon).
- No error balloon raised during the sequence
  (`isErrorBallon` returns `false`).

### Scenario 2 — `Bio | Manage | Monomers` view open

Steps:

1. On the menu ribbon, go to **Bio** > **Manage** > **Monomers**.
   The selector is
   `[name="div-Bio---Manage---Monomers"]` per
   `.claude/skills/grok-browser/references/bio.md#L488` — the
   leaf opens a VIEW (not a `.d4-dialog`), mirroring the shape
   documented for `Bio | Manage | Monomer Libraries` (the
   `manageMonomersView` function, browseable from the top menu).
2. Wait for the view to mount: `grok.shell.v.name` resolves to
   the `Manage Monomers` title (or the runtime-specific
   equivalent; bio.md L488 flags this view as
   out-of-scope-for-selector-recon on Bio 2.26.5, so the
   assertion bar is the view-open contract rather than per-row
   controls).
3. Verify the view exposes at minimum:
   - a root container (`grok.shell.v.root` is a non-null
     `HTMLElement`),
   - the `Manage Monomers` heading or equivalent title-bar
     text inside the view root,
   - ≥1 child element representing the monomer-list surface
     (a list / grid / table container — exact selector left
     bounded per the recon-log note).

Expected:

- The top-menu click opens a view via `grok.shell.v` (NOT a
  `.d4-dialog`).
- `grok.shell.v.name` matches the `Manage Monomers` title.
- The view root is a non-null element with at least one
  child element representing the monomer-list surface.
- No error balloon raised.

### Scenario 3 — `Bio | Manage | Match with Monomer Library...` dispatch + `standardiseMonomerLibrary` normalization

Steps:

1. Open `System.AppData/Bio/tests/filter_HELM.csv` so a HELM
   sequence column is available as the `Molecules` input for the
   Match dialog. The Macromolecule detector classifies the
   sequence column.
2. Dispatch the top menu: click
   `[name="div-Bio---Manage---Match-with-Monomer-Library..."]`
   per `.claude/skills/grok-browser/references/bio.md#L490-L498`.
3. Wait for the dialog to mount:
   `[name="dialog-matchWithMonomerLibrary"]` (title
   `matchWithMonomerLibrary` — note the raw function-name title
   is a documented minor inconsistency per `bio.md#L494`).
4. Verify the dialog exposes three host inputs:
   - `[name="input-host-Table"]` — current table picker.
   - `[name="input-host-Molecules"]` — sequence/molecules
     column picker; defaults to the detected Macromolecule
     column.
   - `[name="input-host-Polymer-Type"]` — polymer type select
     (`PEPTIDE` / `RNA` / `CHEM` per atlas).
5. Confirm the Polymer Type select carries the three atlas-
   declared options (`PEPTIDE`, `RNA`, `CHEM`).
6. Library normalization side: programmatically invoke the
   `standardiseMonomerLibrary` API exposed as
   `grok.functions.call('Bio:standardiseMonomerLibrary',
   {library: <library-JSON-object>})`. Pass a minimal HELM
   JSON-shape monomer-library object (e.g. the normalized form
   of a single PEPTIDE monomer with `symbol`, `name`, `molfile`,
   `polymerType: PEPTIDE`); capture the returned normalized
   object.

Expected:

- The Match dialog opens with the three host inputs present
  (`Table`, `Molecules`, `Polymer-Type`).
- The `Polymer-Type` select carries the three atlas options.
- The `standardiseMonomerLibrary` invocation resolves without
  error and returns a non-null normalized library object.
- The normalized library is structurally compatible with the
  HELM JSON library schema (object with a `monomers` /
  `polymerSchemas` shape per the standardization pipeline; the
  scenario asserts non-null-and-object rather than a per-field
  match, because the standardization output for a minimal
  single-monomer input is canonicalization rather than
  enrichment).
- No error balloon raised during dialog open or the
  standardization invocation.

## Notes

- Complements `manage.md` (a quick smoke check of the
  **Bio | Manage | Monomer Libraries** view) by covering the three
  sibling **Manage** entries it doesn't exercise (**Monomers**,
  **Match with Monomer Library...**) plus the standalone app and
  tree-browser surface.
- Deferrals: the per-library monomer-manager's content (Scenario 1,
  step 4) and the **Manage Monomers** view's per-row controls
  (Scenario 2) are only checked at the "panel opens, no error"
  level, not down to individual row selectors — the underlying view
  shapes weren't available for detailed selector recon at the time
  this scenario was written.

---
{
  "order": 24
}
