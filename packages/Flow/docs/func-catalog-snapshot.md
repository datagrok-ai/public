# Flow function-catalog snapshot (empirical)

Captured from the live `localhost` catalog with a one-off diagnostic test that
dumped every registered function's `nqName`/role/tags/input-output signature and
its computed category (since retired). The behaviour it established is locked in
by the regression tests in
[src/tests/function-browser-tests.ts](../src/tests/function-browser-tests.ts)
and [src/tests/type-map-tests.ts](../src/tests/type-map-tests.ts).

This is the data behind the function-browser classification + exclusion list.
Counts drift run-to-run as the server lazily loads more packages — treat them as
orders of magnitude, not exact.

## Scale of the problem

Latest snapshot (fully-loaded stand):

| | count |
|---|---|
| Raw `DG.Func.find({})` | **~2,357** |
| Excluded by Flow | ~2,129 |
| **Kept (shown as nodes)** | **~228** |

The previous filter kept ~568. The removed ones were almost all
**sketcher / file-handler / cell-renderer / right-click-action / core-viewer**
functions and **internal helper / near-duplicate / demo / plumbing** functions
that a scientist can't wire into a pipeline. **Widget-producing functions are
kept** (they flow to the Widgets pane and can be previewed) — only right-click
(`semantic_value`) widgets and specific dashboard-chrome widgets are dropped.
A later hand-review round (docs/things-to-do-with-functions.md) denylisted ~46
more heavy analyses / dialogs / search panels / scripts, taking ~283 → ~228.

## Exclusion rules (in [node-factory.ts](../src/rete/node-factory.ts) `shouldExcludeFunc`)

A function is dropped from the catalog when **any** holds:

0. **`meta.includeInFlow: false`** declared on the function (author opt-out, checked first; surfaces as `func.options.includeInFlow` — boolean `false` or string `'false'`). First user: `Flow:openCreationScriptFlowDialog`.
1. **Package** ∈ `EXCLUDED_PACKAGES` (`Dbtests, ApiTests, UiTests, DevTools, Tutorials, ApiSamples, UsageAnalysis`).
2. **`nqName`** ∈ the curated denylist [excluded-funcs.ts](../src/rete/excluded-funcs.ts) — individually-assessed helpers, internal twins, demo/test, and plumbing (see below).
3. **Name** is `test` or starts with `test`.
4. **Role** (comma-split, exact) ∈ `EXCLUDED_ROLES`, **or a tag** ∈ `EXCLUDED_TAGS` (case-insensitive). Checking **tags** is the single biggest lever: sketchers, cell renderers, folder viewers, `Internal`/`@editors`/`Viewers`-tagged funcs almost always declare their kind as a **tag**, not in the role field, so the old role-only check let them all through. **`EXCLUDED_TAGS` deliberately omits `panel`/`widget`/`widgets`/`tooltip`** — widget-producing functions (incl. context panels) are usable in Flow now (Widgets pane + preview).
5. Takes a **`funccall` input** — command/dialog wrapper.
6. Takes a **`semantic_value` input** — a right-click / context action that mutates UI or clipboard (this is what keeps context-widgets out while allowing real widgets).
7. Outputs a **`view`/`viewer`** (needs a TableView lifecycle; replaced by the Viewers pane) or a **`tablerowfiltercall`/`colfiltercall`** (filter-DSL builder that only feeds Aggregate/Filter internals).
8. Primitive-only (every input **and** output is a scalar) — pure formula plumbing.

### The `nqName` denylist ([excluded-funcs.ts](../src/rete/excluded-funcs.ts))

Rules 4/6/7 handle whole *categories* generically. The denylist is the
case-by-case residue that rules can't catch — produced by dumping the real
catalog and assessing each surviving function (with agents reading package
source) against *"can a scientist wire this into a data pipeline?"*. Each entry
is an internal engine/helper (`object`/handle/model/`dynamic` output), an
internal near-duplicate of a kept canonical (`Chem:getInchis` over
`Chem:addInchisTopMenu`; `Eda:apply*` over `Eda:train*`), a menu/dialog/editor
action, a demo/test/autocomplete helper, or a state-plumbing RPC. **It's meant
to be reviewed and edited by hand** — add an nqName to hide a function, remove
one to surface it.

## Classification — "what it does" (in [function-browser.ts](../src/panel/function-browser.ts) `categorizeFunc`)

Two axes:

- **Domain (by source package)** wins first, **but only for operations on data**:
  a function from a cheminformatics package (`Chem, Chembl, ChemblApi, PubchemApi,
  Chemspace, Surechembl, Admetica, Docking, Retrosynthesis, Marvin,
  ChemDrawSketcher, KetcherSketcher, HitTriage, Datagrokdsmf, Curves`) groups
  under **Cheminformatics**, a bioinformatics package (`Bio, SequenceTranslator,
  Helm, Proteomics, Bionemo, Biologics, OligoBatchCalculator, Parabilisseq,
  Sequenceutils, BiostructureViewer, PhyloTreeViewer, Peptides`) under
  **Bioinformatics** — **only if it takes a dataframe/column input**
  (`isDomainOperation`). A chem/bio function that merely *produces* a table from
  scalars (a DB fetch, generator, or query) is a **source**, not an operation, so
  it falls back to its signature category (Data Sources) — the domain sections
  stay about *doing something to your table*, never queries/sources.
  `domainSection()` / `domainCategory()` live in
  [type-map.ts](../src/types/type-map.ts) and drive the node title-bar color too
  (pink / deep-purple).
- **`DG.DataQuery` instances** are filtered out of the category tree entirely and
  shown in the **Queries pane**, grouped by connection — so all queries live
  under Queries, never under a domain section.
- **Signature category** for everything else, keyed on input/output types
  (Data Sources → Combine Tables → Transform Tables → Column Operations →
  Compute Values → Visualize → Cheminformatics → Bioinformatics → Other), in
  pipeline-build order.

Latest rendered domain sections: **Cheminformatics ~63**, **Bioinformatics ~19**
function items (queries → Queries pane; sources → Data Sources). Widgets pane
holds ~10 widget-producing functions.

## Dart-proxy note (for future catalog work)
- `getTags(func)` ([dart-proxy-utils.ts](../src/utils/dart-proxy-utils.ts)) returns a real array — `func.tags` itself is `undefined`, so it falls back to `grok_Script_Get_Tags(func.dart)`. **Don't read `func.tags` directly.**
- `func.options` is a proxy that `JSON.stringify`/`Object.entries` render as `{}`/`[object Object]`; read keys via `safeGet` — that's how `getRole` works.
- `func.nqName` is the namespace-qualified name (`Chem:getInchis`; core funcs are `core:JoinTables`) — the key used by the denylist.
