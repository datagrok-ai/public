---
name: dg-task
description: Datagrok end-to-end task orchestrator. Drives discover → plan → implement → test → critic → iterate from a single user intent, anchored on the knowledge graph at .kg/. Use this when the user describes a bug to fix, a feature to add, or behavior to change in a Datagrok plugin/library/JS API and doesn't specify the steps.
when-to-use: User says "fix bug X", "add feature Y", "change how Z works", "investigate W" with minimal detail. Or types `/dg-task <intent>`.
argument-hint: "<free-form intent — bug description, feature request, etc.>"
context: main
effort: high
---

# /dg-task — Datagrok task orchestrator

You are running a closed loop on behalf of the user. The user gave you
an intent (a bug description, a feature ask, a behavior change). Your
job is to land a working, tested, critic-approved change without
pestering them for clarification you can derive yourself.

## Hard rules

- **Read [`.claude/rules/kg-first.md`](../../rules/kg-first.md) and [`.claude/rules/critic-not-slop.md`](../../rules/critic-not-slop.md) before step 1.** They are non-negotiable.
- **Don't commit anything.** Leave a clean working tree change set for the user to review.
- **Don't edit `*.g.ts` / `*-api.ts` / `*.api.g.ts`** (a hook will block you anyway — run `grok api` to regenerate).
- **The graph is the first lookup, grep is the fallback.** If a graph query returns nothing for a question you're confident should match, fall back, but also append the discovered fact to `.kg/.learned/<date>-<slug>.md` (format in the kg-first rule).

## The loop

### 0. Bootstrap + warm the KG server (every /dg-task invocation)

If `.kg/.venv/Scripts/python.exe` (Windows) or `.kg/.venv/bin/python` (POSIX) doesn't exist, run the bootstrap:

```bash
bash .kg/scripts/bootstrap.sh   # POSIX / Git Bash
# or
.kg/scripts/bootstrap.ps1       # PowerShell
```

This creates the venv, installs requirements, and verifies the graph
loads. ~30s if the DB exists, ~3 min on a fresh clone.

Then warm the persistent query server (cheap if already running):

```bash
bash .claude/hooks/kg-warmup.sh
```

The warmup script silently no-ops if the server is already up. Subsequent
`qq.py` calls in this loop will hit the warm server in ~5–50 ms.

### 1. Discover (includes prior-art search)

Spawn the **`kg-explorer`** subagent with the user's intent verbatim.
It returns a structured digest with these mandatory sections:

```
features:    [target area — Feature ids + 1-line descriptions, ≤ 10]
prior_art:   [EXACT MATCH and/or SIMILAR (style donors) from OTHER packages]
files:       [paths most likely to need editing, ≤ 15]
tests:       [existing tests for the target features, ≤ 10]
docs:        [help pages that document the area, ≤ 5]
gaps:        [anything the agent had to grep for because the KG was empty]
notes:       [≤ 3 sentences synthesis, including reuse/mimic recommendation]
```

The `prior_art` section is **load-bearing for step 2**. If it's
missing or empty without the explorer explicitly saying "no prior
art found", reject the digest and re-spawn the explorer with the
directive "PRIOR ART section is mandatory; even 'none' must be
stated".

Do **not** read every file in the digest. The next steps will pull
just what's needed.

### 2. Plan (with prior-art + convention discovery)

This step has **two purposes**, both load-bearing:

1. **Avoid duplicating what already exists** (the reuse case).
2. **Learn how Datagrok does this kind of thing** — even when implementing
   fresh. The codebase has 79 packages and 50+ extension-point roles, and
   the "right way" to do almost any UI / grok.* / DG.* operation is
   already established somewhere. Querying the KG for analogues is how
   you discover that — it's a *learning* tool, not just a *reuse* tool.

**Always run prior-art queries**, regardless of how novel the task feels.
Even genuinely novel work benefits from cross-package convention scan:
naming, file layout, decorator shape, which `ui.*` factories are
canonical, how to structure async init, error handling, etc.

Three outcome cases for prior art:

**A. EXACT MATCH exists** — the user's intent is already implemented
elsewhere. Stop and tell the user. Format:

```
A feature that already does this exists: <feature_id> in package <X>.
Implementation: <file>
Docs: <url>

Want me to (1) just point you at it, (2) extend it for your specific
case, or (3) implement separately because <stated reason>?
```

Wait for the user's choice before proceeding. Don't silently re-implement.

**B. SIMILAR features exist (style donors)** — implement, but mimic
the donor's style. The plan must name at least one donor file the
implementer will read first:

```
Style donor: <donor_file> (from feature <donor_id> in package <X>) — same interaction_kind, similar data shape.
The new code will follow the donor's structure (file layout, decorator
shape, naming, error handling) and only diverge where the domain
genuinely differs.
```

**C. Truly novel** — no comparable feature exists. **Still run the
convention-discovery queries below.** Pick the most idiomatic example
of each cross-cutting concern (UI building, dapi access, event
subscription, semType handling) and cite it as a pattern reference.
The implementer will read those alongside the new code spec.

#### 2.1 Convention discovery — what to query

If the task touches any of these surfaces, query the KG for *how it's
already done* in 3-5 comparable places. This is **separate from** the
exact-match question — it's about adopting the codebase's idioms.

| Task touches… | Query |
|---|---|
| **UI building** (dialogs, inputs, panels, layouts, drag-drop) | `MATCH (m:TsMethod)-[u:USES_UI_COMPONENT]->(c:UiComponent {name:'<factory>'}) RETURN m.id, u.use_count ORDER BY u.use_count DESC LIMIT 10` — find the methods that use this `ui.*` most heavily and read 2-3 of them for the canonical idiom (where to dispose, how to wire async, where to put the OK callback). Then the donor file list goes into the plan. |
| **`grok.dapi.*` calls** (server entities, queries, files, docker) | `MATCH (m:TsMethod)-[c:CALLS_DAPI_ENDPOINT]->(:DapiEndpoint {name:'<endpoint>'}) RETURN m.id, c.call_count ORDER BY c.call_count DESC LIMIT 10` — and if your task is about a specific entity served by a dapi endpoint, also `MATCH (e:DapiEndpoint)-[:TYPED_FOR]->(c:TsClass) WHERE c.name='<EntityClass>' RETURN e.name, e.accessor_path` |
| **`grok.events.*` subscriptions** | `MATCH (m:TsMethod)-[s:SUBSCRIBES_TO_EVENT]->(:JsEventStream {name:'<event>'}) RETURN m.id, s.subscribe_count ORDER BY s.subscribe_count DESC LIMIT 10` — read 2-3 to see how the codebase handles cleanup, scoping, and current-view filtering. |
| **`DG.<class>` API** (DataFrame, Grid, Viewer, GridCellRenderer, etc.) | `MATCH (m:TsMethod)-[u:USES_API_CLASS]->(:TsClass {name:'<X>', is_jsapi:true}) RETURN m.id, u.use_count ORDER BY u.use_count DESC LIMIT 10` — and `MATCH (sub:TsClass)-[:EXTENDS_CLASS]->(:TsClass {name:'<X>', is_jsapi:true}) RETURN sub.id LIMIT 10` if you're subclassing it. |
| **A specific semType** (Molecule, Macromolecule, OligoNucleotide, etc.) | `MATCH (rf:RegisteredFunction)-[:CONSUMES_SEMTYPE]->(:SemanticType {name:'<X>'}) RETURN rf.id, rf.role LIMIT 15` — shows every function that already operates on this semType, so you can see how others handle validation, defaulting, error UI. |
| **Cell renderers / cell editors / panels / file viewers** (any platform-dispatched function) | `MATCH (rf:RegisteredFunction)-[:HAS_TAG]->(:Tag {name:'<role>'}) RETURN rf.id, rf.package_id, rf.meta` — every example with column/file scope. |
| **Tutorials / demos** (the user-visible counterpart of the new feature) | `MATCH (t:Tutorial)-[:DEMONSTRATES]->()<-[:PART_OF_FEATURE]-(:Feature {id:'<related_feature>'}) RETURN t.name, t.help_url` |

The returned method/file IDs are *donor candidates*. Pick the 2-3 with
the highest use_count or that come from packages closest in domain
(Bio for sequences, Chem for molecules, Charts for viewers, PowerGrid
for grid customizations, etc.) and add them to the donor list in the
plan, even when no single feature is a "style donor" for the whole
intent.

#### 2.2 Draft the plan

After the prior-art + convention discovery, draft a 3-6 step plan:

```
1. Reproduce / locate: <which file:line carries the suspected defect, or which file gets the new code>
2. Convention reads (mandatory): <2-3 donor file(s) the implementer must read before editing — pick from the convention queries above>
3. Donor verify (only if delegating across packages): <see step 2.5>
4. Change: <one-sentence description of the edit, referencing the convention donors by name where applicable>
5. Test: <which test file gets a new assertion, naming the assertion's distinguishing inputs per critic-not-slop>
6. Verify: grok test scope (--category or --test prefix)
7. (optional) Doc update: <which help-page section, only if user-facing behavior changed>
```

If the intent is ambiguous in a way you can't resolve from the graph
or 1-2 quick file reads, **ask the user one question** before
proceeding. Otherwise, plow on.

### 2.5. Donor verification (mandatory when case A or B)

If the plan says "delegate to", "reuse", or "call into" an existing
function from another package/library, you MUST read the donor's
**body** (not just its registration metadata or signature) before
the implementer runs. The KG and the explorer's digest only see the
declared signature — runtime preconditions live in the body.

Read the donor file and list every precondition it imposes on its
caller's data:

| Precondition kind | What to look for |
|---|---|
| Column semType / units / tags | `seqHelper.getSeqHandler(col)`, `col.semType !== ...`, `col.meta.units !== ...`, throws / asserts on tag mismatch |
| Required column shape | `col.length === 0` guards, type checks, regex validation |
| Globals it depends on | `_package.helmHelper` (must be initialized), `grok.shell.tv` (requires an open table view), monomer lib loaded |
| Side effects on caller's state | mutates `col.meta`, dirties the dataframe, opens a global dialog |

Then compare against your caller's data shape. Write one of:

```
DONOR VERIFIED: <donor_id>'s preconditions (<list>) are satisfied by our caller's data (<our shape>). Delegation safe.
```

```
DONOR INCOMPATIBLE: <donor_id> requires <precondition>; our caller has <our shape>. Cannot delegate.
Plan downgraded from "reuse" to "implement directly using <shared helper>".
```

A `DONOR INCOMPATIBLE` result requires the plan to be rewritten before
step 3. The shared helper named in the new plan must be one the donor
itself uses internally (e.g. `helmHelper.createWebEditorApp` is what
Helm's `editMoleculeCell` uses under the hood — that's the right level
to reuse from).

This step costs 30–60 seconds. Skipping it cost a full extra iteration
in the OligoNucleotide cellEditor task (May 2026): the donor
`Helm:editMoleculeCell` looked reusable from its signature but its body
called `seqHelper.getSeqHandler(col)` which throws on non-Macromolecule
columns. Test passed on metadata, runtime threw on first double-click.

### 3. Implement

Spawn the **`dg-implementer`** subagent with:
- the plan (verbatim)
- the digest (file paths only, not contents — the implementer reads
  what it needs)
- a directive: "make minimal edits; do not refactor adjacent code; do
  not add error handling for cases that can't happen; preserve existing
  styles per `.claude/rules/code-style.md`"

The implementer returns a unified diff summary (file paths + LOC
changed) plus a list of any new files created.

### 4. Test

Spawn the **`dg-tester`** subagent with the diff summary, the planned
test description, and the package(s) that changed. It:
- writes the test in the right place (`src/tests/<name>.ts`,
  registered in `package-test.ts`)
- runs `grok test --host localhost --skip-build --no-retry --test "<prefix>"` first; if TS changed, drops `--skip-build`
- returns pass/fail + any failing assertions verbatim

**Interaction-test directive.** When the user's intent describes a UI
interaction (verbs like *click*, *double-click*, *drag*, *hover*,
*open*, *select*, *paste*, *resize*) OR registers a platform-dispatched
function (`cellEditor`, `cellRenderer`, `panel`, `fileViewer`,
`fileExporter`, `semTypeDetector`), the orchestrator's directive to the
tester MUST include:

> The test MUST construct the dispatching context (DataFrame +
> TableView + GridCell with the matching column tags / file with
> matching extension / semantic value with matching semType) and
> invoke the entry point through the platform's dispatch path
> (`Func.find(...).apply({...})` for cellEditor / panel / etc.;
> simulated DOM event for click / drag / hover). Then assert the
> user-visible outcome (dialog opens, cell value updates, panel
> renders the expected DOM). A `DG.Func.find` registration check
> alone is NOT acceptable — it is necessary smoke, never sufficient
> proof. See `.claude/rules/critic-not-slop.md` § "Registration-only
> tests".

If the test fails for a reason unrelated to the fix (env, missing
fixture), the tester investigates and reports back; you decide
whether to retry or surface to the user.

### 5. Critique

Spawn the **`dg-critic`** subagent (read-only) with:
- the diff
- the test file(s) it added/changed
- the plan it was meant to fulfil

It returns findings grouped by severity. **Blocking severity** must
include at least:
- a slop-test detection (per `critic-not-slop.md`)
- a generated-file edit
- a violation of an explicit rule in `.claude/rules/`
- a missing test for behavior that's now reachable

### 6. Iterate

If the critic raises blocking findings, loop back to step 3 with the
findings appended to the implementer's directive. **Max 3 iterations.**
If still blocked at iteration 3, write a brief postmortem to the user
explaining what was tried and what's stuck.

### 7. Refresh KG (explicit — runs as part of /dg, not on session Stop)

When all green:

1. Append every touched package folder name to `.kg/.dirty-packages`
   (one per line, deduped). The orchestrator (you) does this in the
   main loop after the implementer reports back.
2. Run the refresh script:
   ```bash
   bash .claude/hooks/refresh-kg.sh
   ```

The script (a) drains `.kg/.learned/*.md` into the canonical
edge JSONL via [`.kg/scripts/learned_to_enrichment.py`](../../../.kg/scripts/learned_to_enrichment.py),
moving applied files into `.learned/applied/`, and (b) runs
`build.py --packages <dirty>` to re-extract just the touched packages
and rebuild Kuzu. Cheap (~10–30s per package) because we're not
re-running enrichers — only the deterministic extractors.

The script is **NOT wired as a session-Stop hook** (intentionally
removed 2026-05-10). It only runs when the `/dg-task` orchestrator calls
it. If a /dg-task run is interrupted before step 7, the queues persist —
the next /dg-task will pick them up.

### 8. Report back

Single message to the user. Format:

```
## What changed
- <file>: <one-line summary>

## How it was tested
- <test file>: <assertion that distinguishes>
- grok test result: <PASS / FAIL summary>

## Critic verdict
- <PASS / list of non-blocking notes>

## Knowledge graph
- Marked dirty: <package list>  (Stop hook will refresh)
- Learned: <count> new fact(s) appended to .kg/.learned/

## Suggested next
- <if any obvious follow-up>
```

Keep it under 30 lines. The user will read the diff themselves if
they want detail.

## Shortcuts you should always take

The KG cookbook lives at [`.kg/docs/QUERYING.md`](../../.kg/docs/QUERYING.md). The high-leverage patterns for `/dg-task`:

### Discover step (kg-explorer)

| Intent | Query |
|---|---|
| "Fix bug in X" — find candidate features | `MATCH (f:Feature) WHERE toLower(f.name + ' ' + coalesce(f.description, '')) CONTAINS 'x' RETURN f.id, f.name, f.description LIMIT 5` |
| "What's already in package P?" | `MATCH (p:Package {name:'P'})-[:HAS_FEATURE]->(f) RETURN f.name, f.interaction_kind` |
| "Where is X tested?" — at the test-block level | `MATCH (t:PackageTest) WHERE toLower(t.name) CONTAINS 'x' OR toLower(t.category) CONTAINS 'x' RETURN t.package_id, t.category, t.name, t.file_path` |
| "Where is class C defined?" | `MATCH (c:TsClass {name:'C'})-[d:DEFINED_IN]->(f:File) RETURN f.relative_path, d.line_start` |

### Donor lookup (Plan step 2 — find prior art)

| Intent | Query |
|---|---|
| "All cellEditors / cellRenderers / panels / fileViewers across packages" | `MATCH (rf:RegisteredFunction)-[:HAS_TAG]->(:Tag {name:'<role>'}) RETURN rf.id, rf.package_id` |
| "Which method opens a `ui.dialog`? — donor candidates" | `MATCH (m:TsMethod)-[u:USES_UI_COMPONENT]->(:UiComponent {name:'dialog'}) RETURN m.id, u.use_count ORDER BY u.use_count DESC LIMIT 20` |
| "Pattern: cellEditor with column scope" | `MATCH (rf:RegisteredFunction)-[:HAS_TAG]->(:Tag {name:'cellEditor'}) OPTIONAL MATCH (rf)-[:REQUIRES_COLUMN_TAG]->(s:SemanticType) RETURN rf.id, collect(s.name)` |

### Donor verification (step 2.5 — does the donor's preconditions match my caller?)

| Intent | Query |
|---|---|
| "What column shape does Helm:editMoleculeCell require?" | `MATCH (rf:RegisteredFunction {id:'func:Helm:editMoleculeCell'})-[r:REQUIRES_COLUMN_TAG]->(s:SemanticType) RETURN r.tag_key, s.name` |
| "What helpers does the donor's body call (for direct reuse)?" | `MATCH (rf:RegisteredFunction {id:'func:Helm:editMoleculeCell'})-[:DEFINED_BY_FUNCTION]->(m:TsMethod)-[u:USES_API_CLASS]->(c:TsClass) RETURN c.name, u.use_count` |

If `REQUIRES_COLUMN_TAG` returns a value your caller doesn't satisfy
(e.g. `quality=Macromolecule` but you have `quality=OligoNucleotide`),
**downgrade the plan from "delegate" to "implement directly using
shared helper"** per the donor-verification rule.

### Test step (where tests live, what they cover)

| Intent | Query |
|---|---|
| "All test() blocks in package P matching keyword" | `MATCH (t:PackageTest) WHERE t.package_id='pkg:P' AND toLower(t.name) CONTAINS 'k' RETURN t.category, t.name, t.file_path` |
| "Tests covering function F" | `MATCH (t:PackageTest)-[:COVERS]->(:RegisteredFunction {name:'F'}) RETURN t.file_path, t.name` |
| "Test file siblings — what's in the same test file?" | `MATCH (other:PackageTest {file_path: '<path>'}) RETURN other.category, other.name` |

### Coupling / blast radius (impact analysis)

| Intent | Query |
|---|---|
| "Who imports this file?" | `MATCH (other:File)-[i:IMPORTS]->(:File {id:'file:<path>'}) RETURN other.relative_path, i.imported_symbols` |
| "Who imports symbol S from anywhere?" | `MATCH (f:File)-[i:IMPORTS]->(t) WHERE 'S' IN i.imported_symbols RETURN f.relative_path, t.id` |
| "Cross-package callers of function F" | `MATCH (caller:Package)-[c:CALLS]->(:RegisteredFunction {id:'func:Pkg:F'}) RETURN caller.name, c.call_count` |
| "Subclasses of DG.X (impact of changing X)" | `MATCH (sub:TsClass)-[:EXTENDS_CLASS]->(:TsClass {name:'X', is_jsapi:true}) RETURN sub.id` |

## When to abort and ask

- The user's intent maps to >1 plausible Feature with very different
  semantics (e.g. "fix the renderer bug" — there are 30 renderers).
- The fix would require touching the JS API surface (`js-api/`)
  or core Dart (`core/`). Those are out of scope for this skill.
- The change requires a new external dependency (npm or pip). Always
  ask first.
