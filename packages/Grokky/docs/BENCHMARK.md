# Grokky benchmark harness — measuring latency & accuracy

A repeatable way to put **numbers** on every latency lever in [LATENCY.md](./LATENCY.md) before and
after you change it, so improvements can be proven to the team instead of asserted.

Companion to LATENCY.md: that doc says *what* to speed up and "instrument first". This is the
instrument.

## What it measures

Each benchmark run drives a fixed suite of prompts through the **real runtime pipeline** — the same
path the panel uses (`buildWorkspaceContext()` prepend, `viewFunctionTools()` declared, real
`datagrok-exec`/`datagrok_verify` tool round-trips) — so the numbers reflect true user-perceived
turns, not a synthetic shortcut.

Per turn it records:

| Metric | Source | Lever it exposes |
|---|---|---|
| **TTFT** (ms) | browser: send → first streamed token | effort / thinking / model / cache |
| **Total** (ms) | browser: send → `final` | everything |
| **Tool round-trips** | count of `tool_activity` events | forced verify / grounding round trips |
| **Input / output / cache-read tokens** | SDK `result.usage`, forwarded on `final` | context diet, prompt-cache stability |
| **Cost (USD)** | SDK `result.total_cost_usd` | context diet, model choice |
| **num_turns**, api duration | SDK `result` | round-trip count, server-side time |
| **Tool names**, in call order | `tool_activity.name` | *which path* a turn took — the sharpest cheap-vs-strong discriminator |
| **Unverified rate** | `final.unverified` | how often the model could not confirm its own action |
| **Workspace-context chars**, **tool-def count** | browser, per turn | context diet (client-visible) |

**Accuracy** is hybrid (chosen to keep authoring effort low):
- **Deterministic** — most viewer/analysis/codegen prompts carry a JS `assert` checked against the
  live workspace after the turn. Reproducible, zero judge cost.
- **LLM-judge** — free-text answers (Help Q&A, "what's the average…") are graded by a **Haiku**
  call against a short `rubric`, returning pass/score/reason. Every rep that produced text is
  judged and the verdict is the majority, with a score floor — judging one rep while timing three
  let the same defect flip pass/fail between runs.

### Assertion policy — assert what was asked for

An **existence** check (*a scatter plot exists*) passes on a half-correct answer. When every prompt
passes at rate 1.0, the accuracy axis carries no information and any comparison between
configurations is meaningless. So asserts **bind the specifics**: which columns a viewer is on,
the exact number of rows selected, that a formula column really holds the requested formula.

- `v.props.<name>` reads a viewer's **live** look including defaults, so it is safe to read even
  when the property was never set explicitly.
- Prefer **self-referential** ground truth — recompute the expected value from `t` rather than
  hardcoding a count — so the suite survives a demo-file refresh.
- **Missing numerics are not `null`.** For a numeric column, `col.toList()[i]` gives `undefined`
  and `col.get(i)` gives a `2.6789344063684636e-34` sentinel. So `v !== null` counts missing rows
  as present and `v === null` never matches one — this silently shifted two counts by exactly the
  number of missing rows and **failed correct model answers**. Use `col.isNone(i)`. A threshold
  comparison (`v > 50`) already excludes `undefined`, so a `v !== null &&` guard next to one is
  dead code that merely looks like it handles missing values. `lint-suite.mjs` flags both.
- **Verify hardcoded constants before trusting them.** A rubric demanding the answer `1274` where
  the true count was `878` marked a correct answer wrong for three runs. Check with
  `dev/harness/introspect.mjs` and prefer recomputing from `t`.
- A turn that fails every rep now scores as **failed**, not "unscored". Leaving it null dropped it
  from the denominator, which is how a 19-prompt run once reported "16/18".
- A failing prompt reports **what the model actually left behind** (viewer types with their bound
  columns, added/removed columns, row/filter/selection counts, grid-hidden columns), so deciding
  "wrong assert or wrong model?" is reading the report rather than re-running the prompt.
- An assert that *throws* is reported separately as a harness bug, not silently counted as a model
  failure.

## Handling noise (so the numbers are credible)

- **Repetitions** — every prompt runs `reps` times (default **3**); the report shows the **median**
  (plus mean/min/max/p90/stddev in the JSON). A single LLM turn is too noisy to publish.
- **Warm-up** — one throwaway turn runs first and is discarded, priming the static-prefix prompt
  cache and the container. Measured reps then reflect steady state.
- **Caching is measured, not hidden** — `cache-read tokens` is reported; when a context-diet or
  cache-stability lever lands, you should see cache-read rise and input tokens fall.
- **Isolation** — each rep opens a **fresh** demo table (a clone) and closes it (plus anything the
  assistant opened) afterward, so viewers/columns from one rep can't leak into the next rep's assert.

## How to run

### Headless (preferred) — `dev/harness/run-benchmark.mjs`

Drives a real browser with Playwright, logging in with the dev key from `~/.grok/config.yaml`
exactly the way `grok test` does, and writes the report into `files/benchmarks/` so results land in
the repo rather than a Downloads folder. A benchmark **cannot** use the browserless T2 driver: its
asserts read live `DG` objects (viewers, filters, grid state) that only exist in a tab.

```bash
cd dev/harness
node run-benchmark.mjs --label sonnet-baseline --reps 3 --model sonnet
node run-benchmark.mjs --label smoke --reps 1 --only trivial,visualization
node run-benchmark.mjs --compare sonnet-clean,haiku-probe
```

`--host` picks a server alias (default `localhost`), `--headed` shows the browser.

### From the browser

```
Grokky:runBenchmark("baseline")                          // reps defaults to 3
Grokky:runBenchmark("pin-haiku", 3, "haiku")             // pin every turn to one model
Grokky:runBenchmark("smoke", 1, null, "trivial,help")    // run part of the suite
Grokky:compareBenchmarks("baseline,pin-haiku")           // 2+ comma-separated labels
```

- **`model`** pins every turn to one tier, producing the **control arms** for a model
  comparison without redeploying the runtime.
- **`only`** narrows the suite to comma-separated categories, difficulties, or prompt
  substrings. A full arm is 44 × reps turns at ~20 s each, so validating a suite change or probing
  one tier should never cost a full run. Partial runs are stamped as such in the report and must
  not be compared against full arms.
- **`compareBenchmarks`** is N-way, not pairwise, so more than two configurations can be read
  together. The first label is the reference. It reports accuracy **per difficulty tier** and per
  prompt shows turn time *with* the pass mark — a faster arm that quietly stopped passing a prompt
  is the failure mode a latency-only table is blind to.

Reports are written to `System:AppData/Grokky/benchmarks/` (JSON + Markdown) and downloaded.

## First result: model tier is not the latency lever

The first thing this instrument was built to answer was whether a cheap-model-first dispatcher
(route to Haiku, escalate to Sonnet when it struggles) would cut latency. Two full arms
(`benchmark-sonnet-clean`, `benchmark-haiku-probe`, in `files/benchmarks/`) said no:

| | Sonnet | Haiku |
|---|---|---|
| Accuracy | **39/44 (89%)** | 27/44 (61%) |
| Median turn | **26.9 s** | 39.6 s — *47% slower* |
| Cost | $5.42 | $2.09 (only **1.8×** cheaper per *passed* prompt) |
| trivial tier | 10.8 s / **1** tool round-trip | 43.1 s / **8** round-trips |
| visualization | 13.9 s / 2 round-trips | 83.8 s / **20** round-trips |

The cheap tier is the *slow* one here. Every tool round-trip is a browser round-trip, and Haiku
explores the API where Sonnet knows it — so it is 4-6× slower on exactly the easy prompts a
dispatcher would route to it. Escalation would make it worse: pay Haiku's flailing, then Sonnet's
work. The router was built, measured, and removed.

Where the time actually goes:

```
sonnet: TTFT 17.9s of 26.6s total (67%), median 2 round-trips
haiku:  TTFT  7.7s of 39.0s total (20%), median 6 round-trips
```

**Two-thirds of a Sonnet turn is spent before the first token.** That points the latency work at
TTFT (context diet, prompt-cache stability, thinking budget) and at eliminating round-trips — not
at model choice. Use `--model` to re-check that conclusion when models change; it is a property of
this architecture, not a permanent law.

## The suite

`files/benchmark/suite.yaml` — 44 prompts across `help`, `visualization`, `analysis`, `codegen`,
`multitool`, `query`. Fields: `category`, `prompt`, optional `difficulty`, `table`, `assert`,
`rubric`. Edit freely — it's plain YAML, no rebuild needed for suite-only changes (it's a static
package file; re-publish to pick it up).

**Assert scope:** `grok, DG, view, t, before, opened, openedViews, tools`, where `before` is
`{cols, rowCount, viewers, tableViews}` (note `viewers` counts the grid, so a fresh view is 1),
`openedViews` are the TableViews the turn opened, and `tools` are the bare tool names it called in
order. `await` is allowed. Asserts run after a short settle delay, because a Filters-panel filter
applies asynchronously — without it every `t.filter` assertion is a race.

**`table:`** resolves to `System:DemoFiles/<table>.csv`, so nested paths are written `chem/SPGI`.
`randomWalk` is special-cased to `grok.data.demo.randomWalk()` (its columns are literally `#0`,
`#1`, `#2` — a useful name-quoting stress test).

**`difficulty:`** (`trivial` | `standard` | `hard`, default `standard`) drives a per-difficulty
rollup in the report. The tiers exist to make model-routing questions answerable: a cheap tier
reveals where configurations differ; a blended number hides it.
Hard prompts deliberately target documented traps — clearing a filter vs a selection uses opposite
polarity; `grid.sortByColumns =` is a read-only getter whose assignment silently no-ops; a
substructure filter given raw SMILES instead of a molblock matches zero rows — so a model that
skips the skill produces a *detectably wrong* state rather than a plausible one.

**Lint before you run.** A run costs ~20 minutes and real quota, so validate first:
```bash
cd dev/harness && node lint-suite.mjs
```
It parses the YAML, compiles every assert exactly the way the harness does, and flags assertions
too weak to carry the accuracy axis. It is also part of `node dev/check.mjs`.

## ⚠️ Deploy note — token/cost metrics need a runtime image rebuild

Latency, tool-round-trip, and context-size numbers work **immediately** (they're measured in the
browser). The **token/cost** columns come from the SDK `result` message, which the runtime now
forwards on the `final` event ([session.ts](../dockerfiles/claude-runtime/src/session.ts),
[types.ts](../dockerfiles/claude-runtime/src/types.ts)). That code compiles **inside** the
`claude-runtime` image — `npm run build` / `grok publish` rebuild only the browser bundle. Until you
rebuild and redeploy the image (`docker build` in `dockerfiles/claude-runtime/`), the report shows a
"token metrics unavailable" note and the token/cost columns are omitted. Every LATENCY.md lever
lives in that same image, so you rebuild it anyway — the metrics forward rides along.

## Caveats / scope

- The runner sends prompts headlessly; **AskUserQuestion** can't be answered, so any prompt that
  triggers it is unblocked with an empty answer and won't score well — keep suite prompts specific.
- Server-side config (effort, thinking budget, model default) is **not** visible to the browser, so
  it isn't auto-captured — the run **label** is how you denote which config produced it. Client-visible
  config (tool-def count, workspace-context size, model the browser sends) *is* recorded.
- `total` includes network + runtime queue + model; it is wall-clock as the user feels it, not pure
  model time (`duration_api_ms` in the JSON isolates the server-side portion).
