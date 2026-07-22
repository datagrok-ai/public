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
| **Workspace-context chars**, **tool-def count** | browser, per turn | context diet (client-visible) |

**Accuracy** is hybrid (chosen to keep authoring effort low):
- **Deterministic** — most viewer/analysis/codegen prompts carry a JS `assert` checked against the
  live workspace after the turn (e.g. *a histogram viewer exists*, *a column was added*, *rows were
  selected*). Reproducible, zero judge cost.
- **LLM-judge** — free-text answers (Help Q&A, "what's the average…") are graded by a **Haiku**
  call against a short `rubric`, returning pass/score/reason.

## Handling noise (so the numbers are credible)

- **Repetitions** — every prompt runs `reps` times (default **3**); the report shows the **median**
  (plus mean/min/max/p90/stddev in the JSON). A single LLM turn is too noisy to publish.
- **Warm-up** — one throwaway turn runs first and is discarded, priming the static-prefix prompt
  cache and the container. Measured reps then reflect steady state.
- **Caching is measured, not hidden** — `cache-read tokens` is reported; when a context-diet or
  cache-stability lever lands, you should see cache-read rise and input tokens fall.
- **Isolation** — each rep opens a **fresh** demo table (a clone) and closes it (plus anything the
  assistant opened) afterward, so viewers/columns from one rep can't leak into the next rep's assert.

## How to run (minimal involvement)

1. Log into Datagrok in the browser as usual.
2. Run the function (console, or the function browser):
   ```
   Grokky:runBenchmark("baseline")          // reps defaults to 3
   Grokky:runBenchmark("baseline", 5)        // or pick reps
   ```
   It streams progress via toast notifications and, when done, **downloads**
   `benchmark-baseline.json` + `benchmark-baseline.md`, and saves a copy to
   `System:AppData/Grokky/benchmarks/` for later comparison.
3. Change one lever (see LATENCY.md), rebuild/redeploy, then run again with a new label:
   ```
   Grokky:runBenchmark("medium-effort")
   ```
4. Diff the two:
   ```
   Grokky:compareBenchmarks("baseline", "medium-effort")
   ```
   → downloads a Markdown delta table (median turn, TTFT, input tokens, cost, accuracy; and
   per-prompt turn-time change with %).

Each lever becomes one labeled row you can show the team.

## The suite

`files/benchmark/suite.yaml` — ~20 prompts across `help`, `visualization`, `analysis`, `codegen`,
`multitool`, `query`. Fields: `category`, `prompt`, optional `table` (demo table to open),
optional `assert` (JS expression, scope `grok, DG, view, t, before, opened`), optional `rubric`
(judge). Edit freely — it's plain YAML, no rebuild needed for suite-only changes (it's a static
package file; re-publish to pick it up).

The default suite runs on `demog` + no-table prompts so it works on any server after login. To
exercise molecule/sequence or real SQL paths, add prompts with `table: SPGI` (etc.) and asserts, or
a `query` prompt against a connection you know exists.

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
