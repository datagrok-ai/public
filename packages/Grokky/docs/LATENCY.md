# Grokky latency — measured findings and plan

Rewritten 2026-07-20 from **measured data** (benchmark harness, [BENCHMARK.md](./BENCHMARK.md)) plus a
full code sweep of the turn path (browser → WS → runtime → SDK → model → tool round-trips). Replaces
the earlier speculative version — several of its guesses are now disproven (see "Disproven" at the end).
Refs are to `dockerfiles/claude-runtime/src/` unless noted. Raw runs: `files/benchmarks/`.

## Measured baseline (Sonnet, effort=high, thinking=1500)

19 prompts × 3 reps through the real pipeline (`Grokky:runBenchmark`):

| | baseline (high) | medium-effort A/B |
|---|---|---|
| Median turn | 19.4 s | 15.9 s |
| Median TTFT | 17.7 s | 12.8 s |
| Help/doc turns | 42–61 s, 8–13 round-trips | 44–85 s, **20–33 round-trips** |
| Simple action turns (viz) | 11–15 s, 2 round-trips | 11–20 s, 2 round-trips |
| Accuracy | 16/18 | 18/19 |
| Cost (Σ per-prompt medians) | $1.62 | **$2.16 (+33%)** |
| Failures | 3 (all 120 s timeouts, doc-search spiral) | 4 (same mechanism) |

Linear fit across prompts:

```
total_ms ≈ 12,400 + 3,600 × tool_round_trips   (high,   R²=0.47)
total_ms ≈ 14,200 + 2,300 × tool_round_trips   (medium, R²=0.72)
```

**The two numbers that matter: ~12–14 s fixed overhead per turn, ~2.3–3.6 s per tool round-trip.**
Round-trip *count* is the multiplier — each hop re-reads the ~35 k-token cached prefix and pays a full
model continuation. On action turns TTFT ≈ total because the model emits text *after* its tool calls
(output ordering, not a transport defect) — so the whole turn is user-visible spinner time.

### The effort A/B result (why "just lower effort" failed)

Lowering effort cut the per-hop cost (3.6 → 2.3 s) but the model compensated for less deliberation
with **more searching** — help turns went from 8–13 to 20–33 round-trips, net slower and +33 % cost.
**Effort stays `high`**: deliberation is what keeps the round-trip count down while the grounding
gate stands. The two levers are coupled — effort can only come down together with Tier-1 fixes below.

## Levers, re-ranked by measured impact

### Tier 1 — cut the round-trip count (~3.6 s each; the only lever for 42–85 s help turns)

**1a. GroundingGate redesign** ([grounding.ts](../dockerfiles/claude-runtime/src/grounding.ts)).
Today: on Stop, if no `Read`/`Grep`/`Glob`/`Bash` whose input contains `'help/'` ran this turn, it
blocks the stop once — *after the model already wrote its answer* (a wasted full answer + a forced
grep chain). The prompt sends it grepping `workspace/help/` — **536 md files / 3.7 MB, no index** —
so it iterates Grep→Read→Grep, a full round-trip each. WebSearch/WebFetch and `js-api/src/` reads
do **not** count as grounded (a trap: legitimate grounding gets re-blocked). This mechanically
produced the 8–33-hop help turns and all the 120 s timeouts.
- ~~Ship a generated **index/TOC of `help/`**~~ — DONE (`help-index.ts`, 405 pages).
- ~~Count WebFetch/WebSearch/js-api reads as grounded~~ — DONE; `Skill` invocations too (2026-07-29).
- ~~Blunt post-answer Stop-block~~ — **narrowed 2026-07-29** after a 36-turn A/B showed 0 of 14
  blocks ever changed an answer (BENCHMARK.md "Third result"): small talk never arms the gate,
  and the block fires only when the visible answer makes ungrounded platform *UI-instruction*
  claims (`makesPlatformClaims`). Contextual data answers end in one API call (measured: 2 → 1
  calls, 7.6 → 3.6 s). Per-turn `gates: {grounding, verify}` overrides exist for future A/Bs.

**1b. Verifier scope-down** ([verify.ts](../dockerfiles/claude-runtime/src/verify.ts)).
Today: **every** `datagrok_exec` arms the verify gate — the exec's own result is never parsed, so
even an exec that already returned `{success:true, returnValue:…}` forces a separate
`datagrok_verify` round-trip (and a failing exec's error report is likewise ignored). Minimum
action turn = 2 round-trips — exactly the measured 11 s histogram. Cap is 3 Stop-blocks.
- Let a successful exec **self-verify**: return machine-checkable post-state from the same browser
  round-trip (or auto-run the assertion browser-side inside exec), keep the Stop-gate only when the
  exec returned nothing checkable.
- Verify once per turn, not per action.

**1c. Keep effort high** (measured, above). Revisit only after 1a/1b land — then re-benchmark
`medium`; with fewer forced hops it may win.

### Tier 2 — cut the ~12–14 s fixed intercept

**2a. Persistent streaming session** ([session.ts:331](../dockerfiles/claude-runtime/src/session.ts#L331)).
Every turn spawns a fresh SDK `query()` = fresh CLI subprocess + **transcript JSONL reload from disk
(grows with history)** + plugin scan + MCP re-handshake. The SDK's streaming-input mode keeps one
process per session and pushes messages into it. Refactor touches the turn queue / abort / fork in
`session.ts`. Est. ~1–2 s/turn plus the history-reload growth.

**2b. mcp-server handshake** ([mcp-server/src/index.ts:345](../dockerfiles/mcp-server/src/index.ts#L345)).
The HTTP `datagrok` MCP server gets `initialize` + `tools/list` **per turn**, through two proxies,
and constructs a fresh `McpServer` + transport per request. Falls out of 2a for free; otherwise
cache the server instance.

**2c. The rest of the intercept is model-side first-segment latency** (thinking + generation at
effort high before the first token/tool). The thinking budget (1500,
[query-options.ts:279](../dockerfiles/claude-runtime/src/query-options.ts#L279)) exists for the
grounding rule; after 1a lands, A/B a reduced/disabled budget **for action turns** via the harness.

### Tier 3 — browser-side (real, smaller; two cliffs)

**3a. Streaming render is O(n²)** ([panel.ts:658](../src/ai/panel.ts#L658) — path from repo:
`src/ai/panel.ts`). Every chunk — and every `tool_activity` — re-renders the **entire accumulated
text** through Dart markdown (parse + sanitizer + `querySelectorAll('*')` walk) and forces a reflow.
Hundreds of ms to seconds per help turn, on the same thread that processes WS messages — it delays
chunk/tool/final handling, inflating real turn time. Fix: append plain text while streaming, render
markdown **once** in `finalizeStreaming` (which already re-renders anyway), or throttle to rAF/150 ms.

**3b. Cliffs**: cold on-demand container → `webSocketProxy` waits up to **60 s**
(js-api `src/dapi.ts:1172`); reconnect is lazy (next send pays it) — auto-reconnect on `onClose`
instead. `runVerification` has a **30 s** hang cap (`src/claude/exec-blocks.ts:27`) — lower it.

**3c. Pre-send serial gauntlet** (`src/ai/ui.ts:276-284`) — `grok.ai.processPrompt` (seven full-DOM
`querySelectorAll` sweeps; purely local) then the usage limiter, serial before send: 10–150 ms. Run
`ensureConnected` in parallel from `handleRun` entry; skip/parallelize `processPrompt` for long
prompts (its DOM-click interpreter targets terse commands).

**Not browser problems** (measured): `buildWorkspaceContext` ≈ 75–90 tokens for demog, sub-ms (only
cap the **ScriptView full-source** inline, `src/claude/exec-blocks.ts:58-63`); usage limiter ~1 ms,
no network; WS double-hop through Datlas is pure passthrough (~ms per leg, only matters × 4 legs ×
many hops); IndexedDB save is post-final jank only.

### Tier 4 — context diet — LANDED 2026-07-28 (27,977 → 17,637 tok/call)

**Measured end to end on the suite** (`benchmark-diet-full` vs `benchmark-sonnet-clean`, 44 shared
prompts): median turn **26.6 s → 17.9 s (−33%)**, TTFT 17.9 s → 14.4 s, median cache-read per turn
92 k → 56 k, accuracy **39/44 in both arms** — the whole saving was free. This also re-tests the
Tier-1 model: cutting the per-hop prefix cuts every hop, so it moves the intercept *and* the slope.


Measured with `dev/harness/context-probe.mjs` (see BENCHMARK.md "Context probe"), which
attributes the prefix by ablation. Three levers landed the same day: restricting the *declared*
built-in tool set (`tools:` — `allowedTools` only pre-approves, so Task/TodoWrite/NotebookEdit
schemas shipped unused, 7,572 tok/call), consolidating the MCP surface (below), and exempting
small talk from the grounding gate's second API call.

- ~~**34** datagrok MCP tool schemas attach to every turn — filter via `allowedTools`~~ — **DONE**,
  differently: filtering by `allowedTools` would have made operations undiscoverable. Instead the 34
  tools became **five domain tools** dispatching on an `op`, whose descriptions carry only operation
  names; a call with no `op` returns the full parameter schema. **3,847 → 1,061 tok/call.**
- ~~`list_functions` & friends are unpaged and pretty-printed~~ — **DONE**: results are compact JSON,
  arrays paged at 50 with `offset`/`limit`, everything capped
  ([mcp-server/src/format.ts](../dockerfiles/mcp-server/src/format.ts)).
- **Still open — the plugin + system-prompt block is 6,363 tok/call**, the largest Datagrok-specific
  contributor, bigger than all tool definitions combined. Trimming skill descriptions is worth
  ~400 tok; the rest is plugin machinery and needs investigation before cutting.
- **Partly done — full mode costs a SECOND API call** (`numTurns` 2 vs 1), because the GroundingGate
  blocks the Stop once on every turn that neither grounded nor acted, and the model burns a whole
  call replying `NO_REVISION`. Small talk is now exempt (gate not armed at all — measured on
  "hello": 2 calls → 1, 9.3 s → 6.0 s). The remaining case — real answers that made no tool call —
  still pays the extra call; that is an accuracy feature and needs a `help`-category A/B (Tier 1a).
- Staged user workspace symlinks the repo's `CLAUDE.md`, `.claude/`, **`.kg/`** into the model's cwd
  ([staged-workspace.ts:22](../dockerfiles/claude-runtime/src/user/staged-workspace.ts#L22)). The
  repo CLAUDE.md instructs a `.kg/scripts/qq.py` query first (self-installs a venv, ~30 s) — a live
  timeout mechanism if the SDK version auto-loads subtree CLAUDE.mds. Exclude those entries; **pin
  the SDK** (currently `>=0.2.0`, [package.json](../dockerfiles/claude-runtime/package.json)).

## Confirmed non-issues / disproven earlier guesses

- ~~"Turn down effort — biggest single win"~~ — **disproven by A/B**: net loss while gates force hops.
- ~~Per-turn `datagrok-view` tool defs bust the prompt cache~~ — defs are module-level constants,
  byte-stable per session (`src/ai/view-tools.ts:131-160`). Cache works: measured input_tokens per
  turn = 4–8, cache-read 55–140 k.
- ~~`buildWorkspaceContext` needs diffing~~ — it's ~80 tokens; only the ScriptView case needs a cap.
- Workspace git sync, user-file sync, `ensureUserDir` — all off the hot path (fire-and-forget/memoized).

## Working method

Change one lever → rebuild the runtime image → `grok publish` →
`Grokky:runBenchmark("<lever>")` → `Grokky:compareBenchmarks("baseline", "<lever>")`.
Suggested order: **1a → 1b → re-A/B effort/thinking → 2a → 3a**, one benchmark row each.
Harness TODO: also aggregate `numTurns`/`durationApiMs` into the report (captured but not yet
summarized) to split model time from browser/exec time per turn.
