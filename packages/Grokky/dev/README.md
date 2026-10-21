# Grokky dev loop

Lets changes to `claude-runtime` be tested without an image rebuild and without a browser.

## One-time setup

```bash
npm --prefix dockerfiles/claude-runtime install
npm --prefix dev/harness install
node dev/runtime.mjs up                        # local container on ws://localhost:5355/ws
```

`runtime.mjs up` reads the Anthropic API key from your environment (`ANTHROPIC_API_KEY`). It hands
the key to the container by name, so it never appears on a command line.

## The loop

```bash
node dev/check.mjs            # compile + unit tests + a live turn   (~15s)
node dev/check.mjs --fast     # compile + unit tests only            (~4s, no model calls)
```

While iterating, `npm --prefix dockerfiles/claude-runtime run watch` keeps `dist/` current;
`node dev/runtime.mjs restart` (~3s) is still needed for the container to pick it up, because
the bind mount refreshes the files but not the modules Node already loaded. Only `dist/` and `plugin/`
are mounted: a change to the agent driver (`dockerfiles/claude-runtime/agent/`) needs an image rebuild.

## Test tiers

| Tier | Command | Container | Browser | Model | Speed |
|---|---|---|---|---|---|
| **T0** unit | `npm --prefix dockerfiles/claude-runtime test` | no | no | no | ~250ms |
| **T2** driver | `dev/harness/drive.mjs` | yes | no | yes | ~2–20s/turn |
| **T3** full E2E | `dev/harness/run-benchmark.mjs` | yes | yes | yes | minutes |

**T3** is [`harness/run-benchmark.mjs`](harness/run-benchmark.mjs) — Playwright drives a real
logged-in tab, because benchmark asserts read live `DG` objects that exist only in a browser.
Its companion [`harness/introspect.mjs`](harness/introspect.mjs) evaluates one expression in that
same tab, which is how you settle *"is the assert wrong or is the model wrong?"* against the real
API instead of by reading source:

```bash
node introspect.mjs "DG.Viewer.fromType(DG.VIEWER.BOX_PLOT, t).props.getProperties().map(p => p.name)"
node introspect.mjs "t.getCol('AGE').toList().filter(v => v > 60).length"   # verify a rubric constant
```

**T0** covers the pure logic: the verify/grounding gate state machines, the fence-aware stream
filter, help-index generation. Anything with no I/O belongs here — it is the tier that runs on
every edit.

**T2** is [`harness/drive.mjs`](harness/drive.mjs): a fake browser speaking the same WebSocket
protocol as `src/claude/runtime-client.ts`. It exercises the real runtime — system prompt,
tools, gates, revision protocol, session resume — and answers browser tool calls from a
per-case stub table. A turn reports which tools ran **in order**, browser round-trips, the code
each `datagrok_exec` asked to run, TTFT/total, and the runtime's token/cost metrics.

```js
import {RuntimeDriver} from './harness/drive.mjs';

const d = await new RuntimeDriver().connect();
const r = await d.turn('Add a histogram of AGE', {
  model: 'claude-sonnet-5-5',
  stubs: {datagrok_exec: () => ({success: true, verified: {passed: true}})},
});
// r.tools -> ['datagrok_exec'],  r.roundTrips -> 1,  r.execCode[0] -> the JS it ran
d.close();
```

Pass `sessionId` to continue a conversation (exercises the resume path), `model` to pin a model,
and `apiKey` + `mcpServerUrl` to the constructor when a case needs the real Datagrok MCP tools.

## Load testing — N concurrent users

[`harness/load-test.mjs`](harness/load-test.mjs) answers "what happens when N users hit the one
runtime container at once". Each virtual user is its own WS connection + sessionId (turns run
concurrently in the runtime exactly like distinct users — only same-session turns queue). All
users start together, then each loops: prompt from the pool → wait for the answer → think → next.
Every event (prompt sent, answer received with TTFT/total, errors, `docker stats` samples, a
5-second status line with the live in-flight count) is printed as it happens and appended to
`dev/harness/results/<label>.jsonl` — full answer texts included, so a run is directly observable.

```bash
cd dev/harness
node load-test.mjs --users 10                     # full Datagrok prompt, pool from suite.yaml
node load-test.mjs --users 100 --mode bash        # cheap echo turns (~$0.07/turn) — infra limits
node load-test.mjs --users 50 --turns 5 --think 0-0 --yes   # max sustained pressure
```

`--prompts <file>` swaps the pool (one prompt per line), `--seed` makes prompt choice and think
times reproducible, `--url` picks the runtime (a platform container works too — find its host
port with `docker port <name>`), `--container` names which container `docker stats` samples,
`--host <alias>` sends the config.yaml key as `apiKey` so per-user file sync is exercised too.
The summary reports peak observed concurrency (so you can verify the run really had N turns in
flight) and scrapes watchdog/reaper kills from the container log — correlating "client saw
timeout" with "watchdog fired at 95% CPU" is the point of the exercise.

## Notes

- The dev container is deliberately separate from the platform's spawned
  `datagrok-dev_cvm_grokky-claude-runtime-*` containers. Neither disturbs the other.
- `apiKey` and `mcpServerUrl` are optional in the wire protocol: without them the runtime skips
  user-file sync and the `datagrok` HTTP MCP server, keeping the browser and view tools, which
  run in the agent and round-trip over the WebSocket. That is what makes browserless turns possible.
- Outside the container, turns need a Hermes install, whose Python `GROKKY_AGENT_PYTHON` names;
  unit tests do not. See Test in the Operations section of `docs/CLAUDE_CODE_FLOW.md`.
