# Grokky dev loop

Lets changes to `claude-runtime` be tested without an image rebuild, without a browser, and
without an Anthropic API key — turns run on the developer's **Claude subscription**.

## One-time setup

```bash
node dev/seed-creds.mjs                        # copy subscription creds (see below)
npm --prefix dockerfiles/claude-runtime install
npm --prefix dev/harness install
node dev/runtime.mjs up                        # local container on ws://localhost:5355/ws
```

## The loop

```bash
node dev/check.mjs            # compile + unit tests + a live turn   (~15s)
node dev/check.mjs --fast     # compile + unit tests only            (~4s, no quota spent)
```

While iterating, `npm --prefix dockerfiles/claude-runtime run watch` keeps `dist/` current;
`node dev/runtime.mjs restart` (~3s) is still needed for the container to pick it up, because
the bind mount refreshes the files but not the modules Node already loaded.

## Credentials

`seed-creds.mjs` copies **only** the `claudeAiOauth` block out of `~/.claude/.credentials.json`
into `dev/.claude/` (gitignored, mode 600). Two reasons it is a copy, not a mount of the real file:

- that file also holds unrelated MCP OAuth tokens (Slack, …) which the container has no business seeing;
- the CLI refreshes the access token in place, so mounting the live file would make the container
  and your own Claude Code session share one credential, where refresh-token rotation can log one
  of them out mid-session.

The copy is mounted read-write so the container renews its own access token. Re-run the script
when the **refresh** token expires (it prints both expiry times).

### ⚠️ The platform container loses its credentials whenever it is recreated

The Dockerfile seeds `~/.claude/.credentials.json` with a literal `{}`. On a subscription, real
credentials only ever arrive through the in-app `claude auth login`, which writes into the
**running container's** filesystem — not the image, not a volume. So *anything* that recreates the
container resets it: rebuilding the image, `docker rm`, or the platform respawning an on-demand
container after it idles out.

It does not fail loudly. Every full-prompt turn returns in ~1 s with **zero tokens and no error**,
which looks exactly like "the model answered instantly and did nothing" — i.e. like a model or
prompt regression, which it is not.

```bash
node dev/seed-platform-creds.mjs --check   # is every running platform runtime authenticated?
node dev/seed-platform-creds.mjs           # seed the empty ones from dev/.claude/
```

Run `--check` before any benchmark: an unauthenticated container produces a full report of
plausible-looking failures. (`runBenchmark`'s readiness probe now catches the dead case, but check
anyway — it is one second and the failure mode is expensive.)

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
MCP servers, gates, revision protocol, session resume — and answers browser tool calls from a
per-case stub table. A turn reports which tools ran **in order**, browser round-trips, the code
each `datagrok_exec` asked to run, TTFT/total, and the SDK's token/cost metrics.

```js
import {RuntimeDriver} from './harness/drive.mjs';

const d = await new RuntimeDriver().connect();
const r = await d.turn('Add a histogram of AGE', {
  model: 'sonnet',
  stubs: {datagrok_exec: () => ({success: true, verified: {passed: true}})},
});
// r.tools -> ['datagrok_exec'],  r.roundTrips -> 1,  r.execCode[0] -> the JS it ran
d.close();
```

Pass `sessionId` to continue a conversation (exercises the resume path), `model` to pin a tier,
and `apiKey` + `mcpServerUrl` to the constructor when a case needs the real Datagrok MCP tools.

## Notes

- The dev container is deliberately separate from the platform's spawned
  `datagrok-dev_cvm_grokky-claude-runtime-*` containers. Neither disturbs the other.
- `apiKey` and `mcpServerUrl` are optional in the wire protocol: without them the runtime skips
  user-file sync and the `datagrok` HTTP MCP server, keeping the in-process browser and view
  tool servers. That is what makes browserless turns possible.
- The container runs Node 22; a typical dev host may be on Node 20, where `require()` of the
  ESM-only Agent SDK throws. Keep unit-testable logic in modules that do not import the SDK
  (this is why the stream filter lives in `src/stream-filter.ts`).
