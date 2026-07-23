# Running ApiTests under Node (`package-test-node`)

`src/package-test-node.ts` is a **headless test runner** that connects to a running Datagrok
server over HTTP and executes the `Dapi:` test categories under Node — no browser, no
Puppeteer. It complements the browser-based `grok test` (driven by `src/package-test.ts`).

It runs the **TypeScript sources directly via [`tsx`](https://tsx.is)** — there is **no build
step** for the runner. Edit a `.ts` test and re-run.

## Prerequisites

1. `npm install` in this package (pulls in `tsx`, `yargs`, `datagrok-api`, `@datagrok-libraries/test`).
2. **Node ≥ 20.6** — the scripts use `--import` plus the `module.register()` hook API.
3. A **running Datagrok server** and a valid **dev key** for it.
4. **ApiTests must be published to that server.** The runner looks the package up via
   `grok.dapi.packages.filter('shortName = "ApiTests"')` and aborts with
   `ApiTests package should be installed.` if it is missing. If the server doesn't have it,
   publish once with `npm run build && grok publish <server>` (or `npm run debug`).
5. **`datagrok-api` must expose the `./datagrok` subpath** (`startDatagrok`). The local
   `js-api` does; if you consume the npm-published `datagrok-api`, you need a version that
   includes that `exports` entry — otherwise link local js-api (`grok link` /
   `npm link datagrok-api`).

## Run

```bash
# defaults: --apiUrl=http://localhost:8888/api --devKey=admin --concurrentRuns=1
npm run start-node

# override server / key, or scope to a category
npm run start-node -- --apiUrl=https://yourserver/api --devKey=<key>
npm run start-node -- -c "Dapi: users"

# concurrency stress sweep (heavy: runs the suite at concurrency 10,15,…,100)
npm run stress-node
```

| Script | What it does | Server? |
|--------|--------------|---------|
| `start-node` | Run all UI-independent suites once (`--mode functional`, the default) | yes |
| `stress-node` | Run the stressTest-marked suite repeatedly at increasing concurrency (`--mode stress --concurrencyRange=10-100 --step=5`) | yes |

The runner loads the UI-independent category folders (`nodeTestDirs` in
`package-test-node.ts`): `dapi`, `dataframe`, `functions`, `bitset`, `valuematcher`,
`property`, `stats`, `shell`. `grid/`, `widgets/`, `ai/`, `packages/`, `utils/` stay
browser-only. A test file that fails to load is reported and skipped — it doesn't
abort the run (but fails the exit code).

CLI flags (see `parseArgs` in `package-test-node.ts`): `--apiUrl`, `--devKey` (required),
`--mode functional|stress` (default `functional`; `stress` = only stressTest-marked tests),
`-c/--categories`, `--concurrentRuns`, `--concurrencyRange "1-10"`, `--step`, `-l/--loop`.

Results are written to `test-report.csv` (gitignored), **overwritten on every run**.

### Via the `grok stresstest` CLI

The same runner is wired into the `grok` CLI (`tools/bin/commands/stress-tests.ts`). Run it
**from this package folder**; it resolves the server/dev key from `~/.grok/config.yaml` (the
`--host` alias), publishes ApiTests to that server, then invokes the tsx runner:

```bash
grok stresstest --host <alias>                       # build + publish ApiTests, then run
grok stresstest --host <alias> --skip-build --skip-publish   # reuse the already-published package
grok stresstest --host <alias> --concurrency-range 10-100 --step 5
grok stresstest --host <alias> --concurrent-runs 4
grok stresstest --host <alias> --loop
```

Flags: `--host` (config alias), `--skip-build`, `--skip-publish`, `--concurrent-runs`,
`--concurrency-range "<start>-<end>"`, `--step`, `--loop`. Unlike the npm scripts, it handles
publishing for you, so it's the easiest way to run against a server that doesn't yet have
ApiTests deployed.

## How it works (and why)

The js-api and `@datagrok-libraries/test` sources are written for the browser (webpack), so a
few adjustments let Node load the same sources unchanged:

- **`tsx` + ESM, dynamic imports.** `startDatagrok` and the test files are loaded with dynamic
  `import()`, so Node resolves their named exports via `cjs-module-lexer` instead of failing at
  static-link time, and `tsx` resolves the extensionless relative imports the libraries emit.
- **`node-test-loader/`** — `register.mjs` (loaded via `--import`) installs:
  - an ESM `load` hook + a CJS `Module._extensions` patch that **stub asset imports**
    (`.css`, images, fonts) the browser sources `import`/`require`;
  - dayjs plugin extensions (`utc`, `advancedFormat`) the sources assume are present.
- **`src/test-package.ts`** — a shared `_package` holder. Both entries (`package-test.ts` for
  the browser, `package-test-node.ts` for Node) call `setTestPackage(...)`. The dapi test files
  read `_package` from here instead of from `package-test.ts`, so the Node runner does **not**
  transitively load the entire browser suite (which would run import-time browser/Dart side
  effects such as `cache.ts`'s top-level `grok.functions.register`).

## UI degradation under Node

The former known failures (`grok.shell.user`, `grok.dapi.groups.currentUserGroups`) are fixed —
the Node Dart bundle now registers `grok_User`, shell vars, custom events, and client build info.
Remaining browser-only APIs degrade gracefully: notifications (`grok.shell.info/warning/error`)
log to the console, DOM builders return inert elements, and Dart-backed UI (views, dialogs,
viewers, grid) throws a descriptive `DGNotSupportedError`. Tests that need them carry
`skipReason: 'NodeJS environment'`.
