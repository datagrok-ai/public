# u2 browser checks

Two lanes over the same selector vocabulary. Design: `core/docs/features/ui2/TESTING.md`.

## Lane 1 — local mode (the iterate loop, seconds per run)

The real, unmodified client with **no server**: `?mode=local` boots from static files with a
constant token, so there is no login, no publish and no stand. Two commands:

```bash
grok-core local watch U2Demo        # re-stages dist/ on every webpack build (+ `npm run watch` in the package)
npm run e2e:local                   # from libraries/u2 — Playwright against http://localhost:63343
npm run e2e:local -- --only editing # one feature (substring of the check-file name); the results file is shared, last run wins
```

`grok-core local up` starts pub serve if nothing is serving yet — **only one pub serve may run**;
if `http://localhost:63343` already answers, use it, never start a second.

`run.mjs` checks both preconditions before opening a browser: pub serve answering, and the package
staged into `core/client/xamgle/web/local/pkg/` (it stages it if it is missing). A warm run is
~25 s including boot; results land in `.artifacts/` (gitignored) as `designer.local.json` plus a
screenshot per check.

Client-side function metadata comes from the captured fixture, so a package function added since the
last `grok-core local fixture` is not in `DG.Func` — `openApp` then loads the staged bundle and calls
the export directly, and says so as a NOTE. Refresh with
`grok-core local fixture --host <alias> --packages U2Demo` (needs a stand) to get the real path back.

## Lane 2 — a stand (the proof loop)

```bash
cd ../../packages/U2Demo && npm run debug-u2demo-<host>     # publish first
cd ../../libraries/u2 && npm run e2e:stand                  # U2_STAND_URL, default http://localhost:8888
```

Stand runs need a stand **nobody else is using**: they log in as admin and drive the shared UI.
The stand lane covers only what local mode cannot answer — that the *published* package registers
the app, the `/apps/...` URL route, and the Browse > Dev card. Behaviour checks belong in Lane 1.

## What lives where

| File | What |
|---|---|
| `local.mjs` | local-mode boot, app open, console/pageerror capture, screenshots, results |
| `lib.mjs` | the shared readers and drivers (`panel`, `selectRow`, `waitStatus`, `openSpec`, `newForm`, `dropControl`, the pickers, `platformDrag`, …) and the specs more than one file opens |
| `checks/<feature>.spec.mjs` | one check file per feature, each exporting `fixture` (its own starting state: `newForm` / `openSpec(sample)` / `reopenApp`) and `checks`; ids are `<feature>/<n> <title>`, numbered per file |
| `run.mjs` | preconditions + one browser + every check file in a fixed order (`leak` last — it closes every view) → `.artifacts/designer.local.json`, nonzero exit on failure; `--only <substring>` runs the matching files alone |
| `stand.mjs` | login and the platform's own routes |
| `designer.stand.spec.mjs` | the stand-only checks, and their runner |

Plain `.mjs` on purpose: no build step between editing a check and running it.

## Conventions

- Address a spec node by `[data-u2-name="…"]` — every named node's root carries it. Platform chrome
  (`.d4-ribbon`, `.grok-prop-panel`) and the designer's own classes are the only other selectors.
- The designer's **status bar** is the authority on what is selected; the context panel is the
  platform rendering that, and it runs one `grok.shell.o` assignment behind, so a check that reads
  the panel re-asserts the selection the way a user would (`selectRow`). Never assert on
  `grok.shell.o` directly.
- Adding a check never renumbers another file: ids restart per file and screenshots are
  `<feature>-<n>-<what>.png`. A check may rely on the order within its own file, never on another
  file's leftovers — that is what each file's `fixture` is for.
- Environment overrides: `U2_LOCAL_URL`, `U2_STAND_URL`, `U2_STAND_LOGIN`, `U2_STAND_PASSWORD`,
  `U2_HEADED=1` (watch the run in a visible browser).
