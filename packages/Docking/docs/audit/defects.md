# Docking — confirmed defects

Adversarially verified against the working tree (2026-09-09). Severity is user-facing. IDs DK-xx.

## High

- **DK-01 — Stale module globals break a grid's rendering after a second run.**
  `src/utils/constants.ts:7,9` export mutable `let` globals set by `processAutodockResults`
  (`src/utils/utils.ts:167,169`, success path only). The `grid.onCellRender` closure installed by
  `runAutodock` (`src/package.ts:153-168`) reads them live: run AutoDock on table A, then on table B —
  A's still-subscribed handler resolves B's column names against A's grid, `grid.col(...)` is null,
  and grid A throws a TypeError on **every cell paint** from then on. On a failed second run the crash
  moves to `addColorCoding` (`src/utils/utils.ts:122`, `isTextColorCoded` on null). Fresh-session
  failure is silent (empty-string guard, `src/package.ts:143-144`).

- **DK-02 — Container: concurrent dockings clobber each other's output capture.**
  `dockerfiles/autodock/autodock.py:37-50`: `run_process` opens bare `out.txt`/`err.txt` in the Flask
  process CWD while `/autodock/dock_ligand_list` (:270-273) docks ligands concurrently in a
  ThreadPoolExecutor — threads truncate, read, and delete the same two files, so stdout/stderr (and
  the `"Error" in output` return-code heuristic, :58) cross-attribute between ligands.
  Non-deterministic wrong results/errors per ligand.

- **DK-03 — Container: command injection through the user-supplied GPF.**
  `get_receptor_name` (autodock.py:117-120) extracts `\S+` from the caller's `autodock_gpf`; the value
  flows into `dlg_path`/`out_path` (:130-131, :191-192) and then into a `shell=True` pipeline
  (:108-114 via :36,:41). `Docking:dockLigandCached` (`src/package.ts:70-85`) forwards an arbitrary
  caller-supplied JSON body, so any authenticated platform user can inject shell syntax inside the
  container. `ligand_format`/`receptor_format` additionally reach filenames (no shell, but path
  traversal inside the container). Mitigated only by container isolation + platform auth.

- **DK-04 — Container: permanently poisoned grid cache after a partial preparation.**
  `prepare_grids` (autodock.py:84-100) creates the cache folder FIRST and gates all work on folder
  existence — a crash/kill between `os.makedirs` and autogrid completion leaves the folder present
  with no `.maps.fld`, so every later request for that receptor+gpf hash skips preparation and fails,
  until the container volume is wiped. The kill endpoint (DK-06) can cause exactly this state.

## Medium

- **DK-05 — Non-200 container responses produce unactionable errors.**
  `fetchAndCheck` (`src/utils/auto-dock-service.ts:215-226`): the status check computes an unused
  errMsg with the real handling commented out; `dockLigandCached` (`src/package.ts:81-84`) has no
  status check at all — a 502/HTML error page surfaces as "Unexpected token '<' … not valid JSON".

- **DK-06 — One user's Terminate kills everyone's dockings.**
  `kill_all_processes` (autodock.py:236-252) terminates all tracked processes with no requester
  scoping; exposed via `AutoDockService.terminate()` and the raw endpoint. (No TS call site today.)

- **DK-07 — `onCellRender` work per painted cell + handler accumulation.**
  `src/package.ts:153-156`: `grid.setOptions({rowHeight})` (a Dart interop with JSON.stringify per
  call) plus two column-width writes run on every rendered cell; the subscription is never disposed
  and a new one is added per `runAutodock` invocation.

- **DK-08 — The `poses` input doesn't affect the output.**
  `runAutoDock` keeps only the min-energy row per ligand (`src/apps/auto-dock-app.ts:175-178`), so the
  user-facing "Number of output conformations" (`src/package.ts:97,135`) only multiplies container
  compute time. Misleading input, silent cost.

- **DK-09 [resolved] — Integration tests silently pass when the service is missing.**
  `src/tests/autodock-tests.ts:20-33`: `before()` swallows init errors, every test early-returns on
  `!adSvc` — 3/3 vacuous green with zero coverage; `ensureContainerRunning` sits after the guard.

- **DK-10 — Panel-open jank: full-column PDB parsing per widget build.**
  `getRemarksFromPdbs` (`src/utils/utils.ts:17-31`) regex-parses every pose in the column and builds an
  N-row frame of which one row is read, synchronously, on every current-cell change
  (`src/package.ts:205`).

- **DK-11 — Python 2 on ubuntu:20.04 in the container.**
  `unicode` builtins (autodock.py:23,30); `Dockerfile:35` `python=2.7` — EOL interpreter serving HTTP.
  Also `processes` dict grows unboundedly (entries removed only by the kill route).

## Low

- **DK-12 [resolved]** — `buildDefaultAutodockGpf` z-dimension uses `npts.y` (`auto-dock-service.ts:48`);
  masked today (all callers cubic; the generated GPF is a dead fallback — `gpfFile ?? autodockGpf`
  with gpfFile always set).
- **DK-13 [resolved]** — Dead exports with latent bugs: `_runAutodock` (`svc.ready` before `init()` →
  TypeError instead of the intended warning, :249-263), `_runAutodock2` (`ET: ${t2 - t2}` always 0,
  `window.alert`, forced download, :265-297) — imported at `package.ts:13`, called nowhere.
- **DK-14 [resolved]** — `autoDockApp` is a registered no-op: loads sample data via `Chem:importSdf` and shows
  nothing (`setRibbonPanels`/`downloadPosesBtn`/`posesGrid` never wired; `auto-dock-app.ts:84-99`,
  `package.ts:42-50`). Not in the app browser (plain func), so impact is wasted I/O.
- **DK-15 [resolved]** — Silent empty receptor: `fetchPdbContent` empty catch → `''` → empty Biostructure viewer
  with no message (`utils.ts:55-68` → `package.ts:213`).
- **DK-16 [resolved]** — `prop()` renames the column before the description lookup — second add gets no
  description; also mutates the column name on icon click (`utils.ts:70-86`).
- **DK-17 [resolved]** — `'Unbound System\s Energy'` → renders "Unbound Systems Energy" (`constants.ts:28`).
- **DK-18 [resolved]** — 11 of 13 registered functions lack descriptions (worst in the audited set);
  6 narration-comment blocks + stale commented-out code (`package.ts:110-121, 252-256`,
  `utils.ts:131-141`, `auto-dock-service.ts:146-151,162-163,220-221`, `demo.ts:22,29-34`);
  `@ts-ignore` at `auto-dock-service.ts:216`; hardcoded `'black'` fill in the error-cell renderer.
  Resolution: descriptions added to every registered function except `info`; stale commented-out code and
  narration blocks removed (two whys kept as one-liners); `@ts-ignore` dropped. The canvas `'black'` fill is
  kept: canvas cannot read CSS tokens and no package in the repo does otherwise.
- **DK-19 [resolved]** — `getAutodockSingle` evaluates `grok.shell.tv` even when `table` is passed; without a
  table view the failure surfaces later as an obscure `.plot` throw (`package.ts:198-212`).
