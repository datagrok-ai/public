---
name: create-package
version: 0.1.0
description: |
  Scaffold a brand-new Datagrok plugin folder, wire the standard
  `build` / `test` / `publish` scripts, and ship the first version
  to a developer-private slot on a Datagrok host. For first-time
  plugin contributors and anyone starting a project from zero who
  needs a template that imports `datagrok-api`, regenerates
  `package.g.ts` from annotations, validates with `grok check`, and
  bundles with webpack. Produces a `<Name>/` directory containing
  `src/package.ts`, `webpack.config.js`, a configured `package.json`,
  and a working `grok publish <host>` round-trip.
  Use when asked to "scaffold a fresh plugin project", "bootstrap a
  new plugin from zero", or "start a brand-new Datagrok project".
triggers:
  - scaffold a plugin from scratch
  - bootstrap a new plugin project
  - start a fresh plugin from zero
  - generate a starter project skeleton
  - first plugin folder setup
  - new datagrok project template
allowed-tools:
  - Read
  - Write
  - Edit
  - Bash
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# create-package

## When to use

You're starting a brand-new Datagrok plugin and need the canonical
folder layout — `src/package.ts`, `webpack.config.js`, `package.json`
with the `grok api && grok check --soft && webpack` build triplet,
and a `grok publish` round-trip — before you write any business code.

## Prerequisites

- An empty target directory: `grok create` aborts with *"The package
  directory should be empty"* otherwise
  (`tools/bin/commands/create.ts:162-164`).

## Steps

1. **Scaffold from the package template.** Run from the directory
   where the new package folder should appear:
   ```bash
   grok create TextStats --test
   ```
   The CLI auto-runs `npm install` inside the new folder (see
   `DG-FACT-463`). `--test` only pulls `@datagrok-libraries/utils`;
   no test files are scaffolded — use `grok add tests` for those
   (see `DG-FACT-464`).

2. **Fix folder casing and `friendlyName` if you want PascalCase /
   Title Case.** `grok create TextStats` writes folder `Textstats/`
   with `friendlyName = "TextStats"`; the CLI normalizes everything
   after the first letter (see `DG-FACT-441`, `DG-FACT-DRIFT-CP-001`).
   ```bash
   mv Textstats TextStats   # only if PascalCase matters to you
   # then edit package.json: "friendlyName": "Text Stats"
   ```

3. **Add a function — header-annotation form.** Drop a panel function
   into `src/package.ts` to confirm wiring end-to-end:
   ```typescript
   //name: TextStats
   //meta.role: panel, widgets
   //input: string str {semType: text}
   //output: widget result
   export function textStats(str: string) {
     const counts = Array.from(str).reduce((acc, ch) => {
       acc[ch] = (acc[ch] || 0) + 1;
       return acc;
     }, Object.create(null));
     return new DG.Widget(ui.divV([
       ui.divText('Counting characters:'),
       ui.divText(JSON.stringify(counts)),
     ]));
   }
   ```
   Header annotations sit directly above `export function`. The
   panel fires for cells whose column `quality` tag equals `text`
   (`quality` IS the semantic-type tag — see `DG-FACT-440`).

4. **Build — do not paraphrase to bare `webpack`.** `npm run build`
   runs the load-bearing triplet `grok api && grok check --soft &&
   webpack` (see `DG-FACT-462`). Skipping the prefix leaves new
   annotations out of the bundle.
   ```bash
   npm run build
   ```
   After: `src/package.g.ts` lists `textStats` with
   `//meta.role: widgets,panel`.

5. **Publish to a dev host in debug mode.** `--debug` is the default;
   visible only to your developer key (see `DG-FACT-156`).
   ```bash
   grok publish dev
   ```
   The package shows up under **Manage → Packages** prefixed
   `v.<your-name>` — that prefix IS the debug-mode marker.

6. **Promote to release when ready.** A release-mode publish makes
   the package visible to every group listed in `canView` /
   `canEdit` in `package.json`.
   ```bash
   grok publish dev --release
   ```
   The `v.` prefix disappears. See `publish-packages` for the full
   ship pipeline.

## Common failure modes

- `grok` not found — `datagrok-tools` not on `PATH`. Install per
  `set-up-environment.md:21`.
- `grok create` aborts with "package directory should be empty" — cd
  to a fresh directory.
- Folder came out lowercase (`Textstats`) — expected, see step 2.
- Build succeeds but panel never registers — you ran `webpack`
  directly; `src/package.g.ts` is stale (`DG-FACT-462`).
- `grok publish` errors "Incompatible options: --debug and --release"
  — flags are mutually exclusive (`DG-FACT-156`); pick one.

## See also

- Source article: `help/develop/how-to/packages/create-package.md`.
- Adjacent articles: `help/develop/dev-process/set-up-environment.md`
  (prerequisites); `help/develop/develop.md#publishing` (`--debug`
  vs `--release`, `canView` / `canEdit`).
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` —
  `DG-FACT-440` (`quality` tag is the semType tag),
  `DG-FACT-441` (name normalization),
  `DG-FACT-462` (build triplet),
  `DG-FACT-463` (auto `npm install`),
  `DG-FACT-464` (`--test` is narrow),
  `DG-FACT-156` (`grok publish` flags),
  plus drifts `DG-FACT-DRIFT-CP-001..004`.
- Reference packages: `tools/package-template/package.json` (source
  template); `packages/Chem/package.json:99,106` (production
  preserves the build triplet); `packages/BenchlingLink`,
  `packages/ChemblAPI` (PascalCase folders preserved by hand-rename).
- Related skills: `add-info-panel` (step 3's function is a panel);
  `add-package-tests` (`grok add tests`, required follow-up for real
  tests); `publish-packages` (full `Packages` workflow + npm publish);
  `define-semantic-type-detectors` (programmatic alternative to the
  manual `quality : text` tag).
