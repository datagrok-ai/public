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
   Expected: a sibling folder is created and the CLI immediately
   runs `npm install` inside it — `Running 'npm install' to get the
   required dependencies...` comes from the CLI, not from you
   (`tools/bin/commands/create.ts:201-205`, `DG-FACT-463`). A manual
   `npm install` is only needed if the auto-install printed errors.
   `--test` adds `@datagrok-libraries/utils` to `dependencies` but
   does NOT scaffold test files — run `grok add tests` later if you
   want them (`tools/bin/commands/create.ts:80-87`, `DG-FACT-464`).

2. **Fix folder casing and `friendlyName` if you want PascalCase /
   Title Case.** `grok create` normalizes the raw arg to
   first-letter-upper + rest-lower for the directory and
   `#{PACKAGE_NAME}` substitution (`tools/bin/commands/create.ts:150`,
   `DG-FACT-441`) — `grok create TextStats` writes folder
   `Textstats/`. `friendlyName` is identical to the raw arg; there
   is no `--friendly-name` flag. Hand-edit both if you want the
   article's convention (`create-package.md:25-29` describes the
   convention, not what the CLI produces;
   `DG-FACT-DRIFT-CP-001`).
   ```bash
   mv Textstats TextStats   # only if PascalCase matters to you
   # then edit package.json: "friendlyName": "Text Stats"
   ```
   Expected: directory name + `jq -r .friendlyName package.json`
   match what you typed.

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
   Expected: header-annotation comments sit immediately above
   `export function`. At runtime the panel fires for cells whose
   column's `quality` tag equals `text` — `quality` IS the
   semantic-type tag (`DG-FACT-440`, `js-api/src/api/ddt.api.g.ts:77`).

4. **Build — do not paraphrase to bare `webpack`.** `npm run build`
   runs the template-defined triplet `grok api && grok check --soft
   && webpack` (`tools/package-template/package.json:21`,
   `DG-FACT-462`). `grok api` regenerates `src/package.g.ts` from
   header annotations; `grok check --soft` warns on package-level
   issues; `webpack` bundles. Skipping the prefix leaves new
   `//name:` / `//meta.role:` annotations out of the bundle.
   ```bash
   npm run build
   ```
   Expected: exit 0; `src/package.g.ts` is rewritten and lists
   `textStats` with `//meta.role: widgets,panel`.

5. **Publish to a dev host in debug mode.** `--debug` is the default;
   the package is visible only to your developer key
   (`tools/bin/commands/publish.js:620-621`, `DG-FACT-156`).
   ```bash
   grok publish dev
   ```
   Expected: exit 0. The package appears under **Manage → Packages**
   on the `dev` host prefixed `v.<your-name>` — that prefix IS the
   debug-mode marker.

6. **Promote to release when ready.** A release-mode publish makes
   the package visible to every group listed in `canView` /
   `canEdit` in `package.json`. See the `publish-packages` skill for
   the full ship pipeline (npm publish via the `Packages` workflow).
   ```bash
   grok publish dev --release
   ```
   Expected: the `v.` prefix disappears; the package is listed for
   the configured groups.

## Common failure modes

- **`grok` not found.** `datagrok-tools` isn't installed globally, or
  the npm global `bin` directory isn't on `PATH`. Fix: install per
  `set-up-environment.md:21`, re-open the shell.
- **`grok create` aborts: "package directory should be empty".** The
  CWD already contains files (`tools/bin/commands/create.ts:162-164`).
  Fix: `cd` to a fresh directory, or pass a name that creates a new
  subfolder.
- **Folder name came out lowercase (`Textstats`).** Expected — the
  CLI capitalizes only the first letter (`DG-FACT-441`). Rename by
  hand if you need PascalCase.
- **`npm run build` succeeds but the panel never registers.** You
  ran `webpack` directly, so `grok api` never ran and
  `src/package.g.ts` is stale (`DG-FACT-462`). Fix: always use
  `npm run build`; verify the new function appears in
  `src/package.g.ts` after the build.
- **`grok publish` exits non-zero: "Incompatible options: --debug
  and --release".** Both flags were passed; they are mutually
  exclusive (`DG-FACT-156`). Fix: drop one — `--release` for shared
  visibility, no flag for per-developer debug.

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
