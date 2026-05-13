---
name: publish-packages
version: 0.1.0
description: |
  Ship a new version of a Datagrok plugin so other users can install
  it — either through the public-repo CI (NPM via the `Packages`
  workflow) or directly to a Datagrok server with `grok publish`.
  For plugin authors who finished local development and need to make
  the build reach teammates or the public platform. Produces a new
  `@datagrok/<name>` version on npm, OR a server-side package record
  on the target host, plus group-level sharing for the latter.
  Use when asked to "ship a plugin update to npm", "release a new
  plugin version to teammates", or "send my plugin build to an
  internal Datagrok instance".
triggers:
  - ship plugin update to npm
  - release plugin version to teammates
  - deploy plugin to internal datagrok
  - roll out a plugin update
  - make plugin available to other users
  - send build to datagrok server
allowed-tools:
  - Read
  - Bash
  - Edit
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# publish-packages

## When to use

You have a finished change in `packages/<Pkg>/` and want users to see
it. Two paths: (a) **public** — bump the version, push to `master`,
let `.github/workflows/packages.yaml` build/test/publish to npm;
(b) **private** — call `grok publish <HOST>` against an internal or
air-gapped Datagrok instance.

## Prerequisites

- A package under `packages/<Pkg>/` of the public repo (or a fork).
  `package.json` `name` MUST be `@datagrok/<lowercase>` —
  `General checks` hard-fails otherwise (`DG-FACT-450`,
  `packages.yaml:613-616`).
- **Public path:** write access to `datagrok-ai/public`'s `master`.
  Other branches and PR events skip publish silently (`DG-FACT-453`).
- **Private path:** `datagrok-tools` installed (`npm i -g
  datagrok-tools`) and `~/.grok/config.yaml` populated via `grok config`
  for the target host (`DG-FACT-456`; see `set-up-environment.md`).
- Locally green tests (`grok test --host <alias>`) — what CI re-runs
  anyway (`test-packages.md#local-testing`).

## Steps

1. **Verify the `package.json` invariants before you bump.**
   ```bash
   cd packages/<Pkg>
   jq -r '.name, .version' package.json
   jq -e 'has("beta")' package.json     # MUST be false (deprecated)
   jq -e '.skipCI // false' package.json
   ```
   `name` starts with `@datagrok/`; no top-level `"beta"` key (CI
   rejects it).

2. **Bump the version (SemVer, ≥ 1.0.0 to actually publish).**
   ```bash
   npm version <patch|minor|major> --no-git-tag-version
   ```
   `0.x.x` is flagged beta and silently skipped (see `DG-FACT-451`).

3. **Confirm no local-path entries remain in `dependencies`.**
   ```bash
   jq '.dependencies | with_entries(select(.value | startswith("../")))' \
      package.json
   ```
   Must be `{}`. A `"../"` entry under `dependencies` blocks publish
   (see `DG-FACT-452`). Local link in `devDependencies` is fine.

4. **Commit and push to `master`.**
   ```bash
   git add packages/<Pkg>/package.json packages/<Pkg>
   git commit -m "<Pkg>: <one-line summary> (vX.Y.Z)"
   git push origin master
   ```
   The `Packages` workflow probes npm; missing version ⇒ publish job
   runs (see `DG-FACT-454`).

5. **Watch the run — Build → Test → grok check → Publish.**
   ```bash
   gh run list --workflow=packages.yaml --limit 5
   gh run watch                              # pick the latest run
   ```
   `grok check` runs except for Tests/Meta/beta (see `DG-FACT-153`);
   publish is `npm publish --access public` (see `DG-FACT-152`).

6. **Re-trigger manually if the auto-run failed.**
   GitHub UI: **Actions → Packages → Run workflow** on `master`. Same
   shape for sibling workflows `Libraries` / `Tools` / `JS-API` (see
   `DG-FACT-150`, `DG-FACT-151`).

7. **Verify the new version is on npm.**
   ```bash
   curl -s https://registry.npmjs.org/<name>/<version> | jq -r .version
   npm view <name> versions --json | jq '.[-3:]'
   ```

8. **(Private path) ship straight to a Datagrok server.**
   For forks or air-gapped instances — bypass CI:
   ```bash
   cd packages/<Pkg>
   webpack --mode=production
   grok publish <HOST> --release   # <HOST> = alias from ~/.grok/config.yaml
   ```
   Drop `--release` for per-developer debug (default). `--debug` and
   `--release` are mutually exclusive (see `DG-FACT-156`).

9. **(Private path) share with a group.** Set `canView` / `canEdit` in
   `package.json` and re-publish `--release`, or share via UI
   (**Manage → Packages → right-click → Share**):
   ```json
   {"canEdit": ["Administrators"], "canView": ["All users"]}
   ```
   Groups see the package on next refresh (see `DG-FACT-458`).

## Common failure modes

- Version `<1.0.0` — workflow green, publish skipped (`DG-FACT-451`).
- Version already on npm — silently skipped (`DG-FACT-454`); increment.
- Wrong branch — only `master` publishes (`DG-FACT-453`).
- Local-path dep in `dependencies` blocks publish (`DG-FACT-452`); move
  to `devDependencies` or pin to a published version.
- `grok check` fails — fix function header drift / manifest issues
  locally first.
- `name` not in `@datagrok` scope — rename in `package.json` and the
  `webpack.config.js` `library` field (`DG-FACT-450`).
- Private publish only visible to developer — re-publish `--release`
  with `canView`/`canEdit`, or share via UI.

## See also

- Source articles: `help/develop/how-to/packages/publish-packages.md`
  (primary), `help/develop/how-to/tests/test-packages.md#local-testing`
  (`grok test --host <alias>`), `help/develop/develop.md#publishing`
  (`--debug` vs `--release`, `canEdit`/`canView`, transitive-dep
  warning at lines `448-454`).
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` —
  `DG-FACT-450` (`@datagrok` scope), `DG-FACT-451..454` (workflow
  gates), `DG-FACT-455` (workflow triggers + siblings),
  `DG-FACT-456..458` (`grok publish <HOST>`, build→test→publish,
  group sharing).
- Reference packages: `packages/Chem/package.json` (canonical scoped
  layout); `packages/ApiTests/package.json` (intentionally pins
  `"datagrok-api": "../../js-api"`, never published).
- Workflow: `.github/workflows/packages.yaml` — gates at lines `84`,
  `128-159`, `613-635`.
- Related skills: `manage-credentials`, `add-package-tests`.
