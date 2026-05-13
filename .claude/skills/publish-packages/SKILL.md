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
   Expected: `name` starts with `@datagrok/`, `version` is SemVer-clean
   and ≥ `1.0.0`, no top-level `"beta"` key. CI rejects a stray `beta`
   key with `"Remove beta property — it is deprecated"`
   (`packages.yaml:617-620`).

2. **Bump the version (SemVer, ≥ 1.0.0 to actually publish).**
   ```bash
   npm version <patch|minor|major> --no-git-tag-version
   ```
   Expected: `package.json` version increments. Any first published
   version must be **≥ 1.0.0** — `0.x.x` is flagged `"beta"` in the
   workflow matrix and skipped with the notice `"Version <v> is under
   1.0.0 and is not going to be published"` (`DG-FACT-451`).

3. **Confirm no local-path entries remain in `dependencies`.**
   ```bash
   jq '.dependencies | with_entries(select(.value | startswith("../")))' \
      package.json
   ```
   Expected: `{}` (empty). A `"datagrok-api": "../../js-api"` under
   `dependencies` (not `devDependencies`) marks the build as
   `unpublished_deps` and the workflow refuses to publish
   (`DG-FACT-452`, `packages.yaml:84,144-147`). Local link in
   `devDependencies` is fine for the test step.

4. **Commit and push to `master`.**
   ```bash
   git add packages/<Pkg>/package.json packages/<Pkg>
   git commit -m "<Pkg>: <one-line summary> (vX.Y.Z)"
   git push origin master
   ```
   Expected: a `Packages` workflow run starts (`on: push: paths:
   packages/**`). The matrix probes
   `https://registry.npmjs.org/<name>/<version>`; missing version ⇒
   publish job runs (`DG-FACT-454`).

5. **Watch the run — Build → Test → grok check → Publish.**
   ```bash
   gh run list --workflow=packages.yaml --limit 5
   gh run watch                              # pick the latest run
   ```
   Expected: green build, green test (unless `skipCI` or `Tests`/`Meta`
   in name), green `grok check` (skipped only for Tests/Meta/beta —
   `DG-FACT-153`, `packages.yaml:624-627`), green `npm publish
   --access public` (`DG-FACT-152`, `packages.yaml:629-634`). On
   success Slack `#build` posts the published version.

6. **Re-trigger manually if the auto-run failed.**
   In the GitHub UI: **Actions → Packages → Run workflow** on branch
   `master`, set `packages` to a space-separated list (e.g.
   `Demo Tutorials`). Same shape applies to sibling workflows —
   `Libraries`, `Tools` (`datagrok-tools`), `JS-API` (`datagrok-api`)
   (`DG-FACT-150`, `DG-FACT-151`).

7. **Verify the new version is on npm.**
   ```bash
   curl -s https://registry.npmjs.org/<name>/<version> | jq -r .version
   npm view <name> versions --json | jq '.[-3:]'
   ```
   Expected: the first command returns `"<version>"` (not `null`); the
   second lists the new version among the most recent.

8. **(Private path) ship straight to a Datagrok server.**
   For forks, internal packages, or air-gapped instances — bypass CI:
   ```bash
   cd packages/<Pkg>
   webpack --mode=production
   grok publish <HOST> --release   # <HOST> = alias from ~/.grok/config.yaml
   ```
   Expected: server returns `200`; package appears under **Manage →
   Packages** on `<HOST>`. Drop `--release` for a per-developer debug
   deploy (default). `--debug` and `--release` are mutually exclusive
   (`DG-FACT-156`, `tools/bin/commands/publish.js:620-621`).

9. **(Private path) share with a group.**
   `--debug` deploys are per-developer only. Either set `canView` /
   `canEdit` in `package.json` and re-publish `--release`, or share
   via the UI (**Manage → Packages → right-click → Share**):
   ```json
   {"canEdit": ["Administrators"], "canView": ["All users"]}
   ```
   Expected: the listed groups see the package on next refresh
   (`DG-FACT-458`).

## Common failure modes

- **Version under 1.0.0 — workflow green, no publish.** Build/test
  pass, then `"Version <v> is under 1.0.0 and is not going to be
  published"` and the publish step skips (`DG-FACT-451`). Fix: bump
  to `1.0.0`+.
- **Same version already on npm — publish silently skipped.** The
  matrix probes `registry.npmjs.org/<name>/<current_version>`; a hit
  means "already published" (`DG-FACT-454`). Fix: increment.
- **Wrong branch — push did not publish.** Only `refs/heads/master`
  publishes; PRs and feature branches build/test only (`DG-FACT-453`).
  Fix: merge to `master`, or **Run workflow** on `master`.
- **Local-path dep in `dependencies` blocks publish.** Workflow emits
  `"Version <v> has unpublished dependencies and is not going to be
  published"` (`DG-FACT-452`, `packages.yaml:84,144-147`). Fix: pin to
  a published version or move the link to `devDependencies`.
- **`grok check` fails — publish blocked.** Runs for every package
  whose name lacks `Tests`/`Meta` and isn't beta-flagged. Fix: run
  `grok check` locally, fix what it reports (function header drift,
  manifest issues), commit, push.
- **`name` not in `@datagrok` scope.** `General checks` hard-fails
  with `"Package should be in '@datagrok' scope"`
  (`packages.yaml:613-616`). Fix: rename in `package.json` (and the
  `library` field in `webpack.config.js`).
- **Private publish only visible to developer.** `--debug` (default)
  deploys are per-developer. Fix: re-publish with `--release` and set
  `canView`/`canEdit`, or share via the UI.

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
