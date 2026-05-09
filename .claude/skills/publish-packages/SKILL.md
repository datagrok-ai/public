---
name: publish-packages
description: Ship a new version of a Datagrok package — either via NPM through the public-repo CI, or directly to a Datagrok server with `grok publish`
---

# publish-packages

## When to use

You have a finished change in `packages/<Pkg>/` and want users to see it.
Two paths: (a) **public** — bump the version, push to `master`, let
`.github/workflows/packages.yaml` build/test/publish to npm; (b)
**private** — call `grok publish <HOST>` against an air-gapped or
internal Datagrok instance. The article calls these "Public packages"
and "Private packages."

## Prerequisites

- A package under `packages/<Pkg>/` of the public repo (or a fork).
  `package.json` `name` SHOULD be `@datagrok/<lowercase>` (knowledge
  `DG-FACT-144`; four legacy packages drift — `DG-FACT-DRIFT-056`).
- **Public path:** write access to `datagrok-ai/public`'s `master` —
  other branches and PR events skip publish silently (`DG-FACT-146`,
  `DG-FACT-148`).
- **Private path:** `datagrok-tools` installed (`npm i -g
  datagrok-tools`) and `~/.grok/config.yaml` populated via `grok config`
  for the target host (`DG-FACT-155`; `set-up-environment.md`).
- Tested locally with `grok test --host <alias>` — what CI runs anyway
  (`test-packages.md#local-testing`).

## Steps

1. **Verify the package.json invariants before you bump.**
   ```bash
   cd packages/<Pkg>
   jq -r '.name, .version' package.json     # name, version
   jq -e '.skipCI // false' package.json     # true ⇒ tests skipped in CI
   jq -e 'has("beta")' package.json          # MUST be false (deprecated)
   ```
   Expected: `name` begins with `@datagrok/`; `version` is SemVer-clean
   and ≥`1.0.0`; no top-level `"beta"` key (CI rejects it,
   `packages.yaml:617-620`); know whether `skipCI` is set
   (knowledge `DG-FACT-154`, `DG-FACT-DRIFT-058`).

2. **Bump the version (SemVer, ≥1.0.0 to actually publish).**
   ```bash
   npm version <patch|minor|major> --no-git-tag-version
   ```
   Expected: `package.json` version increments. The first npm-published
   version of any package must be **≥ 1.0.0** — `0.x.x` is flagged
   `"beta"` in the workflow matrix and never reaches npm with no error
   surfaced (knowledge `DG-FACT-149`, drift `DG-FACT-DRIFT-057`).

3. **Confirm no local-path dependencies remain in `dependencies`.**
   ```bash
   jq '.dependencies | with_entries(select(.value | startswith("../")))' \
      package.json
   ```
   Expected: `{}` (empty). A `"datagrok-api": "../../js-api"` or
   `"@datagrok-libraries/X": "../../libraries/X"` in `dependencies`
   (not `devDependencies`) marks the build as `unpublished_deps` and
   the workflow refuses to publish — see ApiTests/PowerPack/Grokky/
   KnimeLink for the four packages that intentionally hold this state
   (workflow `packages.yaml:84,144-147`).

4. **Commit and push to `master`.**
   ```bash
   git add packages/<Pkg>/package.json packages/<Pkg>/...
   git commit -m "<Pkg>: <one-line summary> (vX.Y.Z)"
   git push origin master
   ```
   Expected: a `Packages` workflow run starts. Trigger is
   `on: push: paths: packages/**` (knowledge `DG-FACT-146`). The
   matrix step diffs changed files, picks up `<Pkg>`, and probes
   `https://registry.npmjs.org/<name>/<version>`; missing version ⇒
   publish job runs (knowledge `DG-FACT-147`).

5. **Watch the run; expect Build → Test → grok check → Publish.**
   ```bash
   gh run list --workflow=packages.yml --limit 5
   gh run watch                              # pick the latest run
   ```
   Expected: green build, green test (unless `skipCI` or `Tests`/`Meta`
   in name), green `grok check` (skipped only for Tests/Meta/beta —
   knowledge `DG-FACT-153`, drift `DG-FACT-DRIFT-059`), green
   `npm publish --access public` (knowledge `DG-FACT-152`). On success
   a Slack `#build` message is posted and `package-lock.json` is
   committed back by `github-actions[bot]`.

6. **Re-trigger manually if the auto-run failed.**
   In the GitHub UI: **Actions → Packages → Run workflow** on branch
   `master`, set `packages` to a space-separated list (e.g.
   `Demo Tutorials`). Same workflow accepts the same input shape for
   sibling pipelines `Libraries`, `Tools` (`datagrok-tools`),
   `JS-API` (`datagrok-api`) (knowledge `DG-FACT-150`, `DG-FACT-151`).

7. **Verify the new version is on npm.**
   ```bash
   curl -s https://registry.npmjs.org/<name>/<version> | jq .version
   npm view <name> versions --json | jq '.[-3:]'
   ```
   Expected: the first command returns `"<version>"`; the second lists
   the new version among the most recent.

8. **(Private path) publish straight to a Datagrok server.**
   When npm is not in the picture — fork, internal package, air-gapped
   instance — bypass CI:
   ```bash
   cd packages/<Pkg>
   webpack --mode=production
   grok publish <HOST> --release   # <HOST> = alias from ~/.grok/config.yaml
   ```
   Expected: server returns `200`; package appears under **Manage →
   Packages** on `<HOST>`. Drop `--release` for a per-developer debug
   deploy (default; `DG-FACT-156`). `--debug` and `--release` are
   mutually exclusive.

9. **(Private path) share with a group.**
   `--debug` deploys are per-developer. Either set `canView`/`canEdit`
   in `package.json` and re-publish `--release`, or share via the UI:
   ```json
   {"canEdit": ["Administrators"], "canView": ["All users"]}
   ```
   Expected: listed groups see the package on next refresh.

## Common failure modes

- **Version under 1.0.0 — workflow runs, no publish, no error.** Build
  and test pass, then a notice "Version <v> is under 1.0.0 and is not
  going to be published" appears and the publish job is skipped
  (knowledge `DG-FACT-149`, drift `DG-FACT-DRIFT-057`). Fix: bump to
  `1.0.0` or higher; the very first published version is the floor.
- **Same version already on npm — publish step silently skipped.** The
  matrix probes `registry.npmjs.org/<name>/<current_version>`; a hit
  means "already published, nothing to do" (knowledge `DG-FACT-147`).
  Fix: increment the version before pushing.
- **Wrong branch — push didn't publish.** `Packages` only publishes
  from `refs/heads/master`; PRs and feature branches build/test only
  (knowledge `DG-FACT-148`, `DG-FACT-146`). Fix: merge to master, or
  trigger manually via `Run workflow` on `master`.
- **Local-path dep in `dependencies` blocks publish.** Workflow marks
  `unpublished_deps` non-empty and emits "Version <v> has unpublished
  dependencies and is not going to be published"
  (`packages.yaml:84,144-147`). Fix: pin to a published version, or move
  the local link to `devDependencies`. The *transitive* variant
  (`develop.md:284-291`) — a published library whose own deps point at
  `../../js-api` — trips the same gate.
- **`grok check` fails — publish blocked.** Runs for every package
  whose name doesn't contain `Tests` or `Meta` and isn't beta-flagged
  (knowledge `DG-FACT-153`, drift `DG-FACT-DRIFT-059`). Fix: run
  `grok check` locally, fix what it reports (function header issues,
  manifest drift), commit, push.
- **`name` not in `@datagrok` scope.** The `General checks` step
  hard-fails with "Package should be in '@datagrok' scope. Change
  package name to '@datagrok/<name>'" (`packages.yaml:613-616`). Fix:
  rename in `package.json` and `webpack.config.js` `library` field.
- **Top-level `"beta": true` left in `package.json`.** Hard-fails
  `General checks` with "Remove beta property — it is deprecated"
  (`packages.yaml:617-620`). Fix: delete the key; the workflow infers
  beta from version <1.0.0 instead.
- **Private publish succeeded but only the developer sees it.**
  `--debug` (default) deploys are per-developer. Fix: re-publish with
  `--release` and set `canView`/`canEdit`, or share via the UI.
  `--debug` and `--release` are mutually exclusive (`DG-FACT-156`).

## Verification

- `gh run list --workflow=packages.yml --limit 1` shows the latest
  publish run as `completed / success`.
- `curl -s https://registry.npmjs.org/<name>/<version> | jq -r .version`
  returns `<version>` (not `null`).
- For the private path: log in to `<HOST>`, **Manage → Packages**,
  filter by package name — the new version is listed and shareable.
- Slack `#build` channel has a `Package <Pkg> version <v> published to
  NPM` message linking to `npmjs.com/package/<name>/v/<version>`.

## See also

- Source articles:
  - `help/develop/how-to/packages/publish-packages.md` (primary).
  - `help/develop/how-to/tests/test-packages.md` (`#local-testing`,
    `--host`, `--skip-build`, `--skip-publish` flags).
  - `help/develop/develop.md` (`#publishing` — `--debug` vs
    `--release`, `canEdit`/`canView`, transitive-dep warning at
    `develop.md:284-291`).
  - `help/collaborate/public-repository.md` (linked from the article).
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-144` … `DG-FACT-156` and drifts
  `DG-FACT-DRIFT-056` … `DG-FACT-DRIFT-059`.
- Reference packages:
  - `packages/Charts/package.json` — canonical scoped layout, `1.7.0`,
    `release-charts` script form.
  - `packages/ApiTests/package.json` — intentionally pins
    `"datagrok-api": "../../js-api"`, so CI never publishes it to npm.
- Workflow: `.github/workflows/packages.yaml` — the matrix and gates
  this skill walks through (lines `84`, `132-156`, `613-635`).
- Related skills:
  - `manage-credentials` (POST credentials after publish).
  - `db-in-plugin` (migrations require `--release`).
