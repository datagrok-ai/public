---
title: "Build system"
sidebar_position: 1
description: How to install, build, and publish packages, libraries, and the JS API in the Datagrok public repository, which is one pnpm workspace built with Turborepo.
keywords:
  - pnpm
  - workspace
  - Turborepo
  - grok build
  - grok setup
  - build packages
  - libraries
  - js-api
---

The [public repository](https://github.com/datagrok-ai/public) is one [pnpm](https://pnpm.io) workspace.
You install dependencies once for the whole repository, and one command builds a package together with
the libraries it depends on. This page covers the daily commands, what a package declares in
`package.json`, and what changed from the earlier npm setup.

:::note

This page applies to code inside the public repository. A standalone package that you create with
[`grok create`](../how-to/packages/create-package.md) outside the repository uses `npm install` and
`npm run build`, and needs nothing from this page.

:::

## Set up

To prepare a checkout, run `grok setup` at the repository root. You need
[Node.js](https://nodejs.org) 20 or later and
[datagrok-tools](set-up-environment.md).

```shell
grok setup
```

`grok setup` does the following:

* Enables pnpm through `corepack` at the version the repository pins
* Removes `node_modules` folders and `package-lock.json` files that npm left inside packages, and the
  `.js`/`.d.ts` files that the npm-era `tsc` emitted next to JS API and library sources
* Runs `pnpm install` for the whole workspace
* Warns you when the global `grok` is older than the one in the repository

Run `grok setup` again after you pull changes that touch `pnpm-lock.yaml`. To see what it would do
without changing anything, run `grok setup --check`. If you prefer plain pnpm, `corepack enable` followed by
`pnpm install` at the root gives the same result.

## Build

To build a package, run `grok build` in its folder. The build also builds the JS API and every library
the package depends on, in the right order. [Turborepo](https://turbo.build) caches the results, so the
next build skips whatever didn't change.

| Command                   | Builds                                                   |
|:--------------------------|:---------------------------------------------------------|
| `grok build`              | The package in the current folder and what it depends on |
| `grok build --all`        | Every package and library                                |
| `grok build --affected`   | Everything your changes against `origin/master` touch    |
| `grok build --typecheck`  | The same as `grok build`, plus a type check              |
| `grok build --local`      | Only the current folder, without its dependencies        |
| `grok publish <server>`   | A debug build, deployed to a [server](set-up-environment.md) |

`grok build` runs Turborepo for you. You can also call it directly:

```shell
pnpm turbo run build --filter=@datagrok/chem...
```

### Package scripts

Every package has the same four scripts, so `npm run <script>` and `pnpm run <script>` behave the same
way everywhere:

| Script      | Does                                                                        |
|:------------|:----------------------------------------------------------------------------|
| `build`     | Bundles the package to `dist/`, generates function metadata, and runs `grok check --soft` |
| `typecheck` | Type-checks with the TypeScript version the workspace provides              |
| `lint`      | Runs ESLint with the repository config                                      |
| `test`      | Runs `grok test` against the [configured server](../how-to/tests/test-packages.md) |

Libraries and the JS API compile to `dist/` with type declarations. Packages consume that output and
never compile library sources.

## Declare dependencies

A package declares only what it imports at runtime. The workspace root provides the toolchain (bundler,
TypeScript, ESLint, and `grok`), so a package needs no build-related `devDependencies`.

```json
"dependencies": {
  "datagrok-api": "workspace:^",
  "@datagrok-libraries/utils": "workspace:^",
  "rxjs": "catalog:",
  "ngl": "^2.0.0"
}
```

Use these version specifiers:

* `workspace:^` for `datagrok-api` and `@datagrok-libraries/*`. The package uses the copy in the
  repository. When the package is published to npm, pnpm replaces the specifier with the real version.
* `catalog:` for libraries that have one shared version in the repository, listed in
  `pnpm-workspace.yaml`
* A regular version range for everything else

To add a dependency, run `pnpm add <name>` in the package folder. Don't use `npm install` inside the
repository. It creates a `package-lock.json` and a separate `node_modules` that conflict with the
workspace.

### Platform libraries

Datagrok serves some libraries at runtime: `cash-dom`, `codemirror`, `dayjs`, `exceljs`, `html2canvas`,
`ngl`, `openchemlib`, `rxjs`, `vue`, and `wu`. The build treats them as externals, and their versions live
in `build-config/platform-deps.json`. Declare them as `catalog:` so that your package compiles against the
version the platform ships. If a package pins its own version instead, `grok check` warns about it.

### Test dependencies

If the tests of your package need another package on the server, list it in `grok.testDependencies`:

```json
"grok": {
  "testDependencies": ["@datagrok/chem", "@datagrok/power-grid"]
}
```

Continuous integration (CI) builds these packages and publishes them to the test server before it runs
your tests. Don't add other packages to `devDependencies` for this purpose.

### Bundler configuration

Most packages need no bundler configuration, and their `tsconfig.json` only extends the shared base. If
your package deviates from the defaults, add an `rspack.config.js`:

```js
module.exports = require('@datagrok/build-config').bundler({externals: {ngl: 'NGL'}, wasm: 'async', jsx: 'react'});
```

## Change a library or the JS API

To try a change in a library or the JS API, edit it and run `grok build` in the package that uses it.
The build compiles the library first because the package depends on it. You don't link anything, and you
don't bump versions for local work. In a workspace, `grok link` does nothing.

## Publish to npm

To release a package or library to npm, bump its version and merge to `master`. CI publishes every
changed package whose version isn't on npm yet. It doesn't publish `0.x` versions. `grok publish` deploys to a
Datagrok server and is unrelated to npm. For details, see [CI flow](ci-flow.mdx) and the
[versioning policy](versioning-policy.md).

## Migrate from npm

If you worked in the repository before it became a workspace, replace your habits as follows:

| Before                                         | Now                                          |
|:-----------------------------------------------|:---------------------------------------------|
| `npm install` in each package                  | `grok setup` once at the root                |
| `grok link` or `npm link`                      | Nothing. `workspace:^` links packages        |
| `npm run build` in js-api, then each library, then the package | `grok build` in the package   |
| `npm install <name>`                           | `pnpm add <name>` in the package folder      |
| `package-lock.json` in each package            | One `pnpm-lock.yaml` at the root             |
| `webpack.config.js` and build `devDependencies` in each package | Shared `@datagrok/build-config` |
| `@datagrok/*` in `devDependencies` for tests   | `grok.testDependencies`                      |

## Troubleshooting

| Problem                                                        | Solution                                           |
|:---------------------------------------------------------------|:---------------------------------------------------|
| `pnpm: command not found`                                      | Run `corepack enable`, or run `grok setup`         |
| `ERR_PNPM_OUTDATED_LOCKFILE` in CI                             | Run `pnpm install` at the root and commit `pnpm-lock.yaml` |
| The build can't resolve a module after you switch branches     | Run `grok setup` again                             |
| A package has its own `package-lock.json` or a large `node_modules` | Run `grok setup`. It removes them             |
| A library change doesn't show up in the package                | Run `grok build` in the package, not in the library. To ignore the cache, add `--force` |

See also:

* [Environment setup](set-up-environment.md)
* [Packages](../develop.md)
* [CI flow](ci-flow.mdx)
