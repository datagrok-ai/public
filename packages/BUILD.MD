# Building plugins, libraries and js-api

`public/` is one pnpm workspace with a Turborepo task graph. There is nothing to link and nothing
to install per package.

```bash
grok setup                        # once per checkout or worktree: pnpm via corepack, install, npm-era clean-up (plain `pnpm install` works too)
grok build                        # this package and what it depends on, in order, cached
grok build --all                  # everything
grok build --affected             # everything my diff against origin/master touches (CI does this)
grok build --typecheck            # bundle and type-check
pnpm turbo run build --filter=@datagrok/chem...   # the same thing, spelled in Turborepo
grok publish dev                  # deploy a debug build to a server alias
```

Inside a package the four scripts are always the same:

| script | does |
|---|---|
| `build` | `grok build`: rspack + swc bundle to `dist/`, function metadata generation, `grok check --soft` |
| `typecheck` | `grok tsc --noEmit -p tsconfig.json`: `tsc --noEmit` with the workspace TypeScript |
| `lint` | `eslint --ext .ts,.tsx src` with the repository config |
| `test` | `grok test` against the configured server |

Libraries and js-api use `grok tsc -p tsconfig.json`, emit to `dist/` with declarations, and expose
their files through `exports` in package.json; consumers never compile library sources.

## What a package declares

Only what it imports at runtime. The toolchain (rspack, swc, TypeScript, eslint, `dg`, `grok`) is
provided by the workspace root through `@datagrok/build-config`, and versions of shared libraries
come from the catalog in `pnpm-workspace.yaml`:

```json
"dependencies": {
  "datagrok-api": "workspace:^",
  "@datagrok-libraries/utils": "workspace:^",
  "rxjs": "catalog:",
  "ngl": "^2.0.0"
}
```

`tsconfig.json` is two lines (`extends` the shared base, `include: ["src"]`). A package needs an
`rspack.config.js` only when it deviates from the defaults:

```js
module.exports = require('@datagrok/build-config').bundler({externals: {ngl: 'NGL'}, wasm: 'async', jsx: 'react'});
```

## Publishing to npm

Bump the version and merge to master; the Workspace workflow publishes every changed package whose
version is not on npm yet (0.x never publishes). `grok publish` deploys to a Datagrok server; it is
unrelated to npm.

## Changing js-api or a library

Edit it, run `grok build` in the plugin that uses it. Turborepo rebuilds the library first because the
plugin depends on it. Nothing to link, nothing to bump for local work.
