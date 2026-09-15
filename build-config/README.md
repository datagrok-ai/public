# @datagrok/build-config

The shared toolchain for Datagrok plugins and libraries: the rspack + swc bundler configuration, the
TypeScript base configs, the eslint rule set, and the TypeScript compiler version. It has no
commands of its own: `grok build` and `grok tsc` (from `datagrok-tools`) use it. Inside the `public/`
workspace it is a root devDependency and packages declare nothing; a standalone package created with
`grok create` gets it as a devDependency.

## Package scripts

```json
"scripts": {
  "build": "grok build",
  "typecheck": "grok tsc --noEmit -p tsconfig.json",
  "lint": "eslint --ext .ts,.tsx src",
  "test": "grok test"
}
```

`grok build` in a package (or inside a Turborepo task) bundles `src/package.ts` into `dist/package.js`,
generates `src/package.g.ts` and `src/package-api.ts`, and runs `grok check --soft`; in a library it
runs `tsc -p tsconfig.json` into `dist/` and mirrors the css/wasm/json assets the sources import.
`grok tsc ...` is the workspace TypeScript compiler with the given arguments.

## tsconfig

```json
{"extends": "@datagrok/build-config/tsconfig.base.json", "include": ["src"]}
```

`tsconfig.base.json`: ES2022 target, `bundler` resolution, strict, legacy decorators with metadata,
`types: ["node", "@datagrok/build-config/css"]` (ambient declarations for `.css`, images and `.wasm`
imports), no emit. Libraries extend `tsconfig.lib.json` instead: emits to `dist/` with declarations
and source maps, ES2020.

## eslint

`eslintrc.json` is the rule set (the Google style the repository always used, TypeScript parser). The
workspace root config extends it; a standalone package's `.eslintrc.json` extends
`./node_modules/@datagrok/build-config/eslintrc.json`.

## rspack.config.js (only when a package deviates from the defaults)

```js
const {bundler} = require('@datagrok/build-config');
module.exports = bundler({
  externals: {ngl: 'NGL'},        // merged with the platform externals; `false` removes one
  wasm: 'async',                  // 'async' | 'sync' | 'asset' (default: .wasm as a URL asset)
  jsx: 'react',                   // swc React automatic runtime for .tsx
  css: false,                     // drop the default style-loader rule to supply your own
  rules: [...], plugins: [...],   // prepended / appended
  resolve: {alias: {...}, fallback: {...}},
  mode: 'development',
});
```

A config may also export a function; `grok build --key=value` arguments arrive as its `env` object
(js-api uses `--only=browser|node`).

Defaults: entries `src/package.ts` and `src/package-test.ts` (or `.js`), output `dist/package.js` and
`dist/package-test.js` as `var` libraries named after the package, source maps, production mode with
swc minification, platform externals (`datagrok-api/*`, rxjs, cash-dom, dayjs, wu, openchemlib,
exceljs, html2canvas), assets as URLs, node core fallbacks off, library `dist/` downlevelled to
Chrome 50 like the package sources. `rspack` and `loaders` (absolute paths of style/css/null loaders)
are exported for hand-written configs.
