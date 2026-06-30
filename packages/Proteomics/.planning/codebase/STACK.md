# Technology Stack

**Analysis Date:** 2026-05-11

## Languages and Runtime

| Language | Version | Use |
|----------|---------|-----|
| TypeScript | 5.6.x | Primary source language for the package |
| JavaScript (ES) | ES2022 target | Compiled output from TypeScript |
| Plain JS | ES module | `detectors.js` (loaded separately by platform — must not use ES module `import`) |
| R | 4.x via Bioconductor | Server-side scripts in `scripts/` |

Browser target: any platform-supported Chromium (Datagrok-supported browsers). Server-side R runs inside the Datagrok platform, not on the developer's machine.

## Frameworks and Platform

| Framework | Version | Role |
|-----------|---------|------|
| `datagrok-api` | ^1.25.0 | Core Datagrok JS/TS API (namespaces: `grok`, `ui`, `DG`) |
| `@datagrok-libraries/statistics` | ^1.2.12 | `tTest` (Welch's t-test) and `fdrcorrection` (Benjamini-Hochberg FDR) |
| `@datagrok-libraries/bio` | ^5.61.3 | Tree-helper APIs; currently imported but the active code paths use `DG.Func.find` for cross-package lookup |
| `@datagrok-libraries/test` | ^1.1.0 | Test framework — `category()`, `test()`, `expect()` |
| `rxjs` | ^6.5.5 | Subscription pattern for cross-DataFrame linkage in `enrichment-viewers.ts` |

## Build Pipeline

Standard Datagrok package build, defined in `package.json` `scripts.build`:

```
grok api && grok check --soft && webpack
```

1. **`grok api`** — Reads `@grok.decorators.*` metadata from `src/package.ts` and regenerates `src/package.g.ts` and `src/package-api.ts`. These files are auto-generated and must not be edited.
2. **`grok check --soft`** — Validates function signatures, input/output declarations, imports, and `package.json` shape. `--soft` downgrades errors to warnings so the build continues.
3. **`webpack`** — Bundles to `dist/package.js` (main bundle) and `dist/package-test.js` (test bundle). Webpack config in `webpack.config.js` defines two entry points.

Other relevant scripts in `package.json`:

| Script | Purpose |
|--------|---------|
| `lint` / `lint-fix` | ESLint over `src/**/*.ts` |
| `test` / `test-local` / `test-dev` | `grok test` against various servers |
| `link-all` | `npm link` to local `datagrok-api` and `@datagrok-libraries/*` |
| `build-all` | Build js-api + libraries + this package in dependency order |
| `debug-proteomics` / `release-proteomics` | `grok publish` aliases |

## TypeScript Configuration

`tsconfig.json`:
- `target`: ES2020
- `strict`: true
- `experimentalDecorators`: true, `emitDecoratorMetadata`: true (for `@grok.decorators.*`)
- Module: ES modules
- Output: `dist/` (cleaned before each build)

## Webpack Externals

Defined in `webpack.config.js`. These are **provided by the platform at runtime** and must NOT be bundled:

```js
externals: {
  'datagrok-api/dg':    'DG',
  'datagrok-api/grok':  'grok',
  'datagrok-api/ui':    'ui',
  'rxjs':               'rxjs',
  'rxjs/operators':     'rxjs.operators'
}
```

Importing any of these in TypeScript pulls the global symbol at runtime — bundling them would conflict with the platform's own copies.

## R Server-Side Environment

Each R script in `scripts/` declares its conda/Bioconductor environment inline via the `#environment:` metadata header. Example from `scripts/limma_de.R` (added in commit `489bb53ec9`):

```r
#environment: channels: [conda-forge, bioconda], dependencies: [bioconductor-limma]
```

This tells the Datagrok platform to pre-warm an R session with the named Bioconductor packages installed. The package currently declares:

| Script | Environment |
|--------|-------------|
| `scripts/limma_de.R` | `bioconductor-limma` |
| `scripts/deqms_de.R` | `bioconductor-deqms` (transitively pulls limma) |
| `scripts/vsn_normalize.R` | `bioconductor-vsn` |

R scripts register as Datagrok functions via the `#name:` / `#language: r` headers. They are invoked from TypeScript via `grok.functions.call('Proteomics:LimmaDE', {...})`.

## Configuration Files

| File | Purpose |
|------|---------|
| `package.json` | npm package metadata, dependencies, scripts |
| `tsconfig.json` | TypeScript compiler config |
| `webpack.config.js` | Two entry points (`src/package.ts`, `src/package-test.ts`), externals, output to `dist/` |
| `.eslintrc.json` | ESLint rules (Google style, 120-char lines, 2-space indent) |
| `.npmignore` | Files excluded from `grok publish` |
| `package-lock.json` | Pinned npm dependency tree |

## Package Identity

- **npm name:** `@datagrok/proteomics`
- **Version:** `0.1.0` (in local source; published debug version is `admin`)
- **Friendly name in platform:** `Proteomics`
- **Category:** `Bioinformatics`

---

*Stack analysis: 2026-05-11*
