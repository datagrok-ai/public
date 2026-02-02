# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is the **Datagrok public repository** - the open-source monorepo for the Datagrok data analytics and visualization platform. It contains the JavaScript/TypeScript API, shared libraries, CLI tools, and 76+ extension packages.

## Repository Structure

```
public/
├── js-api/          # Core JavaScript/TypeScript API (datagrok-api)
├── libraries/       # Shared TypeScript libraries (@datagrok-libraries/*)
├── packages/        # Extension packages (viewers, connectors, scientific tools)
├── tools/           # CLI tool "grok" (datagrok-tools)
├── python-api/      # Python API bindings
├── help/            # Documentation
├── connectors/      # External service connectors
├── docker/          # Docker deployment configs
├── environments/    # Environment configurations
└── misc/            # ESLint config, utilities
```

## js-api (datagrok-api)

The core TypeScript API providing bindings to the Datagrok Dart backend. Packages import three namespaces:

- **`grok`** - High-level APIs: shell, functions, events, data operations, server API, AI
- **`ui`** - UI components: elements, dialogs, inputs, menus, viewers
- **`dg`** - Re-exports all types, constants, and classes

Key modules in `js-api/src/`: `dataframe.ts` (DataFrame, Column, BitSet), `entities.ts` (User, Project, Package), `widgets.ts` (Dialog, Menu, InputBase), `viewer.ts`, `grid.ts`, `dapi.ts` (server HTTP API), `shell.ts`, `functions.ts`, `events.ts`, `const.ts` (enums).

Files ending in `.g.ts` or `.api.g.ts` are auto-generated from the Dart codebase - do not edit manually.

```bash
cd js-api
npm run build           # tsc && webpack
npm run build-ts        # TypeScript + ESLint fix only
npm run build-js-api    # Clean, compile, bundle
```

## Libraries (@datagrok-libraries/*)

Reusable TypeScript libraries published under `@datagrok-libraries` scope in `libraries/`:

| Library | Purpose |
|---------|---------|
| **bio** | Macromolecule & Molecule3D data support |
| **chem-meta** | RDKit JS API, molfile parsing |
| **compute-api** | Compute global API |
| **compute-utils** | Compute-related utilities |
| **cruddy** | CRUD app framework |
| **db-explorer** | Database exploration |
| **dock-spawn-dg** | TypeScript docking manager |
| **gridext** | Grid extensions |
| **math** | High-performance computation |
| **ml** | Machine learning utilities |
| **statistics** | Statistics & aggregation |
| **test** | Common test utilities (Puppeteer) |
| **tutorials** | Tutorial creation helpers |
| **utils** | Common utilities (hashing, encoding) |
| **webcomponents** | Web components wrappers |
| **webcomponents-vue** | Vue.js web components wrappers |

Each library has its own `package.json` and builds independently with `npm run build`.

## Tools (datagrok-tools / `grok` CLI)

The `grok` CLI manages the full package lifecycle. Install globally via `npm install -g datagrok-tools` or link locally from `tools/`.

### Key Commands

```bash
grok create <name>       # Create new package from template
grok add <type>          # Add function/script/query/app/view/viewer to a package
grok check               # Validate package (signatures, imports, package.json, changelog)
grok check --soft        # Validate with warnings only (non-blocking)
grok api                 # Auto-generate TypeScript wrappers from function metadata
grok publish             # Upload package to default server (debug mode)
grok publish --release   # Upload as public release
grok publish --build     # Build webpack before publishing
grok link                # Link all local libraries and js-api for development
grok link --path         # Update package.json dependencies to local paths instead of npm link
grok link --unlink       # Revert to versioned dependencies
grok test                # Run Puppeteer-based tests
grok test --gui          # Run tests with visible browser
grok test --verbose      # Detailed test output
grok config              # Manage server configuration (~/.grok/config.yaml)
```

### Configuration

Server config is stored in `~/.grok/config.yaml`:

```yaml
default: 'public'
servers:
  public:
    url: 'https://public.datagrok.ai/api'
    key: 'developer-key-here'
  dev:
    url: 'https://dev.datagrok.ai/api'
    key: ''
```

## Binding Packages to Local Libraries and js-api

When developing packages, you need local copies of `datagrok-api` and `@datagrok-libraries/*` instead of npm versions.

### Method 1: `grok link` (Recommended)

Run from the package directory. Automatically discovers and links all dependencies found in the repository:

```bash
cd packages/MyPackage
grok link                # npm link all local dependencies
grok link --verbose      # Show detailed linking progress
grok link --dev          # Also link devDependencies
grok link --repo-only    # Only link packages from this repository
grok link --unlink       # Revert to npm versioned dependencies
```

This recursively resolves dependency chains (e.g., if your package depends on `@datagrok-libraries/bio`, which depends on `@datagrok-libraries/utils`, both get linked).

### Method 2: `grok link --path`

Instead of npm symlinks, rewrites `package.json` dependencies to relative paths:

```bash
grok link --path         # Changes "datagrok-api": "^1.26.0" → "datagrok-api": "../../js-api"
grok link --unlink       # Reverts paths back to version numbers
```

### Method 3: Manual npm link

First register packages globally, then link them in your package:

```bash
# Register js-api and libraries globally
cd js-api && npm link
cd libraries/utils && npm link
cd libraries/ml && npm link

# Link into your package
cd packages/MyPackage
npm link datagrok-api @datagrok-libraries/utils @datagrok-libraries/ml
```

### Method 4: npm link-all script

Most packages define a `link-all` script in their `package.json`:

```bash
npm run link-all    # Runs: npm link datagrok-api @datagrok-libraries/utils ...
npm run link-api    # Runs: npm link datagrok-api (just the API)
```

## Building Packages

### Standard Build

```bash
cd packages/MyPackage
npm run build            # Typically: grok api && grok check --soft && webpack
```

The build pipeline:
1. `grok api` - Auto-generates TypeScript wrappers from function metadata comments
2. `grok check --soft` - Validates signatures, imports, package.json
3. `webpack` - Bundles the package for browser

### Build with Local Dependencies

To build a package along with its local dependencies:

```bash
npm run build-all        # Builds js-api → libraries → package (in dependency order)
```

Example from Chem package:
```bash
# build-all: npm --prefix ../../libraries/chem-meta run build &&
#            npm --prefix ../../js-api run build &&
#            npm --prefix ../../libraries/utils run build &&
#            npm --prefix ../../libraries/ml run build &&
#            npm run build
```

## Publishing Packages

```bash
# Debug mode (visible only to developer)
npm run debug-local      # webpack && grok publish local
grok publish             # Publish to default server
grok publish dev         # Publish to 'dev' server

# Release mode (public)
npm run release-local    # grok publish local --release
grok publish --release   # Publish release to default server

# Options
grok publish --build     # Run webpack before publishing
grok publish --rebuild   # Force rebuild
grok publish --skip-check # Skip validation
```

Publishing steps: gather files (respects .npmignore/.gitignore/.grokignore) → validate → process env vars in connection files → create ZIP → upload to server.

## Testing Packages

Tests run in a Puppeteer-controlled browser against a running Datagrok instance:

```bash
grok test                          # Run tests against default server
grok test --host localhost         # Test against local instance
grok test --host dev               # Test against dev server
grok test --gui                    # Visual browser (not headless)
grok test --verbose                # Detailed output
grok test --gui --debug            # Debug breakpoints
grok test --record                 # Record test execution
grok test --csv                    # Export results to CSV
grok test --skip-build             # Skip building before test
grok test --skip-publish           # Skip publishing before test
grok test --category "CAT"         # Run tests only from categories starting with CAT
grok test --test "TEST"             # Run only test starting with TEST
```

Most packages define convenience scripts:

```bash
npm run test             # grok test --host localhost
npm run test-dev         # grok test --host dev
```

## Function Metadata Pattern

Packages declare functions using JSDoc-style comments that `grok api` and `grok check` validate:

```typescript
//name: myFunction
//description: What it does
//input: dataframe df {caption: Input}
//input: column col {type: numerical}
//input: string name {optional: true}
//output: dataframe result
//meta.role: panel
//meta.cache: client
export function myFunction(df: DG.DataFrame, col: DG.Column, name?: string): DG.DataFrame {
  // ...
}
```

## Webpack Externals

Packages must not bundle these (they are provided by the platform at runtime):

- `datagrok-api` (including subpaths `datagrok-api/grok`, `datagrok-api/ui`, `datagrok-api/dg`)
- `rxjs`, `rxjs/operators`
- `cash-dom`
- `dayjs`
- `openchemlib/full`
- `wu`

The `grok check` command validates that imports match webpack externals.

## Code Style

- 2-space indentation
- Single quotes for strings
- Semicolons required
- Windows line endings (CRLF)
- TypeScript strict mode
- File naming: kebab-case (`my-component.ts`)
- Package naming: letters, numbers, underscores, hyphens
- Semantic versioning: `X.Y.Z` or `X.Y.Z-rc` or `X.Y.Z-rc.N`
