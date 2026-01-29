# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is **datagrok-tools**, a CLI utility for creating, validating, testing, and publishing packages to Datagrok. The tool is distributed as `grok` command globally via npm.

## Build & Development Commands

### Build the CLI
```bash
npm run build                    # Transpile TypeScript to JavaScript using Babel
npm run debug-source-map         # Build with source maps for debugging
```

The build process uses Babel with `@babel/preset-typescript` to transpile TypeScript files from `bin/` (TypeScript source) to `bin/` (JavaScript output). The source TypeScript files are transpiled in-place.

### Link for Local Development
```bash
npm link                         # Make 'grok' command available globally for testing
```

### Testing
While this CLI tool manages testing for Datagrok packages, it doesn't have its own test suite. To test changes:
1. Build the tool: `npm run build`
2. Link it: `npm link`
3. Test commands manually: `grok create test-package`, `grok check`, etc.

## Code Architecture

### Command Structure (`bin/commands/`)

The CLI uses a modular command pattern. Each command is a separate module that:
- Exports a default function taking `argv` (parsed arguments)
- Returns `boolean` (true for success, false for failure)
- Is imported and mapped in `bin/grok.js`

**Available commands:**
- `config.ts` - Manage server configuration (`~/.grok/config.yaml`)
- `create.ts` - Create new packages from templates
- `add.ts` - Add entities (functions, scripts, queries, etc.) to packages
- `publish.ts` - Upload and deploy packages to Datagrok servers
- `check.ts` - Validate package structure, signatures, imports
- `test.ts` - Run Puppeteer-based tests for a single package
- `test-all.ts` - Run tests across multiple packages
- `api.ts` - Auto-generate TypeScript wrappers for scripts/queries
- `link.ts` - Link libraries for plugin development
- `migrate.ts` - Update legacy packages
- `init.ts` - Apply configuration to existing packages

### Template System

**Package Templates** (`package-template/`)
- Supports TypeScript (default) and JavaScript packages
- Uses placeholder replacement: `#{PACKAGE_NAME}`, `#{PACKAGE_NAME_LOWERCASE}`, `#{CURRENT_DATE}`, etc.
- Includes webpack config, tsconfig, package.json, and optional IDE/ESLint configs

**Entity Templates** (`entity-template/`)
Supports 8 entity types:
- Code: `function`, `app`, `script`, `detector`, `view`, `viewer`
- Data: `connection`, `query`
- Test template with `package-test.ts`

**Script Templates** (`script-template/`)
Language-specific templates for: JavaScript, Node.js, Python, R, Julia, Octave

### Utilities (`bin/utils/`)

**utils.ts** - Core utilities:
- Case converters: `kebabToCamelCase()`, `camelCaseToKebab()`, `spaceToCamelCase()`
- Validators: `isPackageDir()`, `isValidCron()`, `isEmpty()`
- Template placeholder replacement system
- Constants: `scriptLangExtMap`, `cacheValues`, `propertyTypes`

**test-utils.ts** - Testing infrastructure:
- Puppeteer browser automation: `getBrowserPage()`
- Authentication: `getToken()`, `getDevKey()`
- Test loading and CSV export utilities

**color-utils.ts** - Console styling: `error()`, `warn()`, `success()`, `info()`

**func-generation.ts** - Auto-generates TypeScript wrappers with type signatures

**ent-helpers.ts** - User-friendly help messages for entity types

### Validation System (`bin/validators/`)

**check.ts** (~695 lines) performs comprehensive validation:
- Function signatures match implementation
- Imports match webpack externals (prevent bundle bloat)
- Package.json structure and fields
- CHANGELOG.md format and version matching
- Script naming conventions (alphanumeric only)
- Parameter validation (types, names, optional/nullable)
- Cache and cron schedule validation

**config-validator.ts** validates `~/.grok/config.yaml`:
- Checks for required fields (default server, servers object, URLs, keys)
- Returns structured result: `{value: boolean, message: string, warnings: string[]}`

### Publishing Workflow

The `publish` command executes these steps:
1. Gather files using `ignore-walk` (respects .npmignore/.gitignore/.grokignore)
2. Run validation checks (signatures, imports, package.json, changelog)
3. Process environment variables in `/connections/*.json` files (replace `${VAR}`)
4. Create ZIP archive with archiver-promise
5. Upload to server: `POST ${host}/packages/dev/${devKey}/${packageName}`

**Key flags:**
- `--debug` (default) - Package visible only to developer
- `--release` - Public package
- `--build` / `--rebuild` - Control webpack bundling
- `--skip-check` - Skip validation

### Testing Framework

Tests use Puppeteer for headless browser automation:

**Test execution flow:**
1. Build package (skip with `--skip-build`)
2. Publish package (skip with `--skip-publish`)
3. Load test list from `package-test.ts`
4. Launch Puppeteer browser
5. Execute tests with optional recording (`--record`) or GUI mode (`--gui`)
6. Report results (console, CSV with `--csv`, or submit with `--report`)

**Important flags:**
- `--gui` - Visual debugging (headless=false)
- `--verbose` - Detailed test output
- `--catchUnhandled` - Catch unhandled exceptions
- `--debug` - Debug breakpoints (requires `--gui`)

## Key Patterns and Conventions

### Naming Conventions

**Package names:** Letters, numbers, underscores, hyphens only (`^([A-Za-z\-_\d])+$`)
**Entity names:** Start with letter, alphanumeric thereafter (`^([A-Za-z])+([A-Za-z\d])*$`)
**Developer keys:** Alphanumeric + hyphens (`^([A-Za-z\d-])+$`)

### Versioning

Semantic versioning enforced: `X.Y.Z` or `X.Y.Z-rc` or `X.Y.Z-rc.N`

CHANGELOG.md format:
```markdown
## X.Y.Z (YYYY-MM-DD | WIP)
```

Latest version in CHANGELOG must match package.json version (enforced for packages v1.0+).

### Function Metadata Pattern

Functions use JSDoc-style metadata comments:
```javascript
//name: functionName
//description: What it does
//input: type paramName {optional: true, nullable: true}
//output: type returnName
//meta.role: panel | viewer | app
//meta.cache: global | session | client | none
//meta.invalidateOn: 0 */5 * * * * (cron schedule)
```

### Environment Variable Substitution

Connection files (`/connections/*.json`) support `${VAR_NAME}` placeholders, replaced at publish time with environment variables.

### Configuration File

Location: `~/.grok/config.yaml`

Structure:
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

### Webpack Externals

The CLI validates that imports match webpack externals to prevent bundling datagrok-api and other provided libraries. Common externals:
- `datagrok-api`
- `rxjs`
- `cash-dom`
- `dayjs`
- `openchemlib/full`
- `wu`

## File References

When working with commands, key entry points are:
- `bin/grok.js` - Main CLI router
- `bin/commands/<command>.ts` - Individual command implementations
- `bin/utils/utils.ts` - Core utility functions
- `bin/validators/check.ts` - Package validation logic

When modifying templates:
- `package-template/` - Full package scaffolding
- `entity-template/` - Individual entity templates
- `script-template/` - Language-specific script templates
