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

The build process uses Babel with `@babel/preset-typescript` to transpile TypeScript files from `bin/` to `bin/` (in-place). **Important:** `.ts` source and `.js` output coexist in the same `bin/` directory — `grok.js` requires the transpiled `.js` files, not the `.ts` sources. After editing any `.ts` file, you must run `npm run build` before testing.

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
- `build.ts` - Build one package or recursively build all packages in a directory
- `test.ts` - Run Puppeteer-based tests for a single package
- `test-all.ts` - Run tests across multiple packages
- `stress-tests.ts` - Run stress tests (must be run from ApiTests package)
- `api.ts` - Auto-generate TypeScript wrappers for scripts/queries
- `link.ts` - Link libraries for plugin development
- `claude.ts` - Launch a Dockerized dev environment with Datagrok + Claude Code
- `migrate.ts` - Update legacy packages
- `init.ts` - Apply configuration to existing packages

The commands `api`, `check`, `link`, `publish`, and `test` support the `--all` flag to run recursively across all packages in the current directory.

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

### `grok build` Command

Builds packages with `npm install` + `npm run build`. Supports:
- Single package: `grok build` (from package directory)
- Recursive: `grok build --recursive` (discovers and builds all packages in subdirectories)
- `--filter "name:Chem"` - Filter packages by package.json fields (supports regex, `&&` for multiple conditions)
- `--parallel N` - Max parallel build jobs (default 4)
- `--no-incremental` - Force full rebuild (default uses `--env incremental`)

### `grok claude` Command

Launches a full Dockerized development environment with Datagrok + Claude Code. The compose configuration is embedded in `claude.ts` — no external files needed.

```bash
grok claude <project>                      # Create worktree + start containers + launch Claude
grok claude <project> --in-place           # Use current directory (no worktree)
grok claude <project> --keep               # Leave containers running on exit
grok claude <project> --version 1.22.0     # Pin Datagrok version
grok claude <project> --profile full       # Include spawner, JKG, demo DBs
grok claude <project> --profile scripting  # Include JKG for Python/R/Julia
grok claude <project> --port 8080          # Fix Datagrok port (default: random free port)
grok claude <project> --prompt "fix bug"   # Pass prompt to Claude Code (-p flag)
grok claude destroy <project>              # Tear down containers + worktree + temp files
grok claude destroy-all                    # Destroy all known projects
```

**Project name restrictions:** `master` and `main` are rejected.

**Lifecycle:**
1. Creates git worktree at `~/pkg-worktrees/<project>` (unless `--in-place` or not in a git repo)
2. Writes `docker-compose.yaml` + `.env` + optional `docker-compose.override.yaml` to `$TMPDIR/dg-pkg-<project>`
3. Runs `docker compose up -d --wait`
4. Detects Claude working directory based on repo type (see below)
5. Launches `claude --dangerously-skip-permissions` inside the `tools-dev` container
6. On exit: stops containers (unless `--keep`)

**Version resolution:** `--version` flag > `bleeding-edge` (if inside public repo) > `latest`

**Compose services (embedded template):**
- Always started: `postgres` (pgvector:pg17), `rabbitmq`, `grok_pipe`, `datagrok`, `grok_connect`, `tools-dev`
- Profile `full`: adds `grok_spawner`
- Profile `scripting`/`full`: adds `jupyter_kernel_gateway`
- Profile `demo`/`full`: adds `world`, `test_db`, `northwind` demo databases

**Note:** The embedded template in `claude.ts` differs from `.devcontainer/docker-compose.yaml`:
- Embedded uses separate version vars per service (`DATAGROK_VERSION`, `GROK_CONNECT_VERSION`, `GROK_SPAWNER_VERSION`, `JKG_VERSION`, `TOOLS_DEV_VERSION`); `.devcontainer/` uses a single `DG_VERSION` for all
- Embedded mounts workspace at `/workspace/repo`; `.devcontainer/` mounts at `/workspace`
- Embedded has `grok_connect` always-on (no profile); `.devcontainer/` puts it under `profiles: ["full"]`
- Embedded uses `grok_pipe:latest`; `.devcontainer/` uses `grok_pipe:${DG_VERSION}`
- Embedded doesn't mount `~/.grok` from host; `.devcontainer/` does

**Host config override:** If `~/.claude` (or `CLAUDE_HOME`) is found on the host, a `docker-compose.override.yaml` is generated to bind-mount it into the container. Also mounts `~/.claude.json` if present.

**Working directory detection inside container:**
- Public repo root (has `js-api/`): `/workspace/repo`
- Monorepo (has `public/js-api/`): `/workspace/repo/public`
- External repo: `/workspace/datagrok/packages/<folder-name>` (waits up to 600s for entrypoint to clone public repo)

#### tools-dev Container (`Dockerfile.pkg_dev`)

Based on `node:22-bookworm-slim`. Pre-installed:
- Google Chrome stable (for Puppeteer), Playwright + Chromium
- `datagrok-tools` (grok CLI) and `@anthropic-ai/claude-code` (global npm)
- git, curl, jq, docker CLI
- Runs as `node` user (UID 1000, added to `docker` group)

#### Entrypoint (`entrypoint.sh`)

1. **Repo detection:** checks if `/workspace/repo` is the public repo (has `.git` + `js-api/`) or monorepo (`public/js-api/`)
2. **Auto-clone:** if workspace is not the public repo, sparse-clones it to `$DG_PUBLIC_DIR` (default `/workspace/datagrok`) — excludes `connectors/`, `docker/`, `environments/`, `python-api/`, etc. for speed. Branch resolved as: `DG_PUBLIC_BRANCH` > `DG_VERSION` mapped to branch > `master` fallback
3. **Workspace linking:** for non-public repos, symlinks `/workspace/repo` into the cloned repo's `packages/` dir and links `.claude`/`CLAUDE.md` at `/workspace/` for context discovery
4. **Grok config:** auto-creates `~/.grok/config.yaml` pointing to `http://datagrok:8080/api` with key `admin` (only if config doesn't already exist)

#### Profiles

| Profile | Additional services |
|---------|-------------------|
| (none) | postgres, rabbitmq, grok_pipe, datagrok, grok_connect, tools-dev |
| `scripting` | + jupyter_kernel_gateway |
| `demo` | + world, test_db, northwind |
| `full` | + grok_spawner, JKG, demo DBs |

See `.devcontainer/PACKAGES_DEV.md` for detailed usage docs, architecture diagram, MCP plugin setup (Jira/GitHub), and troubleshooting.

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
