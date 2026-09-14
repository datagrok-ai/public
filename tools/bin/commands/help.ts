import { migrate } from "./migrate";
import { HELP_SERVER } from "./server";
import { testAll } from "./test-all";
import { HOME_ROOTS } from "../utils/kg/homes";

const HELP = `
Usage: grok <command>

Datagrok's package management tool

Commands:
    add         Add an object template
    api         Create wrapper functions
    build       Build a package or multiple packages
    check       Check package content (function signatures, etgc.)
    claude      Launch Claude Code in a Datagrok dev container
    config      Create and manage config files
    create      Create a package
    docker-gen  Generate Celery Docker artifacts from Python functions
    init        Modify a package template
    link        Link \`datagrok-api\` and libraries for local development
    publish     Upload a package
    report      Manage user error reports (fetch, resolve, create ticket)
    run         Build, publish, and open in browser
    test        Run package tests
    testall     Run packages tests
    migrate     Migrate legacy tags to meta.role
    server (s)  Manage a Datagrok server (list/get/delete entities, run functions)
    kg          Check and generate the knowledge graph (type files, home documents)

To get help on a particular command, use:
    grok <command> --help

Read more about the package development workflow:
https://datagrok.ai/help/develop/develop
`;

const HELP_CLAUDE = `
Usage: grok claude <project>              Start or reattach to a project
       grok claude destroy <project>      Stop containers + remove worktree
       grok claude destroy-all            Destroy all projects

Launch Claude Code inside a Datagrok dev container. Creates a git worktree
for the project, starts a full Datagrok stack (postgres, rabbitmq, grok_pipe,
datagrok) and opens Claude Code in a tools-dev container.

Version is auto-detected: bleeding-edge for the public repo,
latest stable release (from Docker Hub) for other repos.

Options:
[--version <tag>] [--profile <name>] [--keep]
[--port <N>] [--prompt <text>] [--in-place]
[--grok-connect-version <tag>] [--grok-spawner-version <tag>]
[--jkg-version <tag>] [--tools-dev-version <tag>]

--version               Datagrok image version (default: bleeding-edge for public repo, latest otherwise)
--profile               Compose profile: demo, scripting, full (default: none)
--keep                  Don't stop containers on exit
--port                  Datagrok host port (default: random available)
--prompt                Pass initial prompt to Claude Code (non-interactive)
--in-place              Use current directory instead of creating a git worktree
--grok-connect-version  grok_connect image version (default: latest)
--grok-spawner-version  grok_spawner image version (default: latest)
--jkg-version           jupyter_kernel_gateway image version (default: latest)
--tools-dev-version     tools-dev image version (default: latest)

Examples:
  grok claude GROK-12345                        Start working on a task
  grok claude GROK-12345 --version 1.22.0       Use specific Datagrok version
  grok claude GROK-12345 --profile full --keep   Start all services, keep running
  grok claude GROK-12345 --prompt "fix the bug"  One-shot command
  grok claude GROK-12345 --in-place              Work in current directory
  grok claude destroy GROK-12345                 Tear down a task
  grok claude destroy-all                        Tear down everything
`;

const HELP_ADD = `
Usage: grok add <entity> <name>

Add an object template to your package:

grok add app <name>
grok add app [name] --domain <schema>[.<table>]
grok add app [name] --domain <path to schema.json>
grok add connection <name>
grok add detector <semantic-type-name>
grok add function [tag] <name>
grok add query <name>
grok add script [tag] <language> <name>
grok add view <name>
grok add viewer <name>
grok add tests

Please note that entity names may only include letters and numbers

--domain scaffolds a working browse/CRUD app over an entity-mapped domain table
(\`grok.dapi.domains\`) from the \`@datagrok-libraries/domain-ui\` defaults. Give it
one table (\`--domain grit.issue\`), a whole schema the package declares in
\`databases/<schema>/schema.json\` (one app per table), or the path to a schema
manifest to copy into the package. A fresh app package is two commands:

grok create MyTracker
cd MyTracker && grok add app --domain grit.issue

Supported languages for scripts:
javascript, julia, node, octave, python, r

Available tags:
panel, init
`;

const HELP_INIT = `
Usage: grok init

Modify a package template by adding config files for linters, IDE, etc.

Options:
[--eslint] [--ide] [--test] [--ts] [--git]

--eslint    Add a configuration for eslint
--ide       Add an IDE-specific configuration for debugging (vscode)
--test      Add tests support (TypeScript packages only)
--ts        Convert a JavaScript package to TypeScript
--git       Configure GIT and install commit linting tools.
            Read more: https://datagrok.ai/help/develop/advanced/git-policy
`;

const HELP_API = `
Usage: grok api

Create wrapper functions for package scripts and queries.
Packages with \`databases/<schema>/schema.json\` manifests also get typed domain
clients in \`src/generated/db.ts\`.

Options:
[-v | --verbose] [--ui]

--verbose         Print detailed output
--ui              Also generate \`src/generated/db-ui.ts\` — typed UI wrappers over
                  \`@datagrok-libraries/domain-ui\` for every domain table. Once the
                  file exists, plain \`grok api\` keeps it up to date; delete it to
                  opt out again
`;

const HELP_CONFIG = `
Usage: grok config

Create or update a configuration file

Options:
[--reset] [--server] [--alias] [-k | --key] [--registry]

--reset     Restore the default config file template
--server    Use to add a server to the config (\`grok config add --alias alias --server url --key key\`)
--alias     Use in conjunction with the \`server\` option to set the server name
--key       Use in conjunction with the \`server\` option to set the developer key
--default   Use in conjunction with the \`server\` option to set the added server as default
--registry  Docker registry URL (default: registry.{server hostname})
`;

const HELP_CREATE = `
Usage: grok create [name]

Create a package:

grok create         Create a package in the current working directory
grok create <name>  Create a package in a folder with the specified name

Please note that the package name may only include letters, numbers, underscores, or hyphens

Options:
[--eslint] [--ide] [--js | --ts] [--test]

--eslint    Add a configuration for eslint
--ide       Add an IDE-specific configuration for debugging (vscode)
--js        Create a JavaScript package
--ts        Create a TypeScript package (default)
--test      Add tests support (TypeScript packages only)
`;

const HELP_PUBLISH = `
Usage: grok publish [host]

Uploads a package
Checks for errors before publishing — the package won't be published if there are any.

Options:
[--all] [--refresh] [--link] [--build] [--release] [--rebuild-docker] [--skip-docker-rebuild] [--skip-check] [-v | --verbose]

--all                  Publish all available packages (run in packages directory)
--refresh              Publish all available already loaded packages (run in packages directory)
--link                 Link the package to local packages
--build                Builds the package
--release              Publish package as release version
--rebuild-docker       Force rebuild Docker images locally before pushing to registry
--skip-docker-rebuild  Skip auto-rebuild when Dockerfile folder has changed
--skip-check           Skip check stage
--verbose              Print detailed output

Running \`grok publish\` is the same as running \`grok publish defaultHost --build --debug\`
`;

const HELP_CHECK = `
Usage: grok check <pluginFolder>

Options:
[-r | --recursive] [-v | --verbose]

--recursive       Check all packages in the current directory
--soft            Even if an error occurs, it doesn't throw an exception
--verbose         Print detailed output

Check package content (function signatures, import statements of external modules, etc.)
`;

const HELP_TEST = `
Usage: grok test

Options:
[--package] [--category] [--test] [--host] [--csv] [--gui] [--skip-build] [--skip-publish] [--link] [--catchUnhandled] [--report] [--record] [--verbose] [--platform] [--benchmark] [--stress-test] [--debug] [--all] [-r | --recursive] [--filter] [--parallel N]

--package           Specify a package name to run tests for
--category          Specify a category name to run tests for
--test              Specify a test name to run
--host              Host alias as in the config file
--csv               Save the test report in a CSV file
--gui               Launch graphical interface (non-headless mode)
--debug             Enables debug point on tests run (useless without gui mode)
--verbose           Show debug information
--retry --no-retry  Enables or disables browser reload after a failed test
--report            Report failed tests to audit, notifies package author (default=false)
--skip-build        Skip the package build step
--skip-publish      Skip the package publication step
--skip-puppeteer    Skip the Puppeteer/DG.Test pass; only run Playwright (for playwright-only test directories)
--skip-playwright   Skip the Playwright pass; only run Puppeteer/DG.Test
--skip-node         Skip the Node (browserless) pass; run all tests in the browser
--node-only         Run only tests annotated {node: true} headless under Node, no browser
--link  	        Link the package to local utils
--record            Records the test execution process in mp4 format
--platform          Runs only platform tests (applicable for ApiTests package only)
--core              Runs package & auto tests & core tests (core tests run only from DevTools package)
--benchmark   	    Runs tests in benchmark mode
--stress-test       Runs shuffled stress-test only
--all               Runs tests for all available packages(run in packages directory)
--recursive         Test all packages in the current directory (parallel, table output)
--filter            Filter packages by package.json fields (e.g. --filter "category:Cheminformatics")
--parallel N        Max parallel test jobs (default: 4)

Run package tests

Examples:
  grok test -r                                       Test all packages
  grok test -r --filter "category:Cheminformatics"   Test matching packages
  grok test -r --parallel 2 --host dev               Test with 2 jobs against dev
  grok test -r --skip-build --skip-publish            Test without rebuilding

See instructions:
https://datagrok.ai/help/develop/how-to/test-packages#local-testing
`;

const HELP_TESTALL = `
Usage: grok testall

Options:
[--packages] [--host] [--csv] [--gui] [--skip-build] [--skip-publish] [--link-package] [--catchUnhandled] [--report] [--record] [--verbose] [--benchmark] [--stress-test] [--order] [--tags] [--testRepeat] [--browsers-count] [--debug]

--packages          Specify a packages names to run tests for
--host              Host alias as in the config file
--csv               Save the test report in a CSV file
--gui               Launch graphical interface (non-headless mode)
--debug             Enables debug point on tests run (useless without gui mode) 
--catchUnhandled    Catch unhandled exceptions during test execution (default=true)
--report            Report failed tests to audit, notifies packages author (default=false)
--skip-build        Skip the packages build step
--skip-publish      Skip the packages publication step
--link-package  	  Link the packages to local utils
--record            Records the test execution process in mp4 format
--verbose           Prints detailed information about passed and skipped tests in the console
--core              Runs packages & core tests (applicable for DevTools packages only)
--benchmark   	    Runs tests in benchmark mode
--stress-test       Runs shuffled stress-test only
--order             Specify order for tests invocation
--tags              Filter tests by tag name for run
--testRepeat        Set amount of tests repeats
--browsers-count    Set amount of browsers for tests run

Run tests of all or specified packages 

See instructions:
https://datagrok.ai/help/develop/how-to/test-packages#local-testing
`;

const HELP_LINK = `
Usage: grok link

Links \`datagrok-api\`, all necessary libraries and packages for local development.
Uses \`npm link\` unless the --path option specified. 
By default, it links packages from the parent directory of the repository's root.

Options:
--dev               Links also dev dependencies
--path              Instead of npm link, sets dependencies in package.json to local
--repo-only         Links packages only from the current repository
--unlink            Unlinks packages and sets last versions instead of local path in package.json dependencies 
--verbose           Prints detailed information about linked packages  
--all               Links all available packages(run in packages directory)
`;

const HELP_DOCKER_GEN = `
Usage: grok docker-gen

Generate Celery Docker artifacts from annotated Python functions in the python/ directory.
Produces Dockerfile, tasks.yaml, and Celery entry point in dockerfiles/<name>/.

Options:
[-v | --verbose]

--verbose         Print detailed output
`;

const HELP_MIGRATE = `
Usage: grok migrate

Migrates legacy function tags into the meta.role field.

Example:
  tags: ['viewer', 'ml']
  ⟶
  meta: { role: 'viewer,ml' }
`;

const HELP_BUILD = `
Usage: grok build

Build a package in the current directory, or recursively build multiple packages.

Options:
[-r | --recursive] [-s | --silent] [--filter] [--no-incremental] [--parallel N] [-v | --verbose]

--recursive       Build all packages in the current directory
--silent          Skip confirmation prompt (for recursive builds)
--filter          Filter packages by package.json fields (e.g. --filter "category:Cheminformatics")
--no-incremental  Run a full build instead of the default incremental build
--parallel N      Max parallel builds (default: 4)
--verbose         Print detailed output

Examples:
  grok build                                          Build the current package
  grok build -r                                       Build all packages in the current directory
  grok build -r -s                                    Build all packages without confirmation
  grok build -r --filter "category:Cheminformatics"   Build only matching packages
  grok build -r --parallel 8                          Build with 8 parallel jobs
`;

// const HELP_MIGRATE = `
// Usage: grok migrate

// Switch to \`grok\` tools by copying your keys to the config
// file and converting your scripts in the \`package.json\` file
// `;

const HELP_RUN = `
Usage: grok run [host]

Build, publish, and open the package in the browser.

Runs \`grok build\`, publishes the package to the server, then opens the server in the default browser.

Options:
[-k | --key] [--release] [-v | --verbose]

--key       Developer key (overrides config)
--release   Publish as a release version (default: debug)
--verbose   Print detailed output

Examples:
  grok run              Build, publish to default server, and open browser
  grok run dev          Build, publish to 'dev' server alias, and open browser
  grok run https://my.datagrok.ai/api --key abc123
`;

const HELP_REPORT = `
Usage: grok report <subcommand> [args]

Manage Datagrok user error reports

Subcommands:
    fetch    Download a report zip from a managed instance (writes <stem>_meta.json next to it)
    read     Normalize a report (zip or json) into one JSON object on stdout
    resolve  Mark a report as resolved (needs the _meta.json written by fetch)
    ticket   Create a JIRA ticket for a report directly in JIRA (no dedup; the key is NOT
             written back to the report — prefer POST /reports/{id}/jira for that)
    comment  Add a comment to a JIRA ticket (--body <text> | --body-file <path>)
    label    Add labels to a JIRA ticket
    attach   Attach a file to a JIRA ticket

JIRA subcommands need JIRA_TOKEN (plus JIRA_USER for a user API token); ticket also needs
--project or $JIRA_PROJECT.

Read flags:
    --extract-screenshot <path>  Write the screenshot binary to <path>
    --extract-d42 <dir>          Unpack .d42 sidecar tables into <dir>
    --extract-client-log         Write a sibling <stem>_client_log.json

Examples:
  grok report fetch dev 1528             Download report #1528 from the 'dev' instance
  grok report read /tmp/report.zip       Print normalized JSON for a local zip
  grok report read /tmp/report.json      Print normalized JSON for a raw report.json
  grok report read dev 1528              Fetch + normalize report #1528 from 'dev'
  grok report read /tmp/report.zip --extract-screenshot ./shot.png
  grok report resolve dev 1528           Resolve report #1528 on the 'dev' instance
  grok report ticket dev <report-uuid>   Create a JIRA ticket for a report

The instance name must match a server alias in ~/.grok/config.yaml.
`;


export const HELP_KG = `
Usage: grok kg <verb> [options]

Validate and generate the knowledge graph: the type files under
core/docs/knowledge-graph (schema.yaml, nodes/*.yaml, edges/*.yaml) and the
home documents, the markdown files whose frontmatter carries a \`feature:\` or a
prefixed \`id:\` key. A home needs a name: \`name:\`, \`title:\`, or the first
\`#\` heading of the body. Type names and files are lower-dash-case (\`part-of\`,
\`customer-contact\`); graph labels are the upper-snake rendering (\`PART_OF\`),
derived by the build and never authored.

A home is a markdown document with frontmatter when the thing is a document
(features, scenarios, initiatives) and a YAML file with the same keys when it is a
record (concepts, people, teams, customers); YAML homes live inside
core/docs/knowledge-graph (concepts/, internal/), put the prose in \`description:\`
and must set \`name:\`. Homes are looked for in
    ${HOME_ROOTS.join('\n    ')}
skipping nodes/, edges/ and schema.yaml, node_modules, dist, build, .dart_tool,
.claude, .git, fixtures and __tests__ folders, and the Test Track files
(packages/UsageAnalysis/files and playwright-public: their \`feature:\` key still
means the area, until it migrates to \`covers:\`).

Verbs:
    check       The validation gate: type files against schema.yaml, home documents
                against their types, frontmatter references (edge keys, reference
                properties) against the other homes, and repo paths cited in
                frontmatter or in the body. Prose \`~id\` mentions and code markers
                are not read yet; they arrive with \`build\`. Exit 1 on errors.
    gen         Write kg.d.ts, the glossary.md tables and feature-tree.md, all inside
                core/docs/knowledge-graph. Refuses while check reports errors.
    build       Run the extractors and write the graph as JSONL into one immutable
                generation, .kg/gen/<batch>/: data/nodes/<type>.jsonl,
                data/edges/<type|property>.jsonl, reports/ (claims.jsonl for membership,
                invalid.jsonl, problems.json), the index as kg.kuzu, and manifest.json
                last of all. Only when the generation is complete does .kg/current, one
                line naming the batch, start pointing at it (with a .kg/gen/current link
                beside it where the platform allows one), so a build that is interrupted
                or fails to load the index leaves the previous generation queryable.
                Earlier generations are never removed by build; grok kg gc removes them.
                Lines are deterministic: two builds of the same inputs are byte-identical
                and share a content-addressed batch id over both revisions, the dirty tree
                of both repositories, the schema and builder versions, the mode, the
                extractor selection, the backlog snapshot and the Dart batch; only the
                manifest carries the time.
                Problems are counted in the manifest, never thrown. Extractors today:
                homes (the home-document layer), ts-packages (packages, libraries,
                semantic types and depends-on from package.json) and ts-functions
                (registered functions, scripts, queries, connections, environments,
                containers and by-name calls under public/packages), ts-declarations,
                ts-imports and ts-uses (source files, declarations, extends/implements,
                resolved imports and JS API usage over js-api, packages and libraries),
                ts-tests (DG and Playwright tests in their suites), ts-samples (ApiSamples),
                ts-changelog (CHANGELOG.md bullets), docs (markdown pages, headings,
                mentions, legacy Test Track scenarios, tutorials), process (backlog
                tickets, release records with their picked commits, and people) and
                membership (which feature owns each file, conventions.md §8; reports
                ownership.json).
                Loading the index is part of build: the JSONL is copied into a Kuzu
                database at .kg/gen/<batch>/kg.kuzu (one node table per root, one rel
                table per edge type and per reference property), and the manifest records
                it as indexed_batch with the memory the load needed (index_memory_mb) and
                the platform it was written on (index_platform). The binding is optional —
                without it build says so in one line and still succeeds. With --no-db the
                generation carries no index, and current stays where it is rather than
                take the index away from query.
    query       Cypher over the built index: grok kg query "MATCH (n:Feature) RETURN n.id",
                or --file q.cypher. Exit 2 when kuzu is not installed.
    impact      What a change reaches: the features that own or take part in a file,
                declaration or feature, their owners, tests, docs and tickets, and for a
                declaration or a file, who calls or imports it.
    tests-for   The tests, scenarios and automations of a feature (with everything under
                it in the tree), or of the feature that owns a path.
    explain     One node: its properties, then every one-hop edge by type and direction.
    find        The vocabulary search over ids, names, aliases, descriptions and keywords
                of all eight tables; features and concepts first.
    report      One maintainer report over the JSONL a build already wrote, never the
                index: orphans (files with no owner, grouped by package or core
                sub-project), stale (citations, tickets, help-urls, specs and
                declarations the graph can no longer reach), coverage (one row per
                feature: owner, tests, scenarios, documents, description), proposed
                (folders with code and no owner, with the id they would take) and
                diff (the features a branch touches and the tests that cover them,
                \`--diff <ref>\`; the md output is the PR comment). \`build\` writes the
                first four to .kg/reports/ as both .json and .md.
    gc          Remove older generations under .kg/gen/, keeping the current one and
                the --keep newest (default 2). A generation whose index a reader holds
                open is reported and left alone.
    help        Show this help

The four operations and query read the generation .kg/current names, next to the type
files, and refuse an index that was loaded from another batch; each of them prints
\`Dart coverage unknown (no kg-dart batch)\` first while no prop_gen batch has been
built. Ids may be written with or without the \`~\` sigil.

\`check\` gates the sources; \`gen --check\` gates the generated files, failing when
kg.d.ts, glossary.md or feature-tree.md on disk differ from what gen would write.

Options:
    --kg <dir>          The knowledge-graph folder (default: found by walking up from
                        the current directory to the monorepo root)
    --types-only        Check the type files only, skip the home documents (check only)
    --check             With gen: fail if the generated files differ from disk, write nothing
    --output <format>   table (default) or json (the check report, or the build manifest);
                        query and the operations also take csv, report takes md
    --file <path>       With query: read the Cypher from a file
    --limit <n>         With the operations: rows per section (default 50)
    --quiet             Print errors only: no warnings, no summary line (check, gen)
    --public            With build: the public projection (public node types, visibility
                        public, no home or owner, edges with both ends public) into public/.kg/
    --only <a,b>        With build: run only the named extractors
    --backlog <dir>     With build: the backlog snapshot repo (used by the process layer)
    --no-db             With build: write the JSONL only, do not load the graph index
    --out <dir>         With build, report and gc: write (or read) under <dir> instead of .kg/
    --keep <n>          With gc: generations to keep besides the current one (default 2)
    --memory <mb>       Kuzu buffer pool in MB: the load takes 2048, a reader 512 or what
                        the manifest says the load needed. KG_KUZU_MEMORY does the same.
    --diff <ref>        With report diff: the revision HEAD is compared against

Examples:
  grok kg check
  grok kg check --types-only --output json
  grok kg gen
  grok kg gen --check
  grok kg build --only homes --no-db
  grok kg build --public --output json
  grok kg gc --keep 1
  grok kg query "MATCH (f:Feature)-[:owner]->(p:Actor) RETURN f.id, p.id"
  grok kg impact core/server/datlas/lib/src/services/spaces_service.dart
  grok kg tests-for ~visualize/viewers/scatter-plot --output json
  grok kg explain ~govern/spaces
  grok kg find scatter
  grok kg report coverage --output json
  grok kg report diff --diff master --output md

The contract is core/docs/knowledge-graph/conventions.md.
`;

export const help = {
  add: HELP_ADD,
  api: HELP_API,
  build: HELP_BUILD,
  check: HELP_CHECK,
  claude: HELP_CLAUDE,
  config: HELP_CONFIG,
  create: HELP_CREATE,
  'docker-gen': HELP_DOCKER_GEN,
  init: HELP_INIT,
  link: HELP_LINK,
  publish: HELP_PUBLISH,
  report: HELP_REPORT,
  run: HELP_RUN,
  test: HELP_TEST,
  testall: HELP_TESTALL,
  migrate: HELP_MIGRATE,
  server: HELP_SERVER,
  s: HELP_SERVER,
  kg: HELP_KG,
  help: HELP,
};
