import { migrate } from "./migrate";
import { testAll } from "./test-all";

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
grok add connection <name>
grok add detector <semantic-type-name>
grok add function [tag] <name>
grok add query <name>
grok add script [tag] <language> <name>
grok add view <name>
grok add viewer <name>
grok add tests

Please note that entity names may only include letters and numbers

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

Create wrapper functions for package scripts and queries

Options:
[-v | --verbose]

--verbose         Print detailed output
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
    fetch    Download a report zip from a managed instance
    read     Normalize a report (zip or json) into one JSON object on stdout
    resolve  Mark a report as resolved
    ticket   Create a JIRA ticket for a report via the Datlas API

Read flags:
    --extract-screenshot <path>  Write the screenshot binary to <path>
    --extract-d42 <dir>          Unpack .d42 sidecar tables into <dir>
    --extract-actions            Write a sibling <stem>_actions.json

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

const HELP_SERVER = `
Usage: grok server <entity> <verb> [args] [options]
       grok s <entity> <verb> [args] [options]

Manage a Datagrok server from the command line.

Entities:
  users, groups, functions, connections, queries, scripts, packages, reports, files, tables

Verbs:
  list      List entities
  get       Get a single entity by ID or name
  delete    Delete an entity by ID

Special commands:
  grok s functions run <Name:func(args)>             Call a function
  grok s functions list [--type <t>] [--language <l>] [--package <p>] [--filter <expr>]
                                                     Type: script|query|function|package; language applies to scripts
  grok s files list [path] [-r]                      List files (recursive with -r)
  grok s shares add <entity> <group>[,<group>...] [--access View|Edit]
                                                     Share an entity with groups
  grok s shares list <entity-id>                     List who an entity (UUID) is shared with
  grok s users save --json user.json                 Create or update a user from JSON
  grok s users block <id-or-login>                   Block a user from the platform
  grok s users unblock <id-or-login>                 Unblock a previously blocked user
  grok s groups save --json group.json [--save-relations]
                                                     Create or update a group from JSON
  grok s connections save --json conn.json [--save-credentials]
                                                     Create or update a connection from JSON
  grok s connections test <id-or-name>               Test connectivity of an existing connection
  grok s connections test --json conn.json           Test connectivity of a connection defined in JSON
  grok s tables upload <name> <file.csv>             Upload a CSV as a Datagrok table
  grok s tables download <name-or-id> [-O <file>]    Download a table as CSV (stdout by default)
  grok s raw <METHOD> <path>                         Hit any API endpoint
  grok s describe <entity-type>                      Show entity JSON schema

Options:
  --host <alias|url>    Server alias from config or full URL
  --output <format>     Output format: table (default), json, csv, quiet
  --filter <text>       Smart filter expression
  --limit <n>           Page size (default: 50)
  --offset <n>          Start offset (default: 0)
  -r, --recursive       Recursive (for files list)
  --json <file>         Read function parameters from JSON file

Examples:
  grok s users list
  grok s connections list --filter "PostgreSQL" --output json
  grok s connections get <id>
  grok s connections delete <id>
  grok s connections save --json conn.json --save-credentials
  grok s connections test "JohnDoe:MyConnection"
  grok s connections test --json conn.json
  grok s users save --json user.json
  grok s groups save --json group.json --save-relations
  grok s shares add "JohnDoe:MyConnection" Chemists,Admins --access Edit
  grok s shares list <entity-uuid>
  grok s functions run 'Chem:smilesToMw("ccc")'
  grok s files list "System:AppData" -r
  grok s raw GET /api/users/current
  grok s describe connections
  grok s users list --host dev
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
  help: HELP,
};
