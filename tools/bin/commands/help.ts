const HELP = `
Usage: grok <command>

Datagrok's package management tool

Commands:
    add       Add an object template
    api       Create wrapper functions
    config    Create and manage config files
    create    Create a package
    publish   Upload a package
    migrate   Switch to \`grok\` tools

To get help on a particular command, use:
    grok <command> --help

Read more about the package development workflow:
https://datagrok.ai/help/develop/develop
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

const HELP_API = `
Usage: grok api

Create wrapper functions for package scripts and queries
`;

const HELP_CONFIG = `
Usage: grok config

Create or update a configuration file

Options:
[--reset] [--server] [--alias] [--key]

--reset     Restore the default config file template
--server    Use to add a server to the config (\`grok config add --alias alias --server url --key key\`)
--alias     Use in conjunction with the \`server\` option to set the server name
--key       Use in conjunction with the \`server\` option to set the developer key
--default   Use in conjunction with the \`server\` option to set the added server as default
`;

const HELP_CREATE = `
Usage: grok create [name]

Create a package:

grok create         Create a package in the current working directory
grok create <name>  Create a package in a folder with the specified name

Please note that the package name may only include letters, numbers, underscores, or hyphens

Options:
[--eslint] [--ide] [--js|--ts] [--jest]

--eslint    Add a configuration for eslint
--ide       Add an IDE-specific configuration for debugging (vscode)
--js        Create a JavaScript package
--ts        Create a TypeScript package (default)
--jest      Add a configuration for jest
`;

const HELP_PUBLISH = `
Usage: grok publish [host]

Upload a package

Options:
[--build|--rebuild] [--debug|--release] [--key] [--suffix]

Running \`grok publish\` is the same as running \`grok publish defaultHost --build --debug\`
`;

const HELP_MIGRATE = `
Usage: grok migrate

Switch to \`grok\` tools by copying your keys to the config
file and converting your scripts in the \`package.json\` file
`;

export const help = {
  add: HELP_ADD,
  api: HELP_API,
  config: HELP_CONFIG,
  create: HELP_CREATE,
  publish: HELP_PUBLISH,
  migrate: HELP_MIGRATE,
  help: HELP
};
