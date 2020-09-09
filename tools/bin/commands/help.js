const HELP = `
Usage: grok <command>

Datagrok's package management tool

Commands:
    add       Add a script, app, or viewer
    config    Create and manage config files
    create    Create a package
    delete    Delete a package
    publish   Upload a package

To get help on a particular command, use:
    grok <command> --help
`;

const HELP_ADD = `
Usage: grok add <entity> <language>* <name>

Add an object template to your package

grok add script r <name>
grok add script python <name>
`;

const HELP_CONFIG = `
Usage: grok config

Create or update a configuration file

Options:
[--reset]

--reset     Restore the default config file template
`;

const HELP_CREATE = `
Usage: grok create <name>

Create a package. Please note that the package name may only include letters, numbers, underscores, or hyphens
`;

const HELP_DELETE = `
Usage: grok delete <name>

Delete a package
`;

const HELP_PUBLISH = `
Usage: grok publish <host>

Upload a package

Options:
[--build|--rebuild] [--debug|--release]

Running \`grok publish\` is the same as running \`grok publish --build --debug\`
`;

module.exports = {
    add: HELP_ADD,
    config: HELP_CONFIG,
    create: HELP_CREATE,
    delete: HELP_DELETE,
    publish: HELP_PUBLISH,
    help: HELP
};
