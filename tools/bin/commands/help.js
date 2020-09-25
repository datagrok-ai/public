const HELP = `
Usage: grok <command>

Datagrok's package management tool

Commands:
    add       Add a script, app, or viewer
    config    Create and manage config files
    create    Create a package
    delete    Delete a package
    publish   Upload a package
    migrate   Switch to \`grok\` tools

To get help on a particular command, use:
    grok <command> --help
`;

const HELP_ADD = `
Usage: grok add <entity> <language>* <name>

Add an object template to your package:

grok add app <name>
grok add function <name>
grok add script <language> <name>
grok add view <name>
grok add viewer <name>

Please note that entity names may only include letters and numbers
`;

const HELP_CONFIG = `
Usage: grok config

Create or update a configuration file

Options:
[--reset]

--reset     Restore the default config file template
`;

const HELP_CREATE = `
Usage: grok create <name>*

Create a package
Please note that the package name may only include letters, numbers, underscores, or hyphens
`;

const HELP_DELETE = `
Usage: grok delete <name>*

Delete a package
`;

const HELP_PUBLISH = `
Usage: grok publish <host>*

Upload a package

Options:
[--build|--rebuild] [--debug|--release] [--key]

Running \`grok publish\` is the same as running \`grok publish defaultHost --build --debug\`
`;

const HELP_MIGRATE = `
Usage: grok migrate

Switch to \`grok\` tools by copying your keys to the config
file and converting your scripts in the \`package.json\` file
`;

module.exports = {
    add: HELP_ADD,
    config: HELP_CONFIG,
    create: HELP_CREATE,
    delete: HELP_DELETE,
    publish: HELP_PUBLISH,
    migrate: HELP_MIGRATE,
    help: HELP
};
