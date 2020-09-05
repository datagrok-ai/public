#! /usr/bin/env node
const argv = require('minimist')(process.argv.slice(2));
const fs = require('fs');
const inquirer = require('inquirer');
const os = require('os');
const path = require('path');
const yaml = require('js-yaml');

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
`;

const HELP_CREATE = `
Usage: grok create <name>

Create a package
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

const CONFIG = {
    servers: {
        dev: { url: 'https://dev.datagrok.ai', key: '' },
        public: { url: 'https://public.datagrok.ai', key: '' },
        local: { url: 'https://localhost:8082', key: '' }
    },
    default: 'public'
};

const GROK_DIR = path.join(os.homedir(), '.grok');
const CONF_PATH = path.join(GROK_DIR, 'config.yaml')

function validateKey(key) {
    if (/^([A-Za-z\d\-])+$/.test(key)) {
        return true;
    } else {
        return 'Developer key may only include letters, numbers or hyphens';
    }
};

function generateKeyQ(server) {
    let question =  {
        name: server,
        type: 'input',
        message: `Developer key (get it from https://${server}.datagrok.ai/u):`,
        validate: validateKey
    };
    if (server === 'local') {
        question.message = 'Developer key for https://localhost:8082:'
    }
    return question;
};

let nOptions = Object.keys(argv).length - 1;
let command = argv['_'][0];

switch (command) {
    case 'add':
        if (false) {
            // If all args and ops are correct, do ...
        } else {
            console.log(HELP_ADD);
        }
        break;
    case 'config':
        if (argv['_'].length === 1 && (nOptions < 1 || nOptions === 1 && argv.reset)) {
            if (!fs.existsSync(GROK_DIR)) {
                fs.mkdirSync(GROK_DIR);
            }
            if (!fs.existsSync(CONF_PATH) || argv.reset) {
                fs.writeFileSync(CONF_PATH, yaml.safeDump(CONFIG));
            }
            let config = yaml.safeLoad(fs.readFileSync(CONF_PATH));
            console.log(`Your config file (${CONF_PATH}):`);
            console.log(config);
            (async () => {
                let answers = await inquirer.prompt({
                    name: 'edit-config',
                    type: 'confirm',
                    message: 'Do you want to edit it?',
                    default: false,
                })
                if (answers['edit-config']) {
                    for (server in config.servers) {
                        let question = generateKeyQ(server);
                        question.default = config['servers'][server]['key'];
                        let devKey = await inquirer.prompt(question);
                        config['servers'][server]['key'] = devKey[server];
                    }
                    let defaultServer = await inquirer.prompt({
                        name: 'default-server',
                        type: 'input',
                        message: 'Your default server:',
                        validate: function (server) {
                            if (server in config.servers) return true;
                            else return 'Only one of the specified servers may be chosen as default';
                        },
                        default: config.default,
                    })
                    config.default = defaultServer['default-server'];
                    fs.writeFileSync(CONF_PATH, yaml.safeDump(config));
                }
            })();
        } else {
            console.log(HELP_CONFIG);
        }
        break;
    case 'create':
        if (false) {
            // If all args and ops are correct, do ...
        } else {
            console.log(HELP_CREATE);
        }
        break;
    case 'delete':
        if (false) {
            // If all args and ops are correct, do ...
        } else {
            console.log(HELP_DELETE);
        }
        break;
    case 'publish':
        if (false) {
            // If all args and ops are correct, do ...
        } else {
            console.log(HELP_PUBLISH);
        }
        break;
    default:
        console.log(HELP);
};
