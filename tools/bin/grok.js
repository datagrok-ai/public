#!/usr/bin/env node
const argv = require('minimist')(process.argv.slice(2), {alias: { k: 'key'}});
const help = require('./commands/help');

const commands = {
    add: require('./commands/add').add,
    config: require('./commands/config').config,
    create: require('./commands/create').create,
    delete: require('./commands/delete').delete,
    publish: require('./commands/publish').publish
};

const command = argv['_'][0];
if (command in commands) {
    try {
        if (!commands[command](argv)) {
            console.log(help[command]);
        }
    } catch (err) {
        console.error(err);
        console.log(help[command]);
    }
} else {
    console.log(help.help);
}
