#!/usr/bin/env node
const argv = require('minimist')(process.argv.slice(2), {
  alias: {k: 'key', h: 'help', r: 'recursive'},
});
const help = require('./commands/help').help;
const runAllCommand = require('./utils/utils').runAll;

const commands = {
  add: require('./commands/add').add,
  api: require('./commands/api').api,
  check: require('./commands/check').check,
  config: require('./commands/config').config,
  create: require('./commands/create').create,
  init: require('./commands/init').init,
  link: require('./commands/link').link,
  publish: require('./commands/publish').publish,
  test: require('./commands/test').test,
  testall: require('./commands/test-all').testAll,
};

const onPackageCommandNames = ['api', 'check', 'link', 'publish', 'test'];

const command = argv['_'][0];
if (command in commands) {
  try {
    if (argv['help']) {
      console.log(help[command]);
      exitWithCode(1);
    } else if (argv.all && onPackageCommandNames.includes(command)) {
      runAllCommand(process.cwd(),
        `grok ${process.argv.slice(2).join(' ')}`.replace('--all', ''), {});
    } else if (!commands[command](argv)) {
      console.log(help[command]);
      exitWithCode(1);
    }
  } catch (err) {
    console.error(err);
    console.log(help[command]);
    exitWithCode(255);
  }
} else
  console.log(help.help);


function exitWithCode(code) {
  console.log(`Exiting with code ${code}`);
  process.exit(code);
}
