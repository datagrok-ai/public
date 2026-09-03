#!/usr/bin/env node
const argv = require('minimist')(process.argv.slice(2), {
  alias: {k: 'key', h: 'help', r: 'recursive'},
  boolean: ['dartium'],
  // keep versions verbatim — minimist would coerce '1.10' to the number 1.1
  string: ['version'],
});
// minimist maps `--no-retry` to `{retry: false}`, so the `args['no-retry']` checks in
// test.ts / playwright-runner.ts never fired and `--no-retry` was silently ignored
// (Playwright kept retrying failed specs). Normalize back to the flag the commands read.
if (argv.retry === false) argv['no-retry'] = true;
const help = require('./commands/help').help;
const runAllCommand = require('./utils/utils').runAll;

const commands = {
  add: require('./commands/add').add,
  api: require('./commands/api').api,
  build: require('./commands/build').build,
  check: require('./commands/check').check,
  claude: require('./commands/claude').claude,
  config: require('./commands/config').config,
  create: require('./commands/create').create,
  'docker-gen': require('./commands/docker-gen').dockerGen,
  init: require('./commands/init').init,
  link: require('./commands/link').link,
  publish: require('./commands/publish').publish,
  report: require('./commands/report').report,
  run: require('./commands/run').run,
  test: require('./commands/test').test,
  testall: require('./commands/test-all').testAll,
  stresstest: require('./commands/stress-tests').stressTests,
  migrate: require('./commands/migrate').migrate,
  schema: require('./commands/schema').schema,
  server: require('./commands/server').server,
  s: require('./commands/server').server,
};

const onPackageCommandNames = ['api', 'check', 'link', 'publish', 'test'];

// A machine-readable run prints its error as JSON on stderr (server.ts) and nothing else:
// a usage dump on stdout would corrupt the output the caller parses.
const outputFormat = argv.output ?? argv.o;
function printUsage(command) {
  if (outputFormat !== 'json')
    process.stderr.write(`${help[command]}\n`);
}

const command = argv['_'][0];
if (command !== 'test' && command !== 'stresstest')
  delete argv.dartium;
if (command in commands) {
  try {
    if (argv['help']) {
      console.log(help[command]);
      exitWithCode(1);
    } else if (argv.all && onPackageCommandNames.includes(command)) {
      runAllCommand(process.cwd(),
        `grok ${process.argv.slice(2).join(' ')}`.replace('--all', ''), {});
    } else {
      const result = commands[command](argv);
      if (result && typeof result.then === 'function') {
        result.then((ok) => {
          if (!ok) {
            printUsage(command);
            exitWithCode(1);
          }
        }).catch((err) => {
          console.error(err);
          printUsage(command);
          exitWithCode(255);
        });
      }
      else if (!result) {
        printUsage(command);
        exitWithCode(1);
      }
    }
  } catch (err) {
    console.error(err);
    printUsage(command);
    exitWithCode(255);
  }
} else
  console.log(help.help);


function exitWithCode(code) {
  if (outputFormat !== 'json')
    console.log(`Exiting with code ${code}`);
  process.exit(code);
}
