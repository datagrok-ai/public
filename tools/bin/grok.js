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
// The help texts are a large module; load them only when one is printed.
let _help;
const help = new Proxy({}, {get: (_, key) => (_help ??= require('./commands/help').help)[key]});
const runAllCommand = require('./utils/utils').runAll;

// Each command module is loaded only when invoked: loading all of them (puppeteer, ts-morph,
// archiver, inquirer) used to cost ~6 s of start-up on every `grok api` / `grok check`.
const lazy = (file, name) => (args) => require(`./commands/${file}`)[name](args);
const commands = {
  add: lazy('add', 'add'),
  api: lazy('api', 'api'),
  build: lazy('build', 'build'),
  check: lazy('check', 'check'),
  claude: lazy('claude', 'claude'),
  config: lazy('config', 'config'),
  create: lazy('create', 'create'),
  'docker-gen': lazy('docker-gen', 'dockerGen'),
  init: lazy('init', 'init'),
  link: lazy('link', 'link'),
  publish: lazy('publish', 'publish'),
  report: lazy('report', 'report'),
  run: lazy('run', 'run'),
  test: lazy('test', 'test'),
  tsc: lazy('tsc', 'tsc'),
  testall: lazy('test-all', 'testAll'),
  stresstest: lazy('stress-tests', 'stressTests'),
  migrate: lazy('migrate', 'migrate'),
  server: lazy('server', 'server'),
  s: lazy('server', 'server'),
  setup: lazy('setup', 'setup'),
};

// `--version` is a string option (grok publish --version 1.10), so a bare `grok --version` parses as ''.
if (argv._.length === 0 && ('version' in argv && argv.version === '' || argv.v === true)) {
  console.log(require('../package.json').version);
  process.exit(0);
}

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
