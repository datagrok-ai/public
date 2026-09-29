/// Command-line parsing for `grok`: minimist, plus the `grok s` rules for values it would misread.
// eslint-disable-next-line @typescript-eslint/no-require-imports
const minimist = require('minimist');

/** Options whose value may start with '-': `--save-levels -audit`, `--since -7d`. */
const DASH_VALUE_FLAGS = ['--print-levels', '--post-levels', '--save-levels', '--debug-flags', '--from', '--to', '--since'];

const BOOLEAN_FLAGS = ['dartium'];

export function parseArgs(args: string[]): any {
  const first = String(minimist(args, {boolean: BOOLEAN_FLAGS})._[0]);
  const server = first === 's' || first === 'server';
  const joined: string[] = [];
  for (const arg of args) {
    const prev = joined[joined.length - 1];
    // minimist reads a value that starts with '-' as more flags: keep it with its option
    if (server && DASH_VALUE_FLAGS.includes(prev) && /^-[^-]/.test(arg))
      joined[joined.length - 1] = `${prev}=${arg}`;
    else
      joined.push(arg);
  }
  return minimist(joined, {
    alias: {k: 'key', h: 'help', r: 'recursive'},
    boolean: BOOLEAN_FLAGS,
    // keep versions and ids verbatim — minimist would coerce '1.10' to 1.1 and the id prefix
    // '08327370' to 8327370; `grok s` positionals are ids, names and paths, never numbers
    string: ['version', 'until-version', 'signature', 'action', 'request', 'session', 'override', 'rule', ...(server ? ['_'] : [])],
  });
}
