#!/usr/bin/env node
/* Launcher: prefers the copy of the library installed in the current package (its bindings and
   generated specs import that one), falls back to this global install. */
import {existsSync} from 'node:fs';
import {createRequire} from 'node:module';
import {dirname, resolve} from 'node:path';
import {fileURLToPath, pathToFileURL} from 'node:url';

// tsx imports tsconfig/package JSON through the ESM loader; Node flags that as experimental on every run
const emitWarning = process.emitWarning.bind(process);
process.emitWarning = (warning, ...rest) => {
  if (!String(warning).includes('Importing JSON modules'))
    emitWarning(warning, ...rest);
};

let lib = resolve(dirname(fileURLToPath(import.meta.url)), '..');
try {
  lib = dirname(createRequire(resolve(process.cwd(), 'package.json')).resolve('@datagrok-libraries/bdd/package.json'));
} catch {
  // no local install — the global one serves
}
const cli = resolve(lib, 'dist', 'src', 'cli.js');
if (!existsSync(cli)) {
  console.error(`${cli} is missing — run \`npm run build\` in ${lib}`);
  process.exit(2);
}
await import(pathToFileURL(cli).href);
