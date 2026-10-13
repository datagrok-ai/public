// EMS's components only; the format's preamble is in u2's manifest. `--check` fails on a stale file.
import {readFileSync, writeFileSync} from 'node:fs';
import {fileURLToPath} from 'node:url';
import '@datagrok-libraries/u2/tests/dom-shim.js';
import {Registry} from '@datagrok-libraries/u2';
import {registerDomainSources} from '../dist/src/core/registrations.js';

const registry = new Registry();
registerDomainSources(registry);
const {components} = registry.manifest();
const manifest = {components};
const path = fileURLToPath(new URL('../manifest.json', import.meta.url));
const json = JSON.stringify(manifest, null, 2) + '\n';

if (!process.argv.includes('--check')) {
  writeFileSync(path, json);
  console.log(`${manifest.components.length} components → ${path}`);
} else {
  const lf = (s) => s.replace(/\r\n/g, '\n');
  const actual = lf(readFileSync(path, 'utf8'));
  if (actual === lf(json))
    console.log(`manifest.json is up to date (${manifest.components.length} components)`);
  else {
    console.error('manifest.json is stale — run `npm run manifest` and commit the result.');
    process.exit(1);
  }
}
