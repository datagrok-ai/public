/* Writes manifest.json from the registry — the machine-readable component list an LLM reads
   before emitting a dg-ui/1 spec. Run `npm run build` first: this loads the compiled sources.
   With --check it compares instead of writing, so a stale checked-in file fails the build. */
import {readFileSync, writeFileSync} from 'node:fs';
import {fileURLToPath} from 'node:url';
import '../tests/dom-shim.js';
import {registry} from '../src/spec/registry.js';
import {registerAll} from '../src/spec/registrations.js';

registerAll();
const manifest = registry.manifest();
const path = fileURLToPath(new URL('../manifest.json', import.meta.url));
const json = JSON.stringify(manifest, null, 2) + '\n';

if (!process.argv.includes('--check')) {
  writeFileSync(path, json);
  console.log(`${manifest.components.length} components → ${path}`);
} else {
  /* the checked-in copy is CRLF in a Windows checkout, the generated string always LF */
  const lf = (s) => s.replace(/\r\n/g, '\n');
  const actual = lf(readFileSync(path, 'utf8'));
  if (actual === lf(json))
    console.log(`manifest.json is up to date (${manifest.components.length} components)`);
  else {
    const a = actual.split('\n');
    const b = lf(json).split('\n');
    let i = 0;
    while (i < a.length && i < b.length && a[i] === b[i]) i++;
    console.error('manifest.json is stale — run `npm run manifest` and commit the result.');
    console.error(`  first difference at line ${i + 1}:`);
    console.error(`    checked in: ${i < a.length ? a[i] : '<end of file>'}`);
    console.error(`    generated:  ${i < b.length ? b[i] : '<end of file>'}`);
    process.exit(1);
  }
}
