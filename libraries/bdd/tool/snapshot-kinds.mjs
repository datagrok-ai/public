/* Writes tests/fixtures/kinds-before.json: every kind key (names and aliases) of the built
   bindings/common/kinds.ts → the entry it resolves to. The compat test in tests/vocabulary.test.ts
   proves a rewritten kind table still has every key, selector, strategy, part and Dart name of this
   one. Run `npm run build` first: this loads dist/. */
import {writeFileSync} from 'node:fs';
import {fileURLToPath} from 'node:url';

const {kindNames, lookupKind} = await import(new URL('../dist/src/registry.js', import.meta.url).href);
await import(new URL('../dist/bindings/common/kinds.js', import.meta.url).href);

const FIELDS = ['name', 'selector', 'match', 'labelSelector', 'dartNames', 'editorSelector', 'parts', 'gestures'];
const table = {};
for (const key of [...kindNames()].sort()) {
  const entry = lookupKind(key);
  table[key] = Object.fromEntries(FIELDS.filter((f) => entry[f] !== undefined).map((f) => [f, entry[f]]));
}

const path = fileURLToPath(new URL('../tests/fixtures/kinds-before.json', import.meta.url));
writeFileSync(path, JSON.stringify(table, null, 2) + '\n');
console.log(`${Object.keys(table).length} kind keys → ${path}`);
