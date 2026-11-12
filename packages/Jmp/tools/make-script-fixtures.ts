// Writes the synthetic saved-analysis fixtures in files/tests/ (JMP_PLUGIN.md §9.2). Each is
// files/jmpreader/compact_UInt8.jmp with its one saved script ("Source") replaced by another of the same
// length, so no offset in the file moves. They test what the package does with a saved analysis, not
// where JMP keeps one: no fixture made by JMP itself has an analysis yet.
//
//   npx tsx tools/make-script-fixtures.ts

import {readFileSync, writeFileSync} from 'node:fs';
import {join} from 'node:path';

const FIXTURES: Record<string, string> = {
  'saved-bar-chart': 'Bar chart of A(Graph Builder(Variables(X(:A)), Elements(Bar(X, Legend(2)))))',
  'saved-fit-model': 'Fit A(Fit Model(Y(:A), Effects(:A), Personality("Nominal Logistic"), Run))',
};

const dir = join(__dirname, '..', 'files');
const source = readFileSync(join(dir, 'jmpreader', 'compact_UInt8.jmp'));
const at = source.indexOf('Source(');
const len = source.readUInt32LE(at - 4);
if (at < 0 || source[at - 5] !== 0x03 || source[at + len - 1] !== 0x29)
  throw new Error('compact_UInt8.jmp no longer has the expected Source script');

for (const [name, script] of Object.entries(FIXTURES)) {
  // pad inside the outer call, which JSL ignores
  const padded = script.slice(0, -1) + ' '.repeat(len - script.length) + ')';
  if (Buffer.byteLength(padded) !== len)
    throw new Error(`${name}: the script is longer than the ${len} bytes available`);
  const out = Buffer.from(source);
  out.write(padded, at, 'utf8');
  writeFileSync(join(dir, 'tests', `${name}.jmp`), out);
  console.log(`files/tests/${name}.jmp`);
}
