// Writes a large synthetic .jmp file and times reading it (JMP_PLUGIN.md §8). The file is
// files/jmpreader/compressed.jmp refilled with varied rows (src/tests/large-jmp.ts); it isn't committed.
//
//   npm run large -- [--mb 100] [--out <file>]

import {readFileSync, writeFileSync} from 'node:fs';
import {join} from 'node:path';

import {makeLargeJmp, rowsForSize} from '../src/tests/large-jmp';
import {readJmp} from '../src/formats/jmp/jmp-parser';

async function main(): Promise<void> {
  const args = process.argv.slice(2);
  const arg = (name: string): string | undefined => {
    const i = args.indexOf(name);
    return i >= 0 ? args[i + 1] : undefined;
  };
  const mb = Number(arg('--mb') ?? 100);
  const out = arg('--out') ?? `large-${mb}mb.jmp`;
  const template = new Uint8Array(readFileSync(join(__dirname, '..', 'files', 'jmpreader', 'compressed.jmp')));

  let t = performance.now();
  const bytes = await makeLargeJmp(template, rowsForSize(mb * 1e6));
  writeFileSync(out, bytes);
  const seconds = (): string => ((performance.now() - t) / 1000).toFixed(2);
  console.log(`wrote ${out}: ${(bytes.length / 1e6).toFixed(1)} MB in ${seconds()} s`);

  global.gc?.();
  const used = (): number => process.memoryUsage().heapUsed + process.memoryUsage().arrayBuffers;
  const before = used();
  t = performance.now();
  const table = await readJmp(bytes);
  const elapsed = seconds();
  const rate = (table.info.nrows * table.info.ncols / (performance.now() - t) / 1000).toFixed(1);
  global.gc?.();
  const grown = ((used() - before) / 1e6).toFixed(0);
  console.log(`read ${table.info.nrows.toLocaleString()} rows × ${table.info.ncols} columns in ${elapsed} s ` +
    `(${rate}M cells/s); memory held by the table (heap + array buffers) +${grown} MB`);
  const bad = table.columns.filter((c) => c.decodeError);
  if (bad.length > 0)
    throw new Error(`columns not read: ${bad.map((c) => `${c.info.name}: ${c.decodeError}`).join('; ')}`);
}

void main();
