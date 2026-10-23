// Reads every .jmp file under a folder with the package's parser and reports what came across
// (JMP_PLUGIN.md §9.4). For corpora that can't be committed: a JMP install's Samples/Data, or customer
// files read on the machine they were shared to.
//
//   npm run corpus -- <folder> [--golden <folder>]
//
// With --golden, each file is also compared with <golden>/<name>.json written by tools/jmp2json.jl
// over the same folder.
//
// It also reports what the format work still needs (JMP_PLUGIN.md §4.5, §11): the saved scripts found,
// which JMP platforms they launch and how many convert to viewers, and row-state columns whose
// undecoded bytes are set, which is where JMP's excluded, hidden and selected flags may be.

import {readdirSync, readFileSync, statSync, existsSync} from 'node:fs';
import {basename, join} from 'node:path';

import {readJmp} from '../src/formats/jmp/jmp-parser';
import {compareToGolden, Golden} from '../src/tests/golden';
import {convertScripts} from '../src/formats/jmp/jsl-to-viewers';
import {isSourceScript} from '../src/formats/jmp/jmp-scripts';

function walk(dir: string): string[] {
  return readdirSync(dir).flatMap((name) => {
    const path = join(dir, name);
    return statSync(path).isDirectory() ? walk(path) : name.toLowerCase().endsWith('.jmp') ? [path] : [];
  });
}

async function main(): Promise<void> {
  const args = process.argv.slice(2);
  const goldenAt = args.indexOf('--golden');
  const goldenDir = goldenAt >= 0 ? args[goldenAt + 1] : undefined;
  const dir = args.find((a, i) => !a.startsWith('--') && !(goldenAt >= 0 && i === goldenAt + 1));
  if (!dir) {
    console.error('usage: npm run corpus -- <folder> [--golden <folder>]');
    process.exit(2);
  }

  const files = walk(dir);
  let read = 0; let partial = 0; let failed = 0; let mismatched = 0;
  const failures = new Map<string, number>();
  const unsupported = new Map<string, number>();
  const platforms = new Map<string, {converted: number, not: number}>();
  const notConverted = new Map<string, number>();
  let scripts = 0; let unreadScripts = 0; let flaggedRowStates = 0;
  const flagValues = new Map<number, number>();
  for (const file of files) {
    try {
      const table = await readJmp(new Uint8Array(readFileSync(file)));
      const bad = table.columns.filter((c) => c.decodeError);
      if (bad.length > 0) {
        partial++;
        for (const c of bad)
          unsupported.set(c.decodeError!, (unsupported.get(c.decodeError!) ?? 0) + 1);
        console.log(`PARTIAL ${file}: ${bad.map((c) => c.info.name).join(', ')}`);
      } else
        read++;
      const analyses = convertScripts(table.scripts.filter((sc) => !isSourceScript(sc)),
        table.columns.filter((c) => !c.decodeError && c.kind !== 'rowstate')
          .map((c) => ({name: c.info.name, numeric: c.kind !== 'string'})));
      scripts += table.scripts.length;
      unreadScripts += table.scripts.filter((sc) => sc.error).length;
      for (const a of analyses) {
        const key = a.platform ?? '(no platform)';
        const entry = platforms.get(key) ?? {converted: 0, not: 0};
        if (a.viewers.length > 0)
          entry.converted++;
        else {
          entry.not++;
          notConverted.set(a.reason!, (notConverted.get(a.reason!) ?? 0) + 1);
        }
        platforms.set(key, entry);
      }
      if (analyses.length > 0) {
        console.log(`SCRIPTS ${file}: ` + analyses.map((a) => `${a.name} [${a.platform ?? '-'}] ` +
          (a.viewers.length > 0 ? a.viewers.map((v) => v.type).join('+') : 'not converted')).join('; '));
      }
      for (const c of table.columns.filter((col) => col.rowStateFlags?.some((f) => f !== 0))) {
        flaggedRowStates++;
        for (const f of c.rowStateFlags!)
          flagValues.set(f, (flagValues.get(f) ?? 0) + 1);
        console.log(`ROWSTATE FLAGS ${file}: column "${c.info.name}"`);
      }
      const goldenFile = goldenDir ? join(goldenDir, basename(file).replace(/\.jmp$/i, '.json')) : undefined;
      if (goldenFile && existsSync(goldenFile)) {
        const diffs = compareToGolden(table, JSON.parse(readFileSync(goldenFile, 'utf8')) as Golden);
        if (diffs.length > 0) {
          mismatched++;
          console.log(`MISMATCH ${file}\n  ${diffs.slice(0, 5).join('\n  ')}`);
        }
      }
    } catch (e) {
      failed++;
      const message = e instanceof Error ? e.message : String(e);
      failures.set(message, (failures.get(message) ?? 0) + 1);
      console.log(`FAIL ${file}: ${message}`);
    }
  }

  const pct = (n: number) => files.length === 0 ? '0' : (100 * n / files.length).toFixed(1);
  console.log(`\n${files.length} files: ${read} read (${pct(read)}%), ${partial} with undecoded columns, ` +
    `${failed} failed${goldenDir ? `, ${mismatched} differ from the reference` : ''}`);
  for (const [m, n] of [...failures].sort((a, b) => b[1] - a[1]))
    console.log(`  ${n} × ${m}`);
  for (const [m, n] of [...unsupported].sort((a, b) => b[1] - a[1]))
    console.log(`  ${n} column(s) × ${m}`);
  console.log(`${scripts} saved script(s)${unreadScripts > 0 ? `, ${unreadScripts} unreadable` : ''}`);
  const total = (e: {converted: number, not: number}): number => e.converted + e.not;
  for (const [p, {converted, not}] of [...platforms].sort((a, b) => total(b[1]) - total(a[1])))
    console.log(`  ${p}: ${converted} converted, ${not} not`);
  for (const [m, n] of [...notConverted].sort((a, b) => b[1] - a[1]).slice(0, 15))
    console.log(`  ${n} × ${m}`);
  if (flaggedRowStates > 0) {
    console.log(`${flaggedRowStates} row-state column(s) with undecoded bytes set; values (byte4 | byte7 << 8):`);
    for (const [f, n] of [...flagValues].sort((a, b) => a[0] - b[0]))
      console.log(`  0x${f.toString(16).padStart(4, '0')}: ${n} row(s)`);
  }
  process.exit(failed + mismatched > 0 ? 1 : 0);
}

void main();
