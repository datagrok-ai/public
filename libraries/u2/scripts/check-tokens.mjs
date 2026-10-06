/* Fails on any var(--dg-…) naming a token no sheet declares. A misspelled token silently drops
   the whole declaration, so without this nothing in the toolchain notices the dead rule. */
import {readdirSync, readFileSync} from 'node:fs';
import {join} from 'node:path';
import {fileURLToPath} from 'node:url';

const root = fileURLToPath(new URL('..', import.meta.url));

function files(dir, ext) {
  const out = [];
  for (const e of readdirSync(join(root, dir), {withFileTypes: true})) {
    const rel = `${dir}/${e.name}`;
    if (e.isDirectory())
      out.push(...files(rel, ext));
    else if (e.name.endsWith(ext) && !e.name.endsWith('.d.ts'))
      out.push(rel);
  }
  return out;
}

const sheets = files('css', '.css');
const defined = new Set();
for (const f of sheets)
  for (const m of readFileSync(join(root, f), 'utf8').matchAll(/(--dg-[\w-]+)\s*:/g))
    defined.add(m[1]);

const offenders = [];
for (const f of [...sheets, ...files('src', '.ts')]) {
  const lines = readFileSync(join(root, f), 'utf8').split(/\r?\n/);
  for (let i = 0; i < lines.length; i++)
    for (const m of lines[i].matchAll(/var\(\s*(--dg-[\w-]+)/g))
      if (!defined.has(m[1]))
        offenders.push(`${f}:${i + 1}  ${m[1]}`);
}

if (offenders.length) {
  console.error(`${offenders.length} undefined design token(s) — fix the name or declare it in css/tokens.css:`);
  for (const o of offenders)
    console.error(`  ${o}`);
  process.exit(1);
}
console.log(`design tokens ok (${defined.size} declared across ${sheets.length} sheets)`);
