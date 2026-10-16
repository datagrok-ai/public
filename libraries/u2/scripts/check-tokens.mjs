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

/* A raw hex colour outside the token sheets bypasses the platform aliases in tokens.css: inside
   Datagrok it would not follow the platform palette. Comments are blanked (lines kept) first; a
   match must sit in a declaration value, so an `#id` selector never counts. */
const RAW_HEX_ALLOWED = ['css/tokens.css', 'css/theme-dark.css'];
const rawHex = [];
for (const f of sheets.filter((s) => !RAW_HEX_ALLOWED.includes(s))) {
  const text = readFileSync(join(root, f), 'utf8')
    .replace(/\/\*[\s\S]*?\*\//g, (c) => c.replace(/[^\n]/g, ' '));
  for (const m of text.matchAll(/:[^;{}]*?(#[0-9a-fA-F]{3,8})\b[^;{}]*[;}]/g)) {
    const at = m.index + m[0].indexOf(m[1]);
    rawHex.push(`${f}:${text.slice(0, at).split('\n').length}  ${m[1]}`);
  }
}

if (offenders.length || rawHex.length) {
  if (offenders.length)
    console.error(`${offenders.length} undefined design token(s) — fix the name or declare it in css/tokens.css:`);
  for (const o of offenders)
    console.error(`  ${o}`);
  if (rawHex.length)
    console.error(`${rawHex.length} raw hex colour(s) — use a --dg-* token from css/tokens.css:`);
  for (const o of rawHex)
    console.error(`  ${o}`);
  process.exit(1);
}
console.log(`design tokens ok (${defined.size} declared across ${sheets.length} sheets, no raw hex colours)`);
