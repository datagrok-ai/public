/* Fails on any var(--dg-…) naming a token no sheet declares. A misspelled token silently drops
   the whole declaration, so without this nothing in the toolchain notices the dead rule. */
import {existsSync, readdirSync, readFileSync} from 'node:fs';
import {join, resolve} from 'node:path';
import {fileURLToPath} from 'node:url';

const root = fileURLToPath(new URL('..', import.meta.url));
const PLATFORM_CSS = join(root, '..', '..', '..', 'core', 'client', 'xamgle', 'web', 'datagrok.css');

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

/** The row and header states a table paints, each to read apart from every other: a hovered row that looks like a
 * group header is B-6. A state with `over` is a tint laid over that state. */
export const STATES = {
  base: {token: '--u2-row-base-bg'},
  hover: {token: '--u2-row-hover-bg', over: 'base'},
  current: {token: '--u2-row-current-bg', over: 'base'},
  selected: {token: '--u2-row-selected-bg', over: 'base'},
  header: {token: '--u2-section-header-bg', over: 'base'},
};

/** CIE76 ΔE of a just-noticeable difference. */
export const MIN_DELTA_E = 2.3;

function declarations(css, selector) {
  const out = new Map();
  const text = css.replace(/\/\*[\s\S]*?\*\//g, '');
  for (const block of text.matchAll(/([^{}]+)\{([^{}]*)\}/g)) {
    if (block[1].split(';').pop().trim() !== selector)
      continue;
    for (const d of block[2].matchAll(/(--[\w-]+)\s*:\s*([^;]+);/g))
      out.set(d[1], d[2].trim());
  }
  return out;
}

/** `value` with every var() replaced by what the first layer that declares it says, else by its fallback. */
function expand(value, layers, depth = 0) {
  const at = value.indexOf('var(');
  if (at < 0)
    return value;
  if (depth > 32)
    throw new Error(`var() nests too deep in ${value}`);
  let end = at + 4;
  let comma = -1;
  for (let open = 1; open > 0; end++) {
    if (end >= value.length)
      throw new Error(`unclosed var() in ${value}`);
    open += value[end] === '(' ? 1 : value[end] === ')' ? -1 : 0;
    if (value[end] === ',' && open === 1 && comma < 0)
      comma = end;
  }
  const name = value.slice(at + 4, comma < 0 ? end - 1 : comma).trim();
  const own = layers.find((l) => l.has(name))?.get(name);
  const fallback = comma < 0 ? undefined : value.slice(comma + 1, end - 1).trim();
  if (own === undefined && fallback === undefined)
    throw new Error(`${name} is not declared`);
  const inner = expand(own ?? fallback, layers, depth + 1);
  return expand(value.slice(0, at) + inner + value.slice(end), layers, depth);
}

const NAMED = {transparent: [0, 0, 0, 0], white: [255, 255, 255, 1], black: [0, 0, 0, 1]};

function rgba(text) {
  const v = text.trim().toLowerCase();
  if (NAMED[v])
    return NAMED[v];
  const hex = v.match(/^#([0-9a-f]{3,8})$/);
  if (hex) {
    const h = hex[1].length <= 4 ? [...hex[1]].map((c) => c + c).join('') : hex[1];
    const n = (i) => parseInt(h.slice(i, i + 2), 16);
    return [n(0), n(2), n(4), h.length === 8 ? n(6) / 255 : 1];
  }
  const fn = v.match(/^rgba?\(([^)]*)\)$/);
  if (fn) {
    const [r, g, b, a = '1'] = fn[1].split(/[\s,/]+/).filter((x) => x);
    return [+r, +g, +b, a.endsWith('%') ? parseFloat(a) / 100 : +a];
  }
  const mix = v.match(/^color-mix\(\s*in srgb\s*,(.*)\)$/);
  if (mix) {
    let depth = 0;
    const parts = [''];
    for (const c of mix[1]) {
      depth += c === '(' ? 1 : c === ')' ? -1 : 0;
      if (c === ',' && depth === 0)
        parts.push('');
      else
        parts[parts.length - 1] += c;
    }
    const [[a, pa], [b, pb]] = parts.map((p) => {
      const m = p.trim().match(/^(.*?)(?:\s+([\d.]+)%)?$/);
      return [rgba(m[1]), m[2] === undefined ? undefined : +m[2] / 100];
    });
    const p = pa ?? (pb === undefined ? 0.5 : 1 - pb);
    return a.map((x, i) => x * p + b[i] * (1 - p));
  }
  throw new Error(`cannot read "${text}" as a colour`);
}

const over = ([r, g, b, a], [R, G, B]) => [r * a + R * (1 - a), g * a + G * (1 - a), b * a + B * (1 - a), 1];

function lab([r, g, b]) {
  const lin = (c) => (c /= 255) <= 0.04045 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4;
  const [R, G, B] = [lin(r), lin(g), lin(b)];
  const f = (t) => t > 216 / 24389 ? Math.cbrt(t) : (24389 / 27 * t + 16) / 116;
  const x = f((0.4124 * R + 0.3576 * G + 0.1805 * B) / 0.95047);
  const y = f(0.2126 * R + 0.7152 * G + 0.0722 * B);
  const z = f((0.0193 * R + 0.1192 * G + 0.9505 * B) / 1.08883);
  return [116 * y - 16, 500 * (x - y), 200 * (y - z)];
}

function deltaE(a, b) {
  const [p, q] = [lab(a), lab(b)];
  return Math.hypot(p[0] - q[0], p[1] - q[1], p[2] - q[2]);
}
const rgb = (c) => `rgb(${c.slice(0, 3).map((v) => Math.round(v)).join(', ')})`;

/** Every pair of {@link STATES} resolved as the standalone sheets set them, under theme-dark.css, and with the
 * platform's palette when core/ sits next to public/: the pairs that do not differ by {@link MIN_DELTA_E}, and the
 * closest pair. */
export function checkStates(read = (f) => readFileSync(join(root, f), 'utf8'),
  platform = existsSync(PLATFORM_CSS) ? readFileSync(PLATFORM_CSS, 'utf8') : undefined) {
  const u2 = declarations(read('css/tokens.css'), ':root');
  const scenarios = [['standalone', [u2]], ['dark', [declarations(read('css/theme-dark.css'), '[data-dg-theme=\'dark\']'), u2]]];
  if (platform !== undefined)
    scenarios.push(['in Datagrok', [u2, declarations(platform, ':root')]]);
  const problems = Object.entries(STATES).filter(([, {token}]) => !u2.has(token))
    .map(([state, {token}]) => `${token} (the ${state} state) is not declared in css/tokens.css`);
  let closest = null;
  if (problems.length)
    return {problems, closest};
  for (const [scenario, layers] of scenarios) {
    const painted = {};
    for (const [state, {token, over: under}] of Object.entries(STATES)) {
      try {
        const colour = rgba(expand(`var(${token})`, layers));
        painted[state] = under === undefined || painted[under] === undefined ? over(colour, [255, 255, 255]) :
          over(colour, painted[under]);
      }
      catch (e) {
        problems.push(`${scenario}: ${token} (the ${state} state): ${e.message}`);
      }
    }
    const states = Object.keys(painted);
    for (let i = 0; i < states.length; i++) {
      for (let j = i + 1; j < states.length; j++) {
        const [a, b] = [states[i], states[j]];
        const d = deltaE(painted[a], painted[b]);
        if (closest === null || d < closest.d)
          closest = {d, a, b, scenario};
        if (d < MIN_DELTA_E) {
          problems.push(`${scenario}: the ${a} row (${STATES[a].token}, ${rgb(painted[a])}) and the ${b} ` +
            `(${STATES[b].token}, ${rgb(painted[b])}) differ by ΔE ${d.toFixed(2)} < ${MIN_DELTA_E}: a ${a} row reads ` +
            `as a ${b}`);
        }
      }
    }
  }
  return {problems, closest};
}

if (resolve(process.argv[1] ?? '') === fileURLToPath(import.meta.url)) {
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

  const states = checkStates();
  if (offenders.length || rawHex.length || states.problems.length) {
    if (offenders.length)
      console.error(`${offenders.length} undefined design token(s) — fix the name or declare it in css/tokens.css:`);
    for (const o of offenders)
      console.error(`  ${o}`);
    if (rawHex.length)
      console.error(`${rawHex.length} raw hex colour(s) — use a --dg-* token from css/tokens.css:`);
    for (const o of rawHex)
      console.error(`  ${o}`);
    if (states.problems.length)
      console.error(`${states.problems.length} row state problem(s) — each state needs its own colour in css/tokens.css:`);
    for (const o of states.problems)
      console.error(`  ${o}`);
    process.exit(1);
  }
  const {d, a, b, scenario} = states.closest;
  console.log(`design tokens ok (${defined.size} declared across ${sheets.length} sheets, no raw hex colours; ` +
    `row states apart, closest ${a}/${b} ΔE ${d.toFixed(2)} ${scenario})`);
}
