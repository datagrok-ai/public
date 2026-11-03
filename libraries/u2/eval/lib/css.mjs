// The page's cascade, for what the DOM shim cannot say; at-rule blocks and pseudo-element rules are left out.

const INHERITED = new Set(['color', 'pointer-events', 'visibility', 'cursor', 'font-size', 'font-weight']);
const HOVER = 'data-eval-hover';
const parsed = new Map();

function skipBlock(css, i) {
  let depth = 0;
  for (; i < css.length; i++) {
    if (css[i] === '{')
      depth++;
    else if (css[i] === '}' && --depth === 0)
      return i + 1;
  }
  return i;
}

function declarations(body) {
  const out = [];
  let depth = 0;
  let start = 0;
  for (let i = 0; i <= body.length; i++) {
    const ch = body[i];
    if (ch === '(')
      depth++;
    else if (ch === ')')
      depth--;
    if (i < body.length && (ch !== ';' || depth > 0))
      continue;
    const part = body.slice(start, i);
    start = i + 1;
    const colon = part.indexOf(':');
    if (colon < 0)
      continue;
    const prop = part.slice(0, colon).trim();
    let value = part.slice(colon + 1).trim();
    const important = /!important\s*$/.test(value);
    value = value.replace(/\s*!important\s*$/, '');
    out.push({prop, value, important});
  }
  return out;
}

function specificity(selector) {
  const s = selector.replace(/:not\(([^()]*)\)/g, ' $1 ');
  const ids = (s.match(/#[\w-]+/g) ?? []).length;
  const classes = (s.match(/\.[\w-]+|\[[^\]]*\]|:(?!not)[\w-]+/g) ?? []).length;
  const tags = (s.replace(/\[[^\]]*\]|[.#:][\w-]+/g, ' ').match(/(^|[\s>+~])[a-zA-Z][\w-]*/g) ?? []).length;
  return ids * 10000 + classes * 100 + tags;
}

/** A sheet's rules in order: {parts: [{selector, specificity}], decls: [{prop, value, important}]}. */
export function parse(css) {
  const rules = [];
  const text = css.replace(/\/\*[\s\S]*?\*\//g, '');
  let i = 0;
  while (i < text.length) {
    const open = text.indexOf('{', i);
    if (open < 0)
      break;
    const head = text.slice(i, open).trim();
    if (head.startsWith('@')) {
      i = skipBlock(text, open);
      continue;
    }
    const close = skipBlock(text, open);
    const parts = head.split(/,(?![^(]*\))/).map((p) => p.trim()).filter((p) => p && !p.includes('::'));
    if (parts.length > 0) {
      rules.push({parts: parts.map((p) => ({selector: p, specificity: specificity(p)})),
        decls: declarations(text.slice(open + 1, close - 1))});
    }
    i = close;
  }
  return rules;
}

export function rules() {
  const all = [];
  for (const style of document.head.querySelectorAll('style')) {
    const text = style.textContent;
    if (!parsed.has(text))
      parsed.set(text, parse(text));
    all.push(...parsed.get(text));
  }
  return all;
}

/** The CSS the solution's `.css` imports installed. */
export function packageCss() {
  return document.head.querySelectorAll('style[data-package-css]').map((s) => s.textContent).join('\n');
}

function matches(el, selector) {
  const s = selector.replace(/:root\b/g, 'html').replace(/:hover\b/g, `[${HOVER}]`);
  try {
    return el.matches(s);
  }
  catch (_) {
    return false;
  }
}

function camel(prop) {
  return prop.startsWith('--') ? prop : prop.replace(/-([a-z])/g, (_, c) => c.toUpperCase());
}

/** A shorthand's longhand that the cascade weighs together with it. */
const SAME = {'background-color': ['background-color', 'background'], 'flex-grow': ['flex-grow', 'flex'],
  'column-gap': ['column-gap', 'gap']};

/** The value the cascade gives `prop` on `el` itself, before inheritance and var(); undefined when nothing sets it. */
export function declared(el, prop) {
  const props = SAME[prop] ?? [prop];
  const inline = props.map((p) => el.style?.[camel(p)]).find((v) => v !== undefined && v !== '');
  let best;
  let order = 0;
  for (const rule of rules()) {
    order++;
    for (const decl of rule.decls) {
      if (!props.includes(decl.prop))
        continue;
      for (const part of rule.parts) {
        if (!matches(el, part.selector))
          continue;
        const rank = [decl.important ? 1 : 0, part.specificity, order];
        if (best === undefined || rank[0] > best.rank[0] || (rank[0] === best.rank[0] &&
          (rank[1] > best.rank[1] || (rank[1] === best.rank[1] && rank[2] >= best.rank[2]))))
          best = {value: decl.value, rank};
      }
    }
  }
  if (inline !== undefined && inline !== '' && !(best?.rank[0] === 1))
    return String(inline);
  return best?.value;
}

function cascaded(el, prop) {
  for (let node = el; node && node.nodeType === 1; node = node.parentNode) {
    const value = declared(node, prop);
    if (value !== undefined && value !== 'inherit')
      return {value, at: node};
    if (!(prop.startsWith('--') || INHERITED.has(prop) || value === 'inherit'))
      return undefined;
  }
  return undefined;
}

function substitute(value, el, depth = 0) {
  if (depth > 20)
    return value;
  return value.replace(/var\(\s*(--[\w-]+)\s*(?:,((?:[^()]|\([^()]*\))*))?\)/g, (_, name, fallback) => {
    const found = cascaded(el, name);
    const raw = found !== undefined ? found.value : fallback?.trim() ?? '';
    return substitute(raw, found?.at ?? el, depth + 1);
  });
}

/** `prop` on `el` as the page computes it, var() resolved; `hover` puts the pointer over that element. */
export function computed(el, prop, {hover} = {}) {
  const marked = [];
  for (let node = hover; node && node.nodeType === 1; node = node.parentNode) {
    node.setAttribute(HOVER, '');
    marked.push(node);
  }
  try {
    const found = cascaded(el, prop);
    return found === undefined ? undefined : substitute(found.value, el).trim();
  }
  finally {
    for (const node of marked)
      node.removeAttribute(HOVER);
  }
}

/** How much of a flex row's free space `el` takes: `flex-grow`, or the grow part of `flex`. */
export function grow(el) {
  const value = computed(el, 'flex-grow');
  if (value === undefined || value === 'none')
    return 0;
  return value === 'auto' ? 1 : Number.parseFloat(value) || 0;
}

export const px = (value) => /^(-?[\d.]+)px$/.test(value ?? '') ? Number.parseFloat(value) : undefined;

const NAMED = {white: [255, 255, 255, 1], black: [0, 0, 0, 1], transparent: [0, 0, 0, 0], none: [0, 0, 0, 0]};

function splitTop(value) {
  const parts = [];
  let depth = 0;
  let start = 0;
  for (let i = 0; i < value.length; i++) {
    if (value[i] === '(')
      depth++;
    else if (value[i] === ')')
      depth--;
    else if (value[i] === ',' && depth === 0) {
      parts.push(value.slice(start, i).trim());
      start = i + 1;
    }
  }
  parts.push(value.slice(start).trim());
  return parts;
}

/** A colour as [r, g, b, alpha] (hex, rgb(a), color-mix in srgb, a uniform gradient, a few names), or null. */
export function color(value) {
  if (value === undefined || value === null)
    return null;
  const v = value.trim().toLowerCase();
  if (v in NAMED)
    return NAMED[v];
  let m = /^#([0-9a-f]{3,8})$/.exec(v);
  if (m) {
    const h = m[1].length < 6 ? [...m[1].slice(0, 3)].map((c) => c + c).join('') : m[1];
    const alpha = h.length === 8 ? parseInt(h.slice(6), 16) / 255 : 1;
    return [...[0, 2, 4].map((i) => parseInt(h.slice(i, i + 2), 16)), alpha];
  }
  m = /^rgba?\(\s*([\d.]+)[\s,]+([\d.]+)[\s,]+([\d.]+)(?:[\s,/]+([\d.]+))?\s*\)$/.exec(v);
  if (m)
    return [+m[1], +m[2], +m[3], m[4] === undefined ? 1 : +m[4]];
  m = /^color-mix\(\s*in srgb\s*,(.*)\)$/.exec(v);
  if (m) {
    const [one, two] = splitTop(m[1]).map((p) => /^(.*?)(?:\s+([\d.]+)%)?$/.exec(p));
    const a = color(one[1]);
    const b = color(two[1]);
    const p = one[2] !== undefined ? +one[2] / 100 : two[2] !== undefined ? 1 - +two[2] / 100 : 0.5;
    return a && b ? a.map((x, i) => i < 3 ? Math.round(x * p + b[i] * (1 - p)) : x * p + b[i] * (1 - p)) : null;
  }
  m = /^linear-gradient\((.*)\)$/.exec(v);
  if (m)
    return color(splitTop(m[1]).find((p) => color(p) !== null) ?? null);
  return null;
}

/** What `el` paints over its parent: the colours of its background layers, top first. */
function layers(el, options) {
  const value = computed(el, 'background-color', options);
  return value === undefined ? [] : splitTop(value).map((p) => color(p) ?? color(p.split(/\s+(?![^(]*\))/)[0]))
    .filter((c) => c !== null);
}

/** The colour `el` shows, as [r, g, b]: its background layers over its ancestors', down to an opaque one. */
export function background(el, options) {
  const stack = [];
  for (let node = el; node && node.nodeType === 1 && !stack.some((c) => c[3] >= 1); node = node.parentNode)
    stack.push(...layers(node, options));
  let [r, g, b] = [255, 255, 255];
  for (const [cr, cg, cb, a] of stack.reverse())
    [r, g, b] = [cr * a + r * (1 - a), cg * a + g * (1 - a), cb * a + b * (1 - a)];
  return [r, g, b].map(Math.round);
}
