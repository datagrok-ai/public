/// The Datagrok function annotation grammar (`//name:`, `#input:`, `--connection:`) as `grok kg build` reads it
/// (build-plan.md WO-3a): check.ts `getFuncMetadata` ported with every key kept verbatim, class-method starts
/// accepted (`detectMolecules(col) {` in detectors.js) and a comment prefix per language.
import {headerTags} from '../../utils';

export interface HeaderParam {
  type: string;
  name: string;
  default?: string;
  /** The `{ semType: Molecule; optional: true }` options, values verbatim. */
  options: Record<string, string>;
  /** The trailing `[description]`. */
  description?: string;
}

export interface Header {
  /** 1-based line of the first header line. */
  line: number;
  /** Every key verbatim in file order; a repeated key (input, output, test) keeps every value. */
  keys: Record<string, string[]>;
  name?: string;
  description?: string;
  inputs: HeaderParam[];
  outputs: HeaderParam[];
  tags: string[];
  /** `meta.*` keys without the prefix, verbatim. */
  meta: Record<string, string>;
}

export interface Declaration {
  line: number;
  name: string;
  kind: 'function' | 'arrow' | 'method';
  /** The enclosing class of a method. */
  owner?: string;
}

export interface HeaderBlock {
  header: Header;
  declaration: Declaration;
  /** The lines after the declaration up to the next header block, for the detector return scan. */
  body: string;
}

const COMMENT_PREFIX: Record<string, string> = {ts: '//', js: '//', python: '#', r: '#', julia: '#', octave: '#', grok: '#', sql: '--'};
const KNOWN_TAGS = new Set([...headerTags, 'feature']);
const FUNCTION_START = /^\s*(?:export\s+)?(?:default\s+)?(?:async\s+)?function\s*\*?\s*([\w$]+)\s*\(/;
const ARROW_START = /^\s*(?:export\s+)?(?:const|let|var)\s+([\w$]+)\s*(?::[^=]+)?=\s*(?:async\s+)?(?:\([^)]*\)|[\w$]+)\s*(?::\s*[^=]+)?=>/;
const METHOD_START = /^\s*(?:(?:public|private|protected|static|async|override)\s+)*([\w$]+)\s*\([^)]*\)?\s*(?::\s*[^{=;]+)?\s*\{?\s*$/;
const CLASS_START = /^\s*(?:export\s+)?(?:default\s+)?(?:abstract\s+)?class\s+([\w$]+)/;
const NOT_A_METHOD = new Set(['if', 'for', 'while', 'switch', 'catch', 'return', 'constructor', 'function', 'typeof', 'await', 'new', 'else', 'do', 'with', 'super', 'this', 'throw', 'yield', 'delete', 'void']);
const lineRegexes = new Map<string, RegExp>();

export function commentPrefix(language: string): string {
  return COMMENT_PREFIX[language] ?? '//';
}

/** `//key: value` for one prefix; `////` and `## ...` lines are not header lines, as in utils.fileParamRegex. */
function lineRegex(prefix: string): RegExp {
  let re = lineRegexes.get(prefix);
  if (!re) {
    const p = prefix.replace(/[.*+?^${}()|[\]\\/]/g, '\\$&');
    lineRegexes.set(prefix, re = new RegExp(`^${p}\\s*(?!\\s*${p})([^:\\s][^:]*?)\\s*:\\s*(.*)$`));
  }
  return re;
}

/** The header the comment [lines] spell, or null when none of them carries a known tag or a `meta.*` key. */
export function parseHeaderLines(lines: string[], firstLine: number, prefix = '//'): Header | null {
  const re = lineRegex(prefix);
  const header: Header = {line: firstLine, keys: {}, inputs: [], outputs: [], tags: [], meta: {}};
  let known = false;
  for (const raw of lines) {
    const m = re.exec(raw.trim());
    if (!m) continue;
    const key = m[1];
    const value = m[2].trim();
    (header.keys[key] ??= []).push(value);
    if (key.startsWith('meta.')) {
      header.meta[key.slice(5)] = value;
      known = true;
      continue;
    }
    if (!KNOWN_TAGS.has(key)) continue;
    known = true;
    switch (key) {
      case 'name': header.name = value.split(/\s[[{]/)[0].trim(); break;
      case 'description': header.description = value; break;
      case 'input': header.inputs.push(parseParam(value)); break;
      case 'output': header.outputs.push(parseParam(value)); break;
      case 'tags': header.tags.push(...value.split(',').map((t) => t.trim()).filter(Boolean)); break;
    }
  }
  return known ? header : null;
}

/** `type name[:sub] [= default] [{ options }] [[description]]`; options split on `;`, or on top-level `,` when there is no `;`. */
export function parseParam(text: string): HeaderParam {
  const s = text.trim();
  const head = /^(\S+)(?:\s+([^\s=[{]+))?/.exec(s);
  const param: HeaderParam = {type: head?.[1] ?? '', name: head?.[2] ?? '', options: {}};
  let rest = s.slice(head?.[0].length ?? 0);
  const desc = /\[([^\]]*)\]\s*$/.exec(rest);
  if (desc && !/^\s*=\s*\[/.test(rest)) {
    param.description = desc[1].trim();
    rest = rest.slice(0, desc.index);
  }
  const open = topLevelIndex(rest, '{');
  if (open >= 0) {
    const close = matching(rest, open);
    for (const piece of splitTopLevel(rest.slice(open + 1, close), ';', ',')) {
      const colon = piece.indexOf(':');
      if (colon > 0) param.options[piece.slice(0, colon).trim()] = piece.slice(colon + 1).trim();
    }
    rest = rest.slice(0, open) + rest.slice(close + 1);
  }
  const eq = rest.indexOf('=');
  if (eq >= 0 && rest.slice(eq + 1).trim()) param.default = rest.slice(eq + 1).trim();
  return param;
}

/** Header blocks directly above function, arrow and class-method starts in a TS/JS source, each with what follows it. */
export function parseFunctionHeaders(text: string, prefix = '//'): HeaderBlock[] {
  const lines = text.split(/\r?\n/);
  const blocks: HeaderBlock[] = [];
  let owner: string | undefined;
  for (let i = 0; i < lines.length; i++) {
    const cls = CLASS_START.exec(lines[i]);
    if (cls) owner = cls[1];
    const declaration = matchStart(lines[i], i + 1, owner);
    if (!declaration) continue;
    let j = i - 1;
    while (j >= 0 && lines[j].trim().startsWith(prefix)) j--;
    const header = parseHeaderLines(lines.slice(j + 1, i), j + 2, prefix);
    if (!header) continue;
    if (blocks.length) blocks[blocks.length - 1].body = lines.slice(blocks[blocks.length - 1].declaration.line, j + 1).join('\n');
    blocks.push({header, declaration, body: lines.slice(i + 1).join('\n')});
  }
  return blocks;
}

/** The header at the top of a script file (leading blank lines skipped), or null. */
export function parseScriptHeader(text: string, language: string): Header | null {
  const prefix = commentPrefix(language);
  const lines = text.split(/\r?\n/);
  let i = 0;
  while (i < lines.length && !lines[i].trim()) i++;
  const first = i;
  while (i < lines.length && lines[i].trim().startsWith(prefix)) i++;
  return parseHeaderLines(lines.slice(first, i), first + 1, prefix);
}

/** Every `--` block of a query file that names a query; a file may hold several, separated by `--end`. */
export function parseQueryHeaders(text: string): Header[] {
  const lines = text.split(/\r?\n/);
  const headers: Header[] = [];
  for (let i = 0; i < lines.length;) {
    if (!lines[i].trim().startsWith('--')) {
      i++;
      continue;
    }
    const start = i;
    while (i < lines.length && lines[i].trim().startsWith('--')) i++;
    const header = parseHeaderLines(lines.slice(start, i), start + 1, '--');
    if (header?.name) headers.push(header);
  }
  return headers;
}

function matchStart(line: string, lineNo: number, owner: string | undefined): Declaration | undefined {
  const fn = FUNCTION_START.exec(line);
  if (fn) return {line: lineNo, name: fn[1], kind: 'function'};
  const arrow = ARROW_START.exec(line);
  if (arrow) return {line: lineNo, name: arrow[1], kind: 'arrow'};
  const method = METHOD_START.exec(line);
  if (method && !NOT_A_METHOD.has(method[1])) return {line: lineNo, name: method[1], kind: 'method', owner};
  return undefined;
}

function topLevelIndex(s: string, ch: string): number {
  let depth = 0;
  let quote = '';
  for (let i = 0; i < s.length; i++) {
    const c = s[i];
    if (quote) {
      if (c === quote) quote = '';
      continue;
    }
    if (c === '\'' || c === '"' || c === '`') quote = c;
    else if (c === ch && depth === 0) return i;
    else if ('([{'.includes(c)) depth++;
    else if (')]}'.includes(c)) depth--;
  }
  return -1;
}

/** Index of the bracket closing the one at [open]; the end of the string when unbalanced. */
function matching(s: string, open: number): number {
  let depth = 0;
  let quote = '';
  for (let i = open; i < s.length; i++) {
    const c = s[i];
    if (quote) {
      if (c === quote) quote = '';
      continue;
    }
    if (c === '\'' || c === '"' || c === '`') quote = c;
    else if ('([{'.includes(c)) depth++;
    else if (')]}'.includes(c) && --depth === 0) return i;
  }
  return s.length;
}

function splitTopLevel(s: string, separator: string, fallback: string): string[] {
  const split = (sep: string) => {
    const parts: string[] = [];
    let depth = 0;
    let quote = '';
    let at = 0;
    for (let i = 0; i < s.length; i++) {
      const c = s[i];
      if (quote) {
        if (c === quote) quote = '';
        continue;
      }
      if (c === '\'' || c === '"' || c === '`') quote = c;
      else if ('([{'.includes(c)) depth++;
      else if (')]}'.includes(c)) depth--;
      else if (c === sep && depth === 0) {
        parts.push(s.slice(at, i));
        at = i + 1;
      }
    }
    parts.push(s.slice(at));
    return parts.map((p) => p.trim()).filter(Boolean);
  };
  const bySeparator = split(separator);
  return bySeparator.length > 1 || s.includes(separator) ? bySeparator : split(fallback);
}
