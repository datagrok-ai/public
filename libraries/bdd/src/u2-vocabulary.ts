import {kind, type ElementGestures, type KindDef, type MatchStrategy} from './registry.js';

export type U2Family = 'input' | 'action' | 'container' | 'collection' | 'navigation' | 'display' | 'overlay' | 'composite';

/** One `kinds` entry of a manifest (CONTRACT §4 `U2KindEntry`; bdd does not import the js-api types),
 * tagged with the manifest it came from. */
export interface U2KindSpec {
  kind: string;
  noun: string;
  aliases?: string[];
  family: U2Family;
  tag?: string;
  role?: string;
  focusRole?: string;
  nameFrom: ('label' | 'aria' | 'text' | 'title' | 'placeholder')[];
  parts: Record<string, string>;
  state: string[];
  props: string[];
  gestures: string[];
  items?: {role: string; noun: string; gestures: string[]; multi?: boolean};
  partOf?: string;
  source: 'u2' | 'ems';
}

/** What only bdd knows about a noun: Dart conventions, legacy names, extra selectors. */
export interface U2Override {
  /** Replaces the `[data-u2=…]` alternatives where the root alone is too wide: a list that is not a
   * tree's inner list, the balloons inside the notification host. */
  u2Selector?: string;
  /** Appended to the `[data-u2=…]` alternatives. */
  selector?: string;
  aliases?: string[];
  match?: MatchStrategy[];
  labelSelector?: string;
  dartNames?: string[];
  editorSelector?: string;
  parts?: Record<string, string>;
  gestures?: ElementGestures;
  description?: string;
  /** Legacy part names of `partOf` kinds: {'fb-clear': 'clear button'}; the part's noun otherwise. */
  partNames?: Record<string, string>;
  /** Whether the noun is an item kind; by default it is when a manifest names it as `items.noun`. */
  items?: boolean;
}

/** Bare words that would take the suffix of existing phrases ("Run history", "console output"), never
 * a u2 alias (CONTRACT CR2). */
const GENERIC_WORDS = ['history', 'status', 'value', 'field', 'output', 'query'];

/** A selector list split at its top-level commas — not those inside `()`, `[]` or quotes. */
export function selectorList(selector: string): string[] {
  const out: string[] = [];
  let depth = 0;
  let quote = '';
  let start = 0;
  for (let i = 0; i < selector.length; i++) {
    const c = selector[i];
    if (quote) {
      if (c === '\\')
        i++;
      else if (c === quote)
        quote = '';
    }
    else if (c === '"' || c === '\'')
      quote = c;
    else if (c === '(' || c === '[')
      depth++;
    else if (c === ')' || c === ']')
      depth--;
    else if (c === ',' && depth === 0) {
      out.push(selector.slice(start, i).trim());
      start = i + 1;
    }
  }
  out.push(selector.slice(start).trim());
  return out.filter((s) => s);
}

const unique = <T>(xs: T[]): T[] => [...new Set(xs)];
const same = (a: string, b: string): boolean => a.trim().toLowerCase() === b.trim().toLowerCase();

/** One kind per noun: the kinds sharing it are alternatives of one selector, a `partOf` kind is a part
 * of its owner, and the override adds to (and for the root selector, `match`, the label and the editor,
 * replaces) what the manifests say. An override for a noun no manifest has is an error: that noun was
 * renamed or removed. */
export function kindDefsFrom(specs: U2KindSpec[], overrides: Record<string, U2Override>): {name: string; def: KindDef}[] {
  const byNoun = new Map<string, U2KindSpec[]>();
  const partsOf = new Map<string, U2KindSpec[]>();
  for (const spec of [...specs].sort((a, b) => a.kind.localeCompare(b.kind))) {
    const [map, key] = spec.partOf === undefined ? [byNoun, spec.noun] : [partsOf, spec.partOf];
    map.set(key, [...map.get(key) ?? [], spec]);
  }
  const unknown = Object.keys(overrides).filter((noun) => !byNoun.has(noun));
  if (unknown.length > 0)
    throw new Error(`u2 vocabulary: no manifest kind has the noun ${unknown.map((n) => `"${n}"`).join(', ')} — an override names it`);
  const itemNouns = new Set(specs.flatMap((s) => s.items ? [s.items.noun] : []));

  return [...byNoun].map(([noun, group]) => {
    const o = overrides[noun] ?? {};
    // the kind the noun names leads: its name sources come first, its family decides the editor
    const lead = group.find((s) => s.kind === noun.replace(/\s+/g, '-')) ?? group[0];
    const ordered = [lead, ...group.filter((s) => s !== lead)];
    const collected: Record<string, string[]> = {};
    const add = (name: string, selector: string) => void (collected[name] = unique([...collected[name] ?? [], selector]));
    for (const spec of ordered) {
      for (const [name, selector] of Object.entries(spec.parts))
        add(name, selector);
      for (const child of partsOf.get(spec.kind) ?? [])
        add(o.partNames?.[child.kind] ?? child.noun, `[data-u2="${child.kind}"]`);
    }
    const own = Object.fromEntries(Object.entries(collected).map(([name, selectors]) => [name, selectors.join(', ')]));
    const parts = {...own, ...o.parts};
    const editor = collected.editor ?? [];
    const editorSelector = o.editorSelector ?? (lead.family === 'input' && editor.length > 0 ?
      editor.flatMap(selectorList).flatMap((e) => [`${e} input`, `${e} select`, `${e} textarea`, `${e} [contenteditable]`]).join(', ') :
      undefined);
    const aliases = unique([...ordered.flatMap((s) => s.aliases ?? []).filter((a) => !GENERIC_WORDS.some((w) => same(a, w))),
      ...o.aliases ?? []]).filter((a) => !same(a, noun));
    const match: MatchStrategy[] = o.match ?? unique<MatchStrategy>([...ordered.flatMap((s) => s.nameFrom), 'name',
      ...o.dartNames ? ['dart' as const] : []]);
    const items = o.items ?? itemNouns.has(noun);
    const def: KindDef = {
      selector: [o.u2Selector ?? group.map((s) => `[data-u2="${s.kind}"]`).join(', '), o.selector].filter((s) => s).join(', '),
      ...aliases.length > 0 ? {aliases} : {},
      match,
      ...items ? {items} : {},
      ...pick({labelSelector: o.labelSelector ?? own.label ?? own.title, dartNames: o.dartNames, editorSelector,
        gestures: o.gestures, description: o.description}),
      ...Object.keys(parts).length > 0 ? {parts} : {},
    };
    return {name: noun, def};
  });
}

function pick<T extends object>(o: T): Partial<T> {
  return Object.fromEntries(Object.entries(o).filter(([, v]) => v !== undefined)) as Partial<T>;
}

/** Registers the kind of every u2 noun; bindings/common/kinds.ts is the one caller (a binding module
 * that registers nouns is imported by every generated spec). */
export function registerU2Kinds(specs: U2KindSpec[], overrides: Record<string, U2Override>): void {
  for (const {name, def} of kindDefsFrom(specs, overrides))
    kind(name, def);
}
