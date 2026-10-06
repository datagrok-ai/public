/* Noun phrases → resolution plans. Pure: the compiler runs it to validate and explain a phrase,
   the runtime runs it again to build the Playwright locator. Resolution order, at every level:
     1. the whole phrase is a registered element (or alias) — the active context's names first,
        then the global ones; an override always wins, so "save button inside toolbar" can be
        pinned to one selector regardless of composition;
     2. a leading ordinal ("second", "last");
     3. composition: "<inner> in|inside|within|on|under|of <outer>" — the inner phrase resolved
        inside the outer one (recursively; "of" also names a registered part);
     4. a registered element, else a generic kind by longest matching suffix
        ("sequence column input" → kind `input`, qualifier "sequence column").
   Inside a context, generic kinds and the context's own names are looked up in the context root
   first (the runtime's job); global names keep their global meaning. */
import {ContextEntry, ElementEntry, KindEntry, kindNames, lookupElement, lookupKind, normalizePhrase} from './registry.js';

export type Ordinal = number | 'last';

export interface KindSplit {
  kind: KindEntry;
  qualifier: string;
}

export type NounPlan =
  | {type: 'entry'; entry: ElementEntry}
  /** The longest kind suffix first; "sequence column input" is also `input` qualified
   * "sequence column", so the runtime tries every split, longest kind first. */
  | {type: 'kind'; kind: KindEntry; qualifier: string; alternatives: KindSplit[]}
  | {type: 'part'; part: string; selector: string};

export interface NounRef {
  raw: string;
  plan: NounPlan;
  /** Zero-based position among the matches, or the last one. */
  ordinal?: Ordinal;
  /** The phrase this one is resolved inside. */
  scope?: NounRef;
  /** The context active where the phrase was written. */
  context?: ContextEntry;
}

export class NounError extends Error {
  constructor(readonly phrase: string, message: string) {
    super(message);
  }
}

export const SCOPE_WORDS = ['inside', 'within', 'in', 'on', 'under', 'of'];

const ORDINAL_WORDS: Record<string, Ordinal> = {
  first: 0, second: 1, third: 2, fourth: 3, fifth: 4, sixth: 5, seventh: 6, eighth: 7, ninth: 8,
  tenth: 9, last: 'last',
};

function splitOrdinal(text: string): [Ordinal | undefined, string] {
  const m = /^(\S+)\s+(.+)$/.exec(text);
  if (!m)
    return [undefined, text];
  const word = m[1];
  if (word in ORDINAL_WORDS)
    return [ORDINAL_WORDS[word], m[2]];
  const n = /^(\d+)(st|nd|rd|th)$/.exec(word);
  if (n)
    return [Number(n[1]) - 1, m[2]];
  return [undefined, text];
}

/** Space-separated tokens with double-quoted runs kept whole. */
function tokenize(text: string): string[] {
  return text.match(/"[^"]*"|\S+/g) ?? [];
}

function splitScope(text: string): {inner: string; word: string; outer: string} | undefined {
  const tokens = tokenize(text);
  for (let i = 1; i < tokens.length - 1; i++) {
    if (SCOPE_WORDS.includes(tokens[i]))
      return {inner: tokens.slice(0, i).join(' '), word: tokens[i], outer: tokens.slice(i + 1).join(' ')};
  }
  return undefined;
}

function lookup(text: string, ctx?: ContextEntry): ElementEntry | undefined {
  return ctx?.lookup(text) ?? lookupElement(text);
}

function parseSimple(raw: string, ctx?: ContextEntry): NounPlan {
  const text = normalizePhrase(raw);
  const entry = lookup(text, ctx);
  if (entry)
    return {type: 'entry', entry};
  const quoted = /^"([^"]+)"\s*(.*)$/.exec(text) ?? /^(.*?)\s*"([^"]+)"$/.exec(text)?.reverse();
  if (quoted) {
    const kind = lookupKind(quoted[2]);
    if (!kind)
      throw new NounError(raw, `"${quoted[2]}" is not a known kind of element`);
    return {type: 'kind', kind, qualifier: quoted[1], alternatives: []};
  }
  const splits: KindSplit[] = [];
  const seen = new Set<KindEntry>();
  let whole: KindSplit | undefined;
  for (const name of kindNames()) {
    const kind = lookupKind(name)!;
    if (seen.has(kind))
      continue;
    if (text === name)
      whole = {kind, qualifier: ''};
    else if (text.endsWith(' ' + name))
      splits.push({kind, qualifier: text.slice(0, -name.length).trim()});
    else
      continue;
    seen.add(kind);
  }
  // "columns input" is the kind and also an input labelled Columns: the qualified reading is tried
  // first, so a label wins whenever one matches and the whole kind stays the fallback
  if (whole)
    splits.push(whole);
  if (splits.length === 0) {
    throw new NounError(raw, `"${text}" is neither a registered element nor "<qualifier> <kind>" ` +
      `for a known kind (${sampleKinds()})`);
  }
  return {type: 'kind', ...splits[0], alternatives: splits.slice(1)};
}

function sampleKinds(): string {
  const names = kindNames().filter((n) => !n.includes(' ')).sort();
  return names.slice(0, 8).join(', ') + (names.length > 8 ? ', …' : '');
}

export function parseNoun(raw: string, ctx?: ContextEntry): NounRef {
  const text = normalizePhrase(raw);
  const whole = lookup(text, ctx);
  if (whole)
    return {raw, plan: {type: 'entry', entry: whole}, context: ctx};
  const [ordinal, rest] = splitOrdinal(text);
  const ordinalHit = ordinal !== undefined ? lookup(rest, ctx) : undefined;
  if (ordinalHit)
    return {raw, plan: {type: 'entry', entry: ordinalHit}, ordinal, context: ctx};
  const split = splitScope(rest);
  if (!split)
    return {raw, plan: parseSimple(rest, ctx), ordinal, context: ctx};
  const scope = parseNoun(split.outer, ctx);
  if (split.word === 'of') {
    const part = normalizePhrase(split.inner);
    const owner = partOwner(scope, part);
    if (owner)
      return {raw, plan: {type: 'part', part, selector: owner.selector}, ordinal, scope: owner.scope, context: ctx};
  }
  return {raw, plan: parseSimple(split.inner, ctx), ordinal, scope, context: ctx};
}

/** The scope's element or kind that defines the part — among every kind reading of the scope
 * ("function form" is `form` qualified "function" before it is the kind `function form`), the
 * one that has it becomes the scope's primary reading. */
function partOwner(scope: NounRef, part: string): {selector: string; scope: NounRef} | undefined {
  const plan = scope.plan;
  if (plan.type === 'entry') {
    const selector = plan.entry.parts?.[part];
    return selector ? {selector, scope} : undefined;
  }
  if (plan.type !== 'kind')
    return undefined;
  const splits: KindSplit[] = [{kind: plan.kind, qualifier: plan.qualifier}, ...plan.alternatives];
  const index = splits.findIndex((s) => s.kind.parts?.[part]);
  if (index < 0)
    return undefined;
  const chosen = splits[index];
  const alternatives = splits.filter((_, i) => i !== index);
  return {selector: chosen.kind.parts![part], scope: {...scope, plan: {type: 'kind', ...chosen, alternatives}}};
}

/** Whether the runtime looks the phrase up inside its context root before the whole page. */
export function contextFirst(ref: NounRef): boolean {
  return ref.context !== undefined && !ref.scope &&
    (ref.plan.type === 'kind' || (ref.plan.type === 'entry' && ref.plan.entry.context !== undefined));
}

/** One line explaining how a phrase resolves — the compiler's diagnostic and the runtime's
 * failure message. */
export function describeNoun(ref: NounRef): string {
  const plan = ref.plan;
  const split = (s: KindSplit) => `kind "${s.kind.name}"` + (s.qualifier ? ` qualified "${s.qualifier}"` : '');
  let s = plan.type === 'entry' ? `element "${plan.entry.name}"${plan.entry.context ? ` of ${plan.entry.context.name}` : ''} [${plan.entry.selector}]` :
    plan.type === 'kind' ? [split(plan), ...plan.alternatives.map(split)].join(' or ') :
      `part "${plan.part}" [${plan.selector}]`;
  if (ref.ordinal !== undefined)
    s += ref.ordinal === 'last' ? ' (last)' : ` (#${ref.ordinal + 1})`;
  if (ref.scope)
    s += ` within ${describeNoun(ref.scope)}`;
  else if (contextFirst(ref))
    s += ` (inside ${ref.context!.name} first)`;
  return s;
}

/** Whether the plan rests on the generic kind matcher anywhere — worth a note in diagnostics. */
export function usesKinds(ref: NounRef): boolean {
  return ref.plan.type === 'kind' || (ref.scope !== undefined && usesKinds(ref.scope));
}
