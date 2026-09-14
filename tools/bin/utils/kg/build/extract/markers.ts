/// What the artifact extractors share (build-plan.md WO-3c): the `~id` and `GROK-n` tokens of a text,
/// their resolution against the home documents, the ticket stubs a mention needs, and the mentions
/// themselves counted per target.
import {loadHomes, HomeSet} from '../../homes';
import {ticketId} from '../ids';
import {Emitter} from '../emitter';
import {Row} from '../normalize';
import {BuildContext} from '../registry';

/** A `~id` in prose or a comment: an optional prefix, kebab segments, an anchor that is dropped; not a `~/` path or a `~~strike~~`. */
export const ID_TOKEN = /(?<![\w~/.\\-])~((?:[A-Z][A-Za-z]{0,5}:)?[a-z][a-z0-9]*(?:-[a-z0-9]+)*(?:\/[a-z0-9]+(?:-[a-z0-9]+)*)*)(?:#[\w-]+)?/g;
/** A line that is a comment and nothing else but a `~id`: the inline participation marker of §6. A marker after code on
 * the same line is prose, and `// ~id` in a doc comment above a declaration is the ownership marker, not this one. */
export const MARKER_LINE = /^\s*\/\/\s*~((?:[A-Z][A-Za-z]{0,5}:)?[a-z][a-z0-9]*(?:-[a-z0-9]+)*(?:\/[a-z0-9]+(?:-[a-z0-9]+)*)*)(?:#[\w-]+)?\s*$/;
export const JIRA_TOKEN = /\bGROK-\d+\b/g;
/** `[#4062](https://github.com/datagrok-ai/public/issues/4062)`: the number named twice, the link decides. */
export const GITHUB_ISSUE_LINK = /\[#(\d+)\]\(https?:\/\/github\.com\/[^)]*\/issues\/\1\)/g;

export interface Resolved {
  id: string;
  type: string;
  root?: string;
}

/** The authored ids and aliases the home documents declare, so a token resolves the way `check` resolves a reference. */
export class HomeIndex {
  private byId = new Map<string, Resolved>();

  constructor(homes: HomeSet) {
    for (const home of homes.homes) {
      const entry = {id: home.id, type: home.type.name, root: home.type.root};
      this.byId.set(home.id, entry);
      for (const alias of home.aliases) this.byId.set(alias, entry);
    }
  }

  /** [raw] with or without the sigil; a bare id is a feature, `#anchor` is dropped. */
  resolve(raw: string): Resolved | undefined {
    const id = raw.trim().replace(/^~/, '').split('#')[0];
    return this.byId.get(id);
  }

  isFeature(raw: string): boolean {
    return this.resolve(raw)?.root === 'feature';
  }
}

/** The home documents, loaded once per build and shared by every extractor that runs after the first to ask. */
export function homesOf(ctx: BuildContext): HomeSet {
  return ctx.homes ??= loadHomes(ctx.system, ctx.repoRoot);
}

/** The home a `~id` an artifact names resolves to; an unresolved token is counted under `unresolved_ids` and never
 * becomes a node, whatever edge the caller was about to draw (conventions.md §6). */
export function resolveMention(emitter: Emitter, index: HomeIndex, token: string, evidence: string): Resolved | undefined {
  const target = index.resolve(token);
  if (target) return target;
  emitter.problem('unresolved_ids', `${evidence}: ~${token} resolves to no home document`);
  return undefined;
}

/** The stub a mentioned ticket needs until the process layer (WO-5) fills it. */
export function ticketStub(emitter: Emitter, id: string): void {
  const github = /^gh:public#(\d+)$/.exec(id);
  emitter.stub(id, 'ticket', id, 'annotation', {tracker: github ? 'github' : 'jira', key: github ? `#${github[1]}` : id, kind: 'unknown', state: 'open'});
}

/** `GROK-n` keys and linked GitHub issues of [text] as ticket ids, each with how often it appears. */
export function ticketTokens(text: string): Map<string, number> {
  const out = new Map<string, number>();
  for (const m of text.matchAll(JIRA_TOKEN)) out.set(m[0], (out.get(m[0]) ?? 0) + 1);
  for (const m of text.matchAll(GITHUB_ISSUE_LINK)) out.set(ticketId(`#${m[1]}`), (out.get(ticketId(`#${m[1]}`)) ?? 0) + 1);
  return out;
}

/** The `~id` a text starts with, without the sigil or anchor, or undefined. */
export function leadingId(text: string): string | undefined {
  return new RegExp(`^${ID_TOKEN.source}`).exec(text.trim())?.[1];
}

/** Every `~id` of [text] with how often it appears, keyed by the id without sigil or anchor. */
export function idTokens(text: string): Map<string, number> {
  const out = new Map<string, number>();
  for (const m of text.matchAll(ID_TOKEN)) out.set(m[1], (out.get(m[1]) ?? 0) + 1);
  return out;
}

export interface MentionSummary {
  /** Resolved `~id` targets, by id. */
  ids: Map<string, Resolved>;
  tickets: string[];
  /** `~id` tokens no home declares. */
  unresolved: string[];
}

/** The mentions of [text] from an artifact: resolved `~id` tokens and ticket keys, each with its count; an unresolved `~id`
 * is counted as a problem and listed, never stubbed. [edgeFor] may replace the mentions row of a resolved `~id` with a
 * stronger edge (a changelog bullet naming a feature is a `changes`). */
export function emitMentions(emitter: Emitter, from: string, text: string, evidence: string, index: HomeIndex, edgeFor?: (target: Resolved) => Row | undefined): MentionSummary {
  const summary: MentionSummary = {ids: new Map(), tickets: [], unresolved: []};
  for (const [token, count] of idTokens(text)) {
    const target = resolveMention(emitter, index, token, evidence);
    if (!target) {
      summary.unresolved.push(token);
      continue;
    }
    summary.ids.set(target.id, target);
    const row = edgeFor?.(target) ?? {type: 'mentions', count};
    emitter.edge({...row, from, to: target.id, derived_by: 'annotation', confidence: 1, evidence: [evidence]});
  }
  for (const [id, count] of ticketTokens(text)) {
    ticketStub(emitter, id);
    emitter.edge({type: 'mentions', from, to: id, count, derived_by: 'annotation', confidence: 1, evidence: [evidence]});
    summary.tickets.push(id);
  }
  return summary;
}

/** `ScatterPlot` -> `scatter-plot`, `MLMethods` -> `ml-methods`, `initial runs` -> `initial-runs`. */
export function kebab(segment: string): string {
  return segment.replace(/([a-z0-9])([A-Z])/g, '$1-$2').replace(/([A-Z]+)([A-Z][a-z])/g, '$1-$2').replace(/[\s_]+/g, '-').toLowerCase().replace(/-+/g, '-').replace(/^-|-$/g, '');
}
