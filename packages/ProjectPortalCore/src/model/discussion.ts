/** Project discussion and Ask SME (F-51), on the platform's own chats. Pure: the messages of a
 * program's chat in, threads, feed lines and the question's text out.
 *
 * The platform keeps one chat per entity and has no threads, anchors or recipients. The portal adds
 * them in the message itself: a first line naming what the post is about and who it asks, with
 * tokens such as `pp:criterion:<id>` that the portal reads back. The platform's own Chats panel
 * shows the same messages, header line included. */

/** What a thread is about: a criterion of the program's gates, a saved run, a document. */
export type AnchorKind = 'criterion' | 'run' | 'document';
export const ANCHOR_KINDS: AnchorKind[] = ['criterion', 'run', 'document'];
export const ANCHOR_LABELS: Record<AnchorKind, string> = {criterion: 'Criterion', run: 'Run', document: 'Document'};

export interface Anchor {
  kind: AnchorKind;
  /** The record's id: gate_criterion, model_run, document. */
  id: string;
  /** Its name as it was when the thread started. */
  label: string;
}

/** One message of the program's chat, as the platform returns it. */
export interface ChatMessage {
  id: string;
  chatId: string;
  /** The author's user id. */
  by: string | null;
  byName?: string | null;
  /** ISO time. */
  on: string;
  text: string;
}

/** A message read back: what the header said, and the body without it. */
export interface Post extends ChatMessage {
  anchor: Anchor | null;
  /** The users a question was put to; empty for a plain post. */
  askedOf: string[];
  /** Whether it is an Ask SME question. */
  question: boolean;
  body: string;
}

export interface Thread {
  /** `general`, or the anchor's `kind:id`. */
  key: string;
  anchor: Anchor | null;
  /** Oldest first, as a conversation reads. */
  posts: Post[];
  /** The newest post's time. */
  last: string;
}

/** Where the portfolio sends Ask SME questions. The question always goes into the program's
 * discussion; `channel` picks who is told: the portfolio's approvers (`discussion`, the default),
 * the members of `group`, or the `users` listed. `external` is a placeholder: nothing is sent
 * outside the platform; the question stays in the discussion and the approvers are told. */
export type SmeChannel = 'discussion' | 'group' | 'user' | 'external';
export const SME_CHANNELS: SmeChannel[] = ['discussion', 'group', 'user', 'external'];

export interface AskSmeConfig {
  channel?: SmeChannel;
  /** A group's name, for `group`. */
  group?: string;
  /** Logins, for `user`. */
  users?: string[];
  /** For `external`: where it would go (email, teams, slack). Not wired. */
  external?: {kind?: string, address?: string};
  /** A line the dialog shows above the question, e.g. who answers and how fast. */
  note?: string;
}

export const GENERAL = 'general';
export const anchorKey = (a: Anchor | null): string => a ? `${a.kind}:${a.id}` : GENERAL;
const ID = /^[A-Za-z0-9-]{1,64}$/;
/** A label goes into the header line: one line, no quote that ends it early, short. */
const cleanLabel = (s: string): string => s.replace(/\s+/g, ' ').replace(/"/g, '\'').trim().slice(0, 120);

/** What a post's header says about its anchor, e.g. `Re: criterion "Potency"`. */
export const describeAnchor = (a: Anchor): string => `${ANCHOR_LABELS[a.kind].toLowerCase()} "${a.label}"`;

/** The text sent to the chat: a header line when the post is about something or asks someone,
 * then the body. [names] are the asked users' display names, shown to people. */
export function composePost(body: string, opts: {anchor?: Anchor | null, askedOf?: string[], names?: string[],
  question?: boolean} = {}): string {
  const text = body.trim();
  const parts: string[] = [];
  const tokens: string[] = [];
  if (opts.question) {
    parts.push(`Question for ${opts.names?.length ? opts.names.join(', ') : 'an expert'}`);
    tokens.push(`pp:sme:${(opts.askedOf ?? []).filter((x) => ID.test(x)).join(',')}`);
  }
  if (opts.anchor) {
    if (!ID.test(opts.anchor.id))
      throw new Error(`"${opts.anchor.id}" is not a record id`);
    parts.push(`${opts.question ? 'about' : 'Re:'} ${describeAnchor({...opts.anchor, label: cleanLabel(opts.anchor.label)})}`);
    tokens.push(`pp:${opts.anchor.kind}:${opts.anchor.id}`);
  }
  return parts.length ? `${parts.join(' ')} · ${tokens.join(' ')}\n${text}` : text;
}

const TOKEN = /\bpp:(criterion|run|document|sme):([A-Za-z0-9,-]*)/g;

/** Reads a message back. A message without the portal's header (one written in the platform's
 * own Chats panel) is a plain post in the general thread. */
export function parsePost(m: ChatMessage): Post {
  const nl = m.text.indexOf('\n');
  const first = nl < 0 ? m.text : m.text.slice(0, nl);
  const tokens = [...first.matchAll(TOKEN)];
  if (!tokens.length || !/ · pp:/.test(first))
    return {...m, anchor: null, askedOf: [], question: false, body: m.text.trim()};
  let anchor: Anchor | null = null;
  let askedOf: string[] = [];
  let question = false;
  for (const [, kind, value] of tokens) {
    if (kind === 'sme') {
      question = true;
      askedOf = value.split(',').filter((x) => x);
    } else if (!anchor && value)
      anchor = {kind: kind as AnchorKind, id: value, label: labelIn(first, kind as AnchorKind) ?? value};
  }
  return {...m, anchor, askedOf, question, body: nl < 0 ? '' : m.text.slice(nl + 1).trim()};
}

/** The label a header gives its anchor: `criterion "Potency"` → Potency. */
function labelIn(header: string, kind: AnchorKind): string | null {
  const m = header.match(new RegExp(`${ANCHOR_LABELS[kind].toLowerCase()} "([^"]*)"`));
  return m ? m[1] : null;
}

/** Posts grouped by what they are about: newest thread first, each oldest post first. */
export function threadsOf(messages: ChatMessage[]): Thread[] {
  const byKey = new Map<string, Thread>();
  for (const p of messages.map(parsePost).sort((a, b) => a.on.localeCompare(b.on))) {
    const key = anchorKey(p.anchor);
    const t = byKey.get(key) ?? {key, anchor: p.anchor, posts: [], last: p.on};
    t.posts.push(p);
    t.last = p.on;
    byKey.set(key, t);
  }
  return [...byKey.values()].sort((a, b) => b.last.localeCompare(a.last));
}

/** A thread's title: what it is about. */
export const threadTitle = (t: Pick<Thread, 'anchor'>): string => t.anchor ?
  `${ANCHOR_LABELS[t.anchor.kind]}: ${t.anchor.label}` : 'General';

/** The address of a thread: the platform's page of the chat, with the thread in the fragment. A
 * decision cites a thread by this address (a link in its evidence). */
export const threadUrl = (origin: string, chatId: string, key: string): string =>
  `${origin.replace(/\/$/, '')}/forums/${chatId}#pp-${key.replace(':', '-')}`;

/** The decisions whose evidence links to [url]: the thread shows what it led to. */
export function decisionsCiting<T extends {target: string}>(url: string, evidence: Map<string, T[]>): string[] {
  return [...evidence.entries()].filter(([, items]) => items.some((e) => e.target === url)).map(([id]) => id);
}

/** Checks a portfolio's askSme block. */
export function validateAskSme(c: AskSmeConfig | null | undefined): string[] {
  if (!c)
    return [];
  const errors: string[] = [];
  const channel = c.channel ?? 'discussion';
  if (!SME_CHANNELS.includes(channel))
    errors.push(`askSme: channel "${channel}" is not one of ${SME_CHANNELS.join(', ')}`);
  if (channel === 'group' && !c.group?.trim())
    errors.push('askSme: channel group needs a group');
  if (channel === 'user' && !(c.users ?? []).some((u) => u.trim()))
    errors.push('askSme: channel user needs users');
  return errors;
}

/** Who a question goes to, in words, and what is resolved where. */
export interface SmeTarget {
  channel: SmeChannel;
  /** Groups whose members are told: the portfolio's approvers group id, or a group's name. */
  groups: string[];
  /** Logins told. */
  logins: string[];
  /** Whether the portfolio's approvers (and the program's approver roles) are told. */
  approvers: boolean;
  /** Said in the dialog: an external channel that is not wired. */
  notice: string | null;
}

/** Resolves the portfolio's config into whom to tell. */
export function smeTargetOf(c: AskSmeConfig | null | undefined): SmeTarget {
  const channel = c?.channel ?? 'discussion';
  if (channel === 'group' && c?.group)
    return {channel, groups: [c.group], logins: [], approvers: false, notice: null};
  if (channel === 'user' && c?.users?.length)
    return {channel, groups: [], logins: c.users.map((u) => u.trim()).filter((u) => u), approvers: false, notice: null};
  if (channel === 'external') {
    return {channel, groups: [], logins: [], approvers: true, notice: `The ${c?.external?.kind ?? 'external'} channel is not ` +
      'connected: nothing leaves the platform. The question goes into the discussion and the approvers are told.'};
  }
  return {channel: 'discussion', groups: [], logins: [], approvers: true, notice: null};
}

/** What the question carries so the expert needs no long explanation. */
export interface SmeContext {
  portfolio: {code: string, name: string};
  program: {code: string, name: string};
  stage: string | null;
  /** The next gate and its readiness, e.g. `Lead optimization: 2 of 5 required met or waived (40%)`. */
  nextGate: string | null;
  /** The open required criteria, one line each. */
  gaps?: string[];
  anchor?: Anchor | null;
  /** Where the anchored item opens. */
  anchorUrl?: string | null;
  /** The program's page in the app. */
  appUrl: string;
}

/** The context from where the program stands: its stage, and the next gate's readiness. */
export function smeContextFrom(x: {portfolio: {code: string, name: string}, program: {code: string, name: string,
  stage_id?: string | null}, stages: {id: string, name: string}[], readiness: {stageId: string | null, empty: boolean,
  requiredMet: number, requiredTotal: number, score: number | null, gaps: string[]} | null, anchor?: Anchor | null,
  anchorUrl?: string | null, appUrl: string}): SmeContext {
  const name = (id?: string | null) => x.stages.find((s) => s.id === id)?.name ?? null;
  const r = x.readiness;
  const gate = r && r.stageId ? `${name(r.stageId)}: ${r.empty ? 'no required criteria defined' :
    `${r.requiredMet} of ${r.requiredTotal} required met or waived (${r.score === null ? 'no criteria' : `${Math.round(r.score * 100)}%`})`}` : null;
  return {portfolio: {code: x.portfolio.code, name: x.portfolio.name}, program: {code: x.program.code, name: x.program.name},
    stage: name(x.program.stage_id), nextGate: gate, gaps: r?.gaps ?? [], anchor: x.anchor ?? null, anchorUrl: x.anchorUrl ?? null,
    appUrl: x.appUrl};
}

/** The context block appended to a question. */
export function smeContextLines(c: SmeContext): string[] {
  return [
    `Portfolio: ${c.portfolio.name} (${c.portfolio.code})`,
    `Program: ${c.program.code}${c.program.name && c.program.name !== c.program.code ? ` — ${c.program.name}` : ''}`,
    `Stage: ${c.stage ?? '—'}`,
    ...(c.nextGate ? [`Next gate: ${c.nextGate}`] : []),
    ...(c.gaps?.length ? [`Open: ${c.gaps.slice(0, 3).join('; ')}${c.gaps.length > 3 ? `; ${c.gaps.length - 3} more` : ''}`] : []),
    ...(c.anchor ? [`About: ${describeAnchor(c.anchor)}${c.anchorUrl ? ` ${c.anchorUrl}` : ''}`] : []),
    `Open in the portal: ${c.appUrl}`,
  ];
}

/** The question's body: what was asked, then the context. */
export function smeQuestionBody(question: string, c: SmeContext): string {
  if (!question.trim())
    throw new Error('Type the question');
  return [question.trim(), '', 'Context:', ...smeContextLines(c)].join('\n');
}

/** A one-line notification for an asked user. */
export const smeNotice = (asker: string, c: Pick<SmeContext, 'program' | 'portfolio' | 'anchor'>, question: string): string =>
  `${asker} asks you about ${c.program.code} (${c.portfolio.code})${c.anchor ? `, ${describeAnchor(c.anchor)}` : ''}: ` +
  `${question.trim().split('\n')[0].slice(0, 160)}`;

/** A feed line for a post. */
export function postFeedText(p: Post): string {
  const about = p.anchor ? ` on ${describeAnchor(p.anchor)}` : '';
  const first = p.body.split('\n')[0];
  const short = first.length > 100 ? `${first.slice(0, 99)}…` : first;
  return p.question ? `Question to an expert${about}: ${short}` : `Discussion${about}: ${short}`;
}

/** An example thread to seed. */
export interface ExampleThread {
  portfolio: string;
  program: string;
  /** The anchor: a criterion by stage and code. */
  criterion?: {stage: string, code: string};
  question?: boolean;
  text: string;
}

/** The example posts not yet in the chat: a post is there when a message's body is its text. */
export function missingExamplePosts<T extends {text: string}>(examples: T[], bodies: string[]): T[] {
  const have = new Set(bodies.map((b) => b.split('\n\nContext:')[0].trim()));
  return examples.filter((e) => !have.has(e.text.trim()));
}
