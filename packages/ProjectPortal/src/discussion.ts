/** Project discussion and Ask SME (F-51) on the platform's chats. A program's discussion is the
 * platform chat attached to its workspace folder: the same chat the platform's Chats panel shows on
 * the folder, read by whoever may view the folder (the portfolio's groups, the program's roles).
 * Viewers post too: the server lets anyone who reads a public chat write to it.
 *
 * The js-api has no chats client, so this calls the server's `/chats` routes. Threads, anchors and
 * the asked users are written into each message (./model/discussion). Asked users are told through
 * the platform's notifications, which the server allows only to callers who may edit users; others
 * find the question in the program's discussion and activity. Nothing leaves the platform from here,
 * but note that the server itself emails a chat's participants and watchers (platform findings). */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {portalDb, PortfolioRow, ProgramRow} from './generated/db';
import {notifyUser} from './notify';
import {programContext} from './portal';
import {programAddress} from './model/pinning';
import {savedRunsOf} from './model-runs';
import {documentsOf} from './documents';
import {textColumn} from './columns';
import {approverUsers, usersOfGroups} from './approvals';
import {
  Anchor, AnchorKind, AskSmeConfig, ChatMessage, composePost, Post, parsePost, plainMentions, smeContextFrom, smeNotice, smeQuestionBody, SmeContext, smeTargetOf,
  SmeTarget, Thread, threadsOf, threadTitle,
} from './model/discussion';

/** Calls one of the server's chat routes as the current user (or the one a test impersonates). */
async function chatApi(method: string, path: string, body?: any): Promise<any> {
  const r = await fetch(`${grok.dapi.root}/chats${path}`, {method, credentials: 'omit',
    headers: {'Content-Type': 'application/json', 'Authorization': grok.dapi.impersonationToken ?? grok.dapi.token},
    body: body === undefined ? undefined : JSON.stringify(body)});
  const text = await r.text();
  if (!r.ok)
    throw new Error(`chats ${method} ${path.split('?')[0]}: ${r.status} ${text.slice(0, 160)}`);
  try {
    return text ? JSON.parse(text) : null;
  } catch (_) {
    return text;
  }
}

/** The chats attached to [entityId] that the user may read. */
async function chatsOn(entityId: string): Promise<{id: string, name: string}[]> {
  const list = await chatApi('GET', `?entityId=${encodeURIComponent(entityId)}`);
  return (Array.isArray(list) ? list : []).map((c: any) => ({id: c.id, name: c.friendlyName ?? c.name ?? ''}));
}

/** Every message of one chat, oldest first. */
async function messagesOfChat(chatId: string): Promise<ChatMessage[]> {
  const list = await chatApi('GET', `/comments/${encodeURIComponent(chatId)}`);
  return (Array.isArray(list) ? list : []).map((c: any) => ({id: c.id, chatId, by: c.user?.id ?? null,
    byName: c.user?.friendlyName ?? c.user?.login ?? null, on: String(c.posted ?? ''), text: String(c.text ?? '')}));
}

export interface ProgramDiscussion {
  /** The program's workspace folder, which the chat hangs on; null when it has none. */
  folderId: string | null;
  /** The chat the portal posts to; null until the first post. */
  chatId: string | null;
  threads: Thread[];
  /** Why there is no discussion here, when there is none. */
  reason: string | null;
}

/** The program's discussion: the folder's chats (normally one), their messages as threads. */
export async function discussionOf(program: Pick<ProgramRow, 'space_id' | 'code'>): Promise<ProgramDiscussion> {
  const folderId = program.space_id ?? null;
  if (!folderId)
    return {folderId, chatId: null, threads: [], reason: `${program.code} has no workspace folder yet: the discussion lives on it`};
  const chats = await chatsOn(folderId);
  const messages = (await Promise.all(chats.map((c) => messagesOfChat(c.id)))).flat();
  return {folderId, chatId: chats[0]?.id ?? null, threads: threadsOf(messages), reason: null};
}

/** The program's posts, read back; empty where it has no folder or chat. */
export async function postsOf(program: Pick<ProgramRow, 'space_id' | 'code'>): Promise<Post[]> {
  return (await discussionOf(program)).threads.flatMap((t) => t.posts);
}

/** Posts [text] to the program's discussion as the current user, starting the chat on its folder
 * when it has none. Returns the post as read back. */
export async function postToProgram(program: ProgramRow, body: string, opts: {anchor?: Anchor | null, question?: boolean,
  askedOf?: string[], names?: string[]} = {}): Promise<{post: Post, chatId: string}> {
  if (!body.trim())
    throw new Error('Type something to post');
  const d = await discussionOf(program);
  if (!d.folderId)
    throw new Error(d.reason!);
  const text = composePost(body, opts);
  const chatId = d.chatId ?? crypto.randomUUID();
  const comment: any = {id: crypto.randomUUID(), text, chat: d.chatId ? {id: chatId} :
    {id: chatId, friendlyName: `Re: ${program.code}`, entities: [{id: d.folderId}]}};
  const saved = await chatApi('POST', `/comments/${encodeURIComponent(chatId)}`, comment);
  const me = DG.User.current();
  return {chatId, post: parsePost({id: saved?.id ?? comment.id, chatId, by: me.id, byName: me.friendlyName,
    on: String(saved?.posted ?? new Date().toISOString()), text})};
}

/** The user ids [target] names for [program]: approvers, a group's members, listed logins. */
export async function smeUsersOf(target: SmeTarget, portfolio: PortfolioRow, programId: string): Promise<string[]> {
  const ids: string[] = [];
  if (target.approvers)
    ids.push(...await approverUsers(portfolio, programId));
  for (const name of target.groups) {
    const g = await grok.dapi.groups.filter(`friendlyName = "${name}"`).first().catch(() => null) ??
      await grok.dapi.groups.filter(`name = "${name}"`).first().catch(() => null);
    if (!g)
      throw new Error(`askSme: no group "${name}"`);
    ids.push(...await usersOfGroups([g.id]));
  }
  for (const login of target.logins) {
    const u = await grok.dapi.users.filter(`login = "${login}"`).first().catch(() => null);
    if (!u)
      throw new Error(`askSme: no user "${login}"`);
    ids.push(u.id);
  }
  return [...new Set(ids)];
}

/** The portfolio's Ask SME config; the default when it declares none. */
export const askSmeOf = (portfolio: PortfolioRow): AskSmeConfig | null => (portfolio.ask_sme as AskSmeConfig | undefined) ?? null;

/** What a question about [programId] carries, read from where the program stands now. */
export async function smeContextOf(programId: string, anchor: Anchor | null = null): Promise<SmeContext> {
  const ctx = await programContext(programId);
  const origin = typeof window !== 'undefined' ? window.location.origin : '';
  return smeContextFrom({...ctx, anchor, appUrl: `${origin}${programAddress(ctx.portfolio.code, ctx.program.code)}`});
}

export interface AskSmeResult {
  post: Post;
  chatId: string;
  /** Who the question was put to, by id. */
  asked: string[];
  /** How the telling went, in one line. */
  line: string;
}

/** Asks an expert: the question with its context goes into the program's discussion, and the
 * people asked are told in the platform's notifications where the server lets the user notify
 * others. [to] overrides the portfolio's choice: logins, or user ids (the users mentioned in the question). */
export async function askSme(programId: string, question: string, context: SmeContext,
  to?: {logins?: string[], userIds?: string[]}): Promise<AskSmeResult> {
  const program = await portalDb.programs.get(programId);
  const portfolio = await portalDb.portfolios.get(program.portfolio_id);
  const named = !!to?.logins?.length || !!to?.userIds?.length;
  const target = named ? {channel: 'user' as const, groups: [], logins: to?.logins ?? [], approvers: false, notice: null} :
    smeTargetOf(askSmeOf(portfolio));
  const me = DG.User.current();
  let unresolved = '';
  const resolved = await smeUsersOf(target, portfolio, programId).catch((e: any) => {
    unresolved = `; who to ask could not be read (${e?.message ?? e})`;
    return [] as string[];
  });
  const asked = [...new Set([...(to?.userIds ?? []), ...resolved])].filter((id) => id !== me.id);
  const names = await Promise.all(asked.map(async (id) => (await grok.dapi.users.find(id).catch(() => null))?.friendlyName ?? id));
  const body = smeQuestionBody(question, context);
  const {post, chatId} = await postToProgram(program, body, {anchor: context.anchor ?? null, question: true, askedOf: asked,
    names: names.length <= 3 ? names : [...names.slice(0, 3), `${names.length - 3} more`]});
  if (!asked.length)
    return {post, chatId, asked, line: `Posted to the discussion; nobody else to tell${unresolved}`};
  const notice = smeNotice(me.friendlyName, context, question);
  const results = await Promise.all(asked.map((id) => notifyUser(id, notice)));
  const told = results.filter((r) => r.ok).length;
  const line = told === asked.length ? `Posted to the discussion; ${told} notified` :
    `Posted to the discussion; ${told} of ${asked.length} notified in the platform` +
    (told ? '' : ' (the server lets only user managers notify others): they find it in the program\'s discussion');
  return {post, chatId, asked, line: target.notice ? `${line}. ${target.notice}` : line};
}

/** An anchor by what a script knows: a criterion's code, a run's run_id or title, a document's
 * title (the newest version), or any of their ids. */
export async function resolveAnchor(programId: string, kind: AnchorKind, ref: string): Promise<Anchor> {
  const r = ref.trim();
  if (kind === 'criterion') {
    const c = (await programContext(programId)).criteria.find((x) => x.id === r || x.code === r || x.name === r);
    if (!c)
      throw new Error(`no criterion "${r}" on this program`);
    return {kind, id: c.id, label: c.name};
  }
  if (kind === 'run') {
    const run = (await savedRunsOf(programId)).map((x) => x.run).find((x) => x.id === r || x.run_id === r || x.title === r);
    if (!run)
      throw new Error(`no saved run "${r}" on this program`);
    return {kind, id: run.id, label: run.title ?? run.run_id ?? run.model_name};
  }
  const docs = (await documentsOf(programId)).filter((d) => d.id === r || d.title === r)
    .sort((a, b) => (b.doc_version ?? '').localeCompare(a.doc_version ?? '', undefined, {numeric: true}));
  if (!docs.length)
    throw new Error(`no document "${r}" on this program`);
  return {kind, id: docs[0].id, label: `${docs[0].title}${docs[0].doc_version ? ` v${docs[0].doc_version}` : ''}`};
}

/** The discussion as a table: thread, when, who, question, to, text, id. */
export function discussionFrame(threads: Thread[], name: string): DG.DataFrame {
  const rows = threads.flatMap((t) => [...t.posts].reverse().map((p) => ({t, p})));
  const df = DG.DataFrame.fromColumns([
    textColumn('thread', rows.map((x) => threadTitle(x.t))),
    DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, 'on', rows.map((x) => new Date(x.p.on))),
    textColumn('by', rows.map((x) => x.p.byName ?? '')),
    DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'question', rows.map((x) => x.p.question)),
    textColumn('text', rows.map((x) => plainMentions(x.p.body))),
    textColumn('id', rows.map((x) => x.p.id)),
  ]);
  df.name = name;
  return df;
}
