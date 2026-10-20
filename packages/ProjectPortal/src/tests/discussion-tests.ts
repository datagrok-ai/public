import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, expectExceptionAsync, test} from '@datagrok-libraries/test/src/test';
import {portalDb, PortfolioRow, ProgramRow} from '../generated/db';
import {loadAllPortfolios, portfolioByCode} from '../portal';
import {
  ChatMessage, composePost, decisionsCiting, mentionedUserIds, missingExamplePosts, parsePost, plainMentions, postFeedText,
  smeContextFrom, smeContextLines, smeNotice, smeQuestionBody, smeTargetOf, threadsOf, threadTitle, threadUrl,
  validateAskSme,
} from '../model/discussion';
import {FeedRecords, programFeedEntries} from '../model/feed';
import {PAGE_PANELS, panelLayout} from '../model/personas';
import {askSme, discussionOf, postToProgram, resolveAnchor, smeContextOf} from '../discussion';
import {programFeedOf} from '../feed';
import {openAskSmeDialog} from '../views/discussion-dialog';
import {filterValue, ROW_LIMIT} from '../rows';
import {withStandUser} from './stand-users';
import {_package} from '../package-test';

const ID = '6e5f14a5-9603-4eea-80c9-d265e2610ec3';
const throws = (f: () => any): boolean => {
  try {
    f();
    return false;
  } catch (_) {
    return true;
  }
};
const msg = (id: string, text: string, on: string, by = 'u1'): ChatMessage => ({id, chatId: 'c1', by, on, text});

category('Portal: discussion rules', () => {
  test('a post says what it is about and whom it asks in its first line, and reads back', async () => {
    const anchor = {kind: 'criterion' as const, id: ID, label: 'hERG "margin"\n10 uM'};
    const text = composePost('  Which assay?  ', {anchor});
    expect(text, `Re: criterion "hERG 'margin' 10 uM" · pp:criterion:${ID}\nWhich assay?`);
    const p = parsePost(msg('m1', text, '2026-10-01T10:00:00Z'));
    expect(`${p.anchor?.kind} ${p.anchor?.id === ID} ${p.anchor?.label} ${p.question} ${p.body}`,
      'criterion true hERG \'margin\' 10 uM false Which assay?');
    const q = composePost('Rat only?', {anchor: {kind: 'run', id: 'r-1', label: 'Run 7'}, question: true, askedOf: ['u2', 'u3'],
      names: ['Ann', 'Bo']});
    expect(q.split('\n')[0], 'Question for Ann, Bo about run "Run 7" · pp:sme:u2,u3 pp:run:r-1');
    const qp = parsePost(msg('m2', q, '2026-10-01T11:00:00Z'));
    expect(`${qp.question} ${qp.askedOf.join()} ${qp.anchor?.kind}:${qp.anchor?.id} ${qp.anchor?.label}`, 'true u2,u3 run:r-1 Run 7');
    expect(composePost('plain'), 'plain', 'nothing to say: no header');
    // a message written in the platform's own Chats panel is a plain post
    const own = parsePost(msg('m3', 'Re: a thought\nwith two lines', '2026-10-01T12:00:00Z'));
    expect(`${own.anchor} ${own.question} ${own.body}`, 'null false Re: a thought\nwith two lines');
    expect(throws(() => composePost('x', {anchor: {kind: 'document', id: 'not an id!', label: 'D'}})), true);
  });

  test('mentions: the platform\'s user markup gives the asked ids and reads as @Name', async () => {
    const a = '<span>#{x.11111111-aaaa-bbbb-cccc-000000000001."Ann Lee"}</span>';
    const b = '<span>#{x.11111111-aaaa-bbbb-cccc-000000000002."Bo Wu"}</span>';
    const markup = `Is the panel enough, ${a} and ${b}? ${a} ran it`;
    expect(mentionedUserIds(markup).map((id) => id.slice(-1)).join(','), '1,2');
    expect(plainMentions(markup), 'Is the panel enough, @Ann Lee and @Bo Wu? @Ann Lee ran it');
    expect(mentionedUserIds('no one').length, 0);
    const notice = smeNotice('Kim', {program: {code: 'P1'}, portfolio: {code: 'EX'}, anchor: null} as any, markup);
    expect(notice.includes('@Ann Lee') && !notice.includes('#{x.'), true, notice);
  });

  test('threads: one per anchor and one general, newest thread first, posts oldest first', async () => {
    const a = composePost('first', {anchor: {kind: 'criterion', id: 'c-1', label: 'Potency'}});
    const b = composePost('reply', {anchor: {kind: 'criterion', id: 'c-1', label: 'Potency'}});
    const t = threadsOf([msg('1', a, '2026-10-01T10:00:00Z'), msg('2', 'hello', '2026-10-02T10:00:00Z'),
      msg('3', b, '2026-10-03T10:00:00Z')]);
    expect(t.map((x) => `${threadTitle(x)}:${x.posts.map((p) => p.body).join('+')}`).join(' | '), 'Criterion: Potency:first+reply | General:hello');
    const url = threadUrl('http://host/', 'chat-1', t[0].key);
    expect(url, 'http://host/forums/chat-1#pp-criterion-c-1');
    const evidence = new Map([['d1', [{target: url}]], ['d2', [{target: 'https://example.org'}]]]);
    expect(decisionsCiting(url, evidence).join(), 'd1', 'a decision citing the thread shows on it');
  });

  test('the portfolio picks who is asked; external is a placeholder that sends nothing', async () => {
    expect(validateAskSme(null).length, 0);
    expect(validateAskSme({channel: 'group'}).join(), 'askSme: channel group needs a group');
    expect(validateAskSme({channel: 'user', users: [' ']}).join(), 'askSme: channel user needs users');
    expect(validateAskSme({channel: 'fax' as any}).join(), 'askSme: channel "fax" is not one of discussion, group, user, external');
    const d = smeTargetOf(undefined);
    expect(`${d.channel} ${d.approvers} ${d.notice}`, 'discussion true null');
    expect(smeTargetOf({channel: 'group', group: 'Chemists'}).groups.join(), 'Chemists');
    expect(smeTargetOf({channel: 'user', users: ['ann ', 'bo']}).logins.join(), 'ann,bo');
    const ext = smeTargetOf({channel: 'external', external: {kind: 'teams'}});
    expect(`${ext.approvers} ${ext.notice?.startsWith('The teams channel is not connected: nothing leaves the platform')}`, 'true true');
  });

  test('a question carries the portfolio, program, stage, next gate, readiness, the item and the address', async () => {
    const c = smeContextFrom({portfolio: {code: 'EX-SMD', name: 'Example SM'}, program: {code: 'P1', name: 'Kinase', stage_id: 's1'},
      stages: [{id: 's1', name: 'Hit-to-Lead'}, {id: 's2', name: 'Lead Optimization'}],
      readiness: {stageId: 's2', empty: false, requiredMet: 1, requiredTotal: 2, score: 0.5, gaps: ['hERG not met']},
      anchor: {kind: 'criterion', id: 'c-1', label: 'hERG margin'}, appUrl: 'http://host/apps/x/P1'});
    expect(smeContextLines(c).join('\n'), ['Portfolio: Example SM (EX-SMD)', 'Program: P1 — Kinase', 'Stage: Hit-to-Lead',
      'Next gate: Lead Optimization: 1 of 2 required met or waived (50%)', 'Open: hERG not met',
      'About: criterion "hERG margin"', 'Open in the portal: http://host/apps/x/P1'].join('\n'));
    expect(smeQuestionBody('Is rat enough?', c).startsWith('Is rat enough?\n\nContext:\nPortfolio: Example SM (EX-SMD)'), true);
    expect(throws(() => smeQuestionBody(' ', c)), true, 'an empty question');
  });

  test('in the feed as discussion and questions; examples posted once; a panel shown in Simple mode', async () => {
    const q = parsePost(msg('m', composePost('Rat only?\n\nContext:\n...', {question: true, askedOf: ['u2'], names: ['Ann']}),
      '2026-10-01T10:00:00Z'));
    expect(postFeedText(q), 'Question to an expert: Rat only?');
    const r: FeedRecords = {program: {id: 'p1', code: 'P1'}, stages: [], decisions: [], evidence: [], runs: [], selections: [],
      documents: [], snapshots: [], roleChanges: [], syncs: [], news: [], discussion: [
        {id: 'a', on: '2026-10-01T10:00:00Z', by: 'u1', text: 'Discussion: x', question: false},
        {id: 'b', on: '2026-10-01T11:00:00Z', by: 'u1', text: postFeedText(q), question: true}]};
    expect(programFeedEntries(r).map((e) => e.kind).sort().join(), 'discussion,question');
    expect(missingExamplePosts([{text: 'A'}, {text: 'B'}], ['B\n\nContext:\nPortfolio: x', 'C']).map((x) => x.text).join(), 'A');
    expect(PAGE_PANELS.includes('discussion'), true);
    expect(panelLayout(null, null, 'simple', false).shown.includes('discussion'), true);
  });
});

/** Removes a post as an admin: the platform has no js-api for it. */
async function removePost(id: string): Promise<void> {
  await fetch(`${grok.dapi.root}/chats/comments/${id}`, {method: 'DELETE', credentials: 'omit',
    headers: {Authorization: grok.dapi.token}}).catch(() => undefined);
}

category('Portal: discussion', () => {
  let smd: PortfolioRow;
  let kin1: ProgramRow;
  const made: string[] = [];

  before(async () => {
    await loadAllPortfolios(_package);
    smd = await portfolioByCode('EX-SMD');
    kin1 = (await portalDb.programs.query({filter: `portfolio_id = ${filterValue(smd.id)} and code = "EXS-KIN1"`, limit: ROW_LIMIT}))[0];
  });

  after(async () => {
    for (const id of made.splice(0))
      await removePost(id);
    document.querySelectorAll('[data-u2-name="ask-sme-dialog"]').forEach((d) => d.remove());
  });

  test('a post about a criterion: on the folder\'s chat, its thread, the function\'s table and the feed', async () => {
    expect(!!kin1?.space_id, true, 'EXS-KIN1 needs a workspace folder');
    const text = `Test post ${Date.now()}`;
    const report: string = await grok.functions.call('ProjectPortal:postToDiscussion', {portfolio: 'EX-SMD', program: 'EXS-KIN1',
      text, aboutKind: 'criterion', about: 'herg_margin'});
    const id = /\(([^)]+)\)$/.exec(report)![1];
    made.push(id);
    expect(report.startsWith('EXS-KIN1: posted to the discussion about criterion "hERG IC50'), true, report);
    const d = await discussionOf(kin1);
    const t = d.threads.find((x) => x.posts.some((p) => p.id === id))!;
    expect(`${t.anchor?.kind} ${t.anchor?.label.startsWith('hERG')} ${t.posts.find((p) => p.id === id)!.by === DG.User.current().id}`,
      'criterion true true');
    const df: DG.DataFrame = await grok.functions.call('ProjectPortal:programDiscussion', {portfolio: 'EX-SMD', program: 'EXS-KIN1'});
    expect(df.col('text')!.toList().includes(text), true);
    const feed = await programFeedOf(smd, kin1, {kinds: ['discussion']});
    expect(feed.some((e) => e.text.endsWith(text)), true, 'in the activity');
    await expectExceptionAsync(async () => {
      await resolveAnchor(kin1.id, 'criterion', 'no_such_criterion');
    });
  });

  test('a viewer posts and reads; the admin reads it; a user of another portfolio reads nothing', async () => {
    const text = `Viewer post ${Date.now()}`;
    const viewer = await withStandUser('viewer', (u) => u.asUser(async () => {
      const {post} = await postToProgram(kin1, text);
      made.push(post.id);
      return (await discussionOf(kin1)).threads.flatMap((t) => t.posts).some((p) => p.id === post.id);
    }));
    expect(viewer, true, 'the viewer reads their own post');
    const seen = (await discussionOf(kin1)).threads.flatMap((t) => t.posts).find((p) => p.body === text);
    expect(!!seen && seen.by !== DG.User.current().id, true, 'the admin reads it, by the viewer');
    const outsider = await withStandUser('outsider', (u) => u.asUser(async () => (await discussionOf(kin1)).threads.length));
    expect(outsider, 0, 'EX-BIOD\'s viewer does not read EX-SMD\'s discussion');
  });

  test('Ask SME: the question and its context posted, the people asked told where the server allows', async () => {
    const context = await smeContextOf(kin1.id, await resolveAnchor(kin1.id, 'criterion', 'adme_panel'));
    expect(context.nextGate?.startsWith('Lead Optimization: '), true, context.nextGate ?? '');
    const r = await askSme(kin1.id, `Rat only? ${Date.now()}`, context, {logins: ['ppviewer']});
    made.push(r.post.id);
    expect(`${r.post.question} ${r.post.anchor?.kind} ${r.post.body.includes('Context:\nPortfolio: ')}`, 'true criterion true');
    expect(r.line, 'Posted to the discussion; 1 notified', 'an admin may notify others');
    // the users @mentioned in the question are the ones asked
    const viewer = await grok.dapi.users.filter('login = "ppviewer"').first();
    const markup = `Rat only, <span>#{x.${viewer.id}."${viewer.friendlyName}"}</span>? ${Date.now()}`;
    const m = await askSme(kin1.id, markup, context, {userIds: mentionedUserIds(markup)});
    made.push(m.post.id);
    expect(m.asked.join(','), viewer.id);
    expect(m.post.body.startsWith(markup), true, 'the mention is posted as the platform\'s markup');
    // a viewer asks the portfolio's approvers: posted; the server refuses to notify others for a non-manager
    const v = await withStandUser('viewer', (u) => u.asUser(() => askSme(kin1.id, `Viewer asks ${Date.now()}`, context)));
    if (v) {
      made.push(v.post.id);
      expect(v.line.startsWith('Posted to the discussion'), true, v.line);
    }
  });

  test('the Ask SME dialog shows the context and the portfolio\'s note', async () => {
    const context = await smeContextOf(kin1.id);
    const dlg = openAskSmeDialog({portfolio: smd, program: kin1, context});
    try {
      const shown = dlg.root.querySelector('[data-u2-name="sme-context"]')?.textContent ?? '';
      expect(shown.includes('Program: EXS-KIN1') && shown.includes('Next gate: ') && shown.includes('Open in the portal: '), true, shown);
      expect((dlg.root.querySelector('[data-u2-name="sme-note"]')?.textContent ?? '').startsWith('Example: an EX-SMD approver'), true);
      expect(!!dlg.root.querySelector('[data-u2-name="sme-ask-button"]'), true);
      // the question is u2's compose box with @-mentions; no logins field
      const box = dlg.root.querySelector('[data-u2-name="sme-question-input-composer"]');
      expect(box?.getAttribute('data-u2'), 'message-input');
      expect(dlg.root.querySelector('[data-u2-name="sme-question-input"]')?.getAttribute('contenteditable'), 'true');
      expect(!!dlg.root.querySelector('[data-u2-name="sme-login-input"]'), false);
    } finally {
      dlg.close();
    }
  });
});
