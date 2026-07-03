// Comments / Chat collaboration on Projects — two real users.
// Manual scenario: UsageAnalysis/files/TestTrack/Chat/projects-chat-collaboration.md
//
// Single-page, token-switch two-user pattern (same as
// Projects/complex-share-second-user-spec.ts): User A drives the project's
// Context Panel "Chats" pane, shares to User B via grok.dapi.permissions.grant,
// then the page re-auths as B (loginAsSecondUser) to verify cross-user
// visibility / reply / per-author permissions, and switches back to A.
//
// Second user is MANDATORY — resolved from DATAGROK_AUTH_TOKEN_2 or the matching
// server's `key2:` in ~/.grok/config.yaml (spec-login.ts). The spec is skipped
// (not failed) when no second user is configured.
import {test, expect, Page} from '@playwright/test';
import {
  loginToDatagrok, loginAsSecondUser, getSecondUserLogin, specTestOptions,
} from '../spec-login';
import {
  gotoApp, setupSession, createProject, shareProjectWithUser, deleteProjectById,
  selectProjectAsCurrentObject, openChatsPane, postComment, commentByText, commentTexts,
  isCommentActionsHidden, getEntityChats, getChatCommentTexts,
} from './_helpers';

test.use(specTestOptions);

const readLogin = (page: Page): Promise<string | null> =>
  page.evaluate(() => (window as any).grok?.shell?.user?.login ?? null);

// Author label / avatar / timestamp of the rendered comment at `index`
// (chat.dart: .grok-comments-message-user, .grok-comments-message-picture, the
// .fa-clock timestamp icon whose title holds the formatted post time).
const commentRowMeta = (page: Page, index: number): Promise<{author: string; hasAvatar: boolean; timestamp: string}> =>
  page.evaluate((i) => {
    const q = (sel: string) => Array.from(document.querySelectorAll('.grok-comments-message-panel ' + sel));
    const authors = q('.grok-comments-message-user');
    const avatars = q('.grok-comments-message-picture');
    const stamps = q('.fa-clock');
    return {
      author: (authors[i]?.textContent ?? '').trim(),
      hasAvatar: !!avatars[i],
      timestamp: ((stamps[i] as HTMLElement)?.title ?? '').trim(),
    };
  }, index);

test('Chat / Projects chat collaboration: post, share, second-user reply, persistence', async ({page}) => {
  test.setTimeout(180_000);

  const COMMENT_A = 'First comment from A';
  const REPLY_B = 'Reply from B';
  const projectName = `ChatCollab${Date.now()}`;
  let authorA = '';
  let authorB = '';

  let secondLogin: string;
  try {
    secondLogin = await getSecondUserLogin();
  } catch (e: any) {
    test.skip(true, `No second user configured (${e?.message ?? e}). Set DATAGROK_AUTH_TOKEN_2 / DATAGROK_DEV_KEY_2 or a server \`key2:\`.`);
    return;
  }

  await gotoApp(page);
  await setupSession(page);
  const ownerLogin = await readLogin(page);
  expect(secondLogin, 'second user must differ from owner').not.toBe(ownerLogin);

  let projectId = '';
  let chatId = '';
  try {
    // Steps 1-2: A opens the project's Chats pane and posts the first comment.
    await test.step('A posts the first comment on the project chat', async () => {
      const proj = await createProject(page, projectName);
      projectId = proj.id;
      expect(await selectProjectAsCurrentObject(page, projectId)).toBe(true);
      await openChatsPane(page);
      await postComment(page, COMMENT_A);
      const chats = await getEntityChats(page, projectId);
      expect(chats.length, 'chat is linked to the project entity').toBe(1);
      chatId = chats[0].id;
      // Author owns their comment → edit/delete controls available to A.
      expect(await isCommentActionsHidden(page, 0)).toBe(false);
      // Feed row shows author, avatar and timestamp (.md step 2).
      const metaA = await commentRowMeta(page, 0);
      expect(metaA.author, 'A comment shows an author label').not.toBe('');
      expect(metaA.hasAvatar, 'A comment shows an avatar').toBe(true);
      expect(metaA.timestamp, 'A comment shows a timestamp').not.toBe('');
      authorA = metaA.author;
    });

    // Step 3: A shares the project with B (View access).
    await test.step('A shares the project with B', async () => {
      await shareProjectWithUser(page, projectId, secondLogin, /* edit */ false);
    });

    // Steps 4-7 + 11: page re-auths as B; B sees A's comment, replies; controls per-author.
    await test.step('B sees A comment, replies; per-author controls', async () => {
      await loginAsSecondUser(page);
      await setupSession(page);
      expect(await readLogin(page)).toBe(secondLogin);

      expect(await selectProjectAsCurrentObject(page, projectId),
        'B should have access to the shared project').toBe(true);
      await openChatsPane(page);
      await expect(commentByText(page, COMMENT_A)).toHaveCount(1);
      // B is not the author of A's comment → its edit/delete controls are hidden.
      expect(await isCommentActionsHidden(page, 0)).toBe(true);
      // B sees A's comment with A's author label (.md step 4).
      expect((await commentRowMeta(page, 0)).author, 'B sees A comment authored by A').toBe(authorA);

      await postComment(page, REPLY_B);
      expect(await commentTexts(page)).toEqual([COMMENT_A, REPLY_B]);
      expect(await isCommentActionsHidden(page, 0)).toBe(true);  // A's comment
      expect(await isCommentActionsHidden(page, 1)).toBe(false); // B's own reply
      const metaB = await commentRowMeta(page, 1);
      expect(metaB.hasAvatar, 'B reply shows an avatar').toBe(true);
      expect(metaB.timestamp, 'B reply shows a timestamp').not.toBe('');
      expect(metaB.author, 'B reply authored by B, distinct from A').not.toBe(authorA);
      expect(metaB.author, 'B reply shows an author label').not.toBe('');
      authorB = metaB.author;
      // .md step 5 (notification/inbox for the shared project / new comment) is
      // NOT asserted: notification generation is async and not reliably queryable
      // on the dev server; deferred for human review rather than a flaky green.
    });

    // Steps 8 + 12: back as A — sees B's reply after refresh; history persisted.
    await test.step('A sees B reply after refresh; history persisted', async () => {
      await loginToDatagrok(page); // re-auth as A (reloads the page)
      await setupSession(page);
      expect(await readLogin(page)).toBe(ownerLogin);

      expect(await selectProjectAsCurrentObject(page, projectId)).toBe(true);
      await openChatsPane(page);
      await expect(commentByText(page, REPLY_B)).toHaveCount(1);
      expect(await commentTexts(page)).toEqual([COMMENT_A, REPLY_B]);
      expect(await isCommentActionsHidden(page, 0)).toBe(false); // A's own comment
      expect(await isCommentActionsHidden(page, 1)).toBe(true);  // B's comment
      // Authors stay correct after the round-trip (.md step 7).
      expect((await commentRowMeta(page, 0)).author, 'comment 0 authored by A').toBe(authorA);
      expect((await commentRowMeta(page, 1)).author, 'comment 1 authored by B').toBe(authorB);

      // Persistence (server truth) — both comments stored; raw endpoint is unsorted.
      const persisted = await getChatCommentTexts(page, chatId);
      expect(persisted.slice().sort()).toEqual([COMMENT_A, REPLY_B].slice().sort());
    });
  } finally {
    // Ensure cleanup runs as the owner (A), who can delete the project.
    await loginToDatagrok(page).catch(() => {});
    if (projectId) await deleteProjectById(page, projectId);
  }
});
