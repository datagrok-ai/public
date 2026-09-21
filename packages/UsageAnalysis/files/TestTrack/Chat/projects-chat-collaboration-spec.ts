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

test('Projects / Chat collaboration: post, share, second-user reply, persistence', async ({page}) => {
  test.setTimeout(300_000);

  const COMMENT_A = 'First comment from A';
  const REPLY_B = 'Reply from B';
  const projectName = `ChatCollab${Date.now()}`;

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

    await test.step('A posts the first comment on the project chat', async () => {
      const proj = await createProject(page, projectName);
      projectId = proj.id;
      expect(await selectProjectAsCurrentObject(page, projectId)).toBe(true);
      await openChatsPane(page);
      await postComment(page, COMMENT_A);
      const chats = await getEntityChats(page, projectId);
      expect(chats.length, 'chat is linked to the project entity').toBe(1);
      chatId = chats[0].id;

      expect(await isCommentActionsHidden(page, 0)).toBe(false);
    });

    await test.step('A shares the project with B', async () => {
      await shareProjectWithUser(page, projectId, secondLogin,  false);
    });

    await test.step('B sees A comment, replies; per-author controls', async () => {
      await loginAsSecondUser(page);
      await setupSession(page);
      expect(await readLogin(page)).toBe(secondLogin);

      expect(await selectProjectAsCurrentObject(page, projectId),
        'B should have access to the shared project').toBe(true);
      await openChatsPane(page);
      await expect(commentByText(page, COMMENT_A)).toHaveCount(1);

      expect(await isCommentActionsHidden(page, 0)).toBe(true);

      await postComment(page, REPLY_B);
      expect(await commentTexts(page)).toEqual([COMMENT_A, REPLY_B]);
      expect(await isCommentActionsHidden(page, 0)).toBe(true);  
      expect(await isCommentActionsHidden(page, 1)).toBe(false); 
    });

    await test.step('A sees B reply after refresh; history persisted', async () => {
      await loginToDatagrok(page); 
      await setupSession(page);
      expect(await readLogin(page)).toBe(ownerLogin);

      expect(await selectProjectAsCurrentObject(page, projectId)).toBe(true);
      await openChatsPane(page);
      await expect(commentByText(page, REPLY_B)).toHaveCount(1);
      expect(await commentTexts(page)).toEqual([COMMENT_A, REPLY_B]);
      expect(await isCommentActionsHidden(page, 0)).toBe(false); 
      expect(await isCommentActionsHidden(page, 1)).toBe(true);  

      const persisted = await getChatCommentTexts(page, chatId);
      expect(persisted.slice().sort()).toEqual([COMMENT_A, REPLY_B].slice().sort());
    });
  } finally {

    await loginToDatagrok(page).catch(() => {});
    if (projectId) await deleteProjectById(page, projectId);
  }
});
