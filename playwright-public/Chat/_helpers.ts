// Helpers for the Chat collaboration spec.
// Auth + second-user resolution live in ../spec-login.ts (token injection,
// the path `grok test` supports); session bootstrap mirrors Projects/_helpers.
// Chat-widget selectors are transcribed from the client source
// core/client/xamgle/lib/src/features/chat/chat.dart (grok-comments-* widget).
import {Page, expect} from '@playwright/test';
import {loginToDatagrok} from '@datagrok-libraries/test/src/playwright/spec-login';

export async function evalJs<T = any>(page: Page, script: string): Promise<T> {
  return page.evaluate(script) as Promise<T>;
}

/** Token-injection login (DATAGROK_AUTH_TOKEN) — same path `grok test` uses. */
export async function gotoApp(page: Page) {
  await loginToDatagrok(page);
}

export async function setupSession(page: Page) {
  await evalJs(page, `(() => {
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = true;
  })()`);
}

/** Create an empty server-side project owned by the current user; returns id + name. */
export async function createProject(page: Page, name: string): Promise<{id: string; name: string}> {
  return evalJs(page, `(async () => {
    const p = DG.Project.create();
    p.name = '${name}';
    await grok.dapi.projects.save(p);
    return {id: p.id, name: p.name};
  })()`);
}

/** Grant a user's personal group VIEW (or Edit) on a project — JS-API share path. */
export async function shareProjectWithUser(page: Page, projectId: string, login: string, edit = false): Promise<void> {
  await evalJs(page, `(async () => {
    const target = await grok.dapi.users.filter('login = "${login}"').first();
    const p = await grok.dapi.projects.find('${projectId}');
    await grok.dapi.permissions.grant(p, target.group, ${edit});
  })()`);
}

export async function deleteProjectById(page: Page, projectId: string): Promise<void> {
  await evalJs(page, `(async () => {
    const tok = window.localStorage.getItem('auth');
    const chats = await (await fetch('/api/chats?entityId=${projectId}', {headers: {Authorization: tok}})).json();
    for (const c of chats) await fetch('/api/chats/' + c.id, {method: 'DELETE', headers: {Authorization: tok}});
    const p = await grok.dapi.projects.find('${projectId}');
    if (p) await grok.dapi.projects.delete(p);
  })()`).catch(() => {});
}

/** Make a project the current object so its Context Panel (incl. Chats) renders. Returns false if no access. */
export async function selectProjectAsCurrentObject(page: Page, projectId: string): Promise<boolean> {
  return evalJs(page, `(async () => {
    const p = await grok.dapi.projects.find('${projectId}');
    if (!p) return false;
    grok.shell.windows.showContextPanel = true;
    grok.shell.o = p;
    return true;
  })()`);
}

/**
 * The chat compose box: PowerPack's contenteditable `.u2-msg-editor` when its chat input is
 * installed, else `textarea.grok-comments-post-input`. Scoped to `.grok-comments-post`, because
 * each rendered comment carries its own hidden inline-edit box.
 */
export const POST_BOX = '.grok-comments-post .u2-msg-editor, .grok-comments-post textarea.grok-comments-post-input';

/**
 * Expand the Context Panel "Chats" pane and wait for its post box.
 * The pane content is a lazy `ui.wait` and may already be expanded, so poll:
 * each header click toggles, re-check until the post box is actually visible.
 */
export async function openChatsPane(page: Page): Promise<void> {
  await evalJs(page, 'grok.shell.windows.showContextPanel = true');
  const header = page.locator('.d4-accordion-pane-header').filter({hasText: /^Chats$/}).first();
  await header.waitFor({state: 'visible', timeout: 20_000});
  const postBox = page.locator(POST_BOX).first();
  for (let attempt = 0; attempt < 4; attempt++) {
    if (await postBox.isVisible().catch(() => false))
      return;
    await header.click();
    await postBox.waitFor({state: 'visible', timeout: 6_000}).catch(() => {});
  }
  await postBox.waitFor({state: 'visible', timeout: 8_000});
}

/** Type a comment into the chat post box and send it (the pane uses send-on-Enter). */
export async function postComment(page: Page, text: string): Promise<void> {
  const before = await page.locator('.grok-comments-message-text').count();
  const postBox = page.locator(POST_BOX).first();
  await postBox.click();
  await postBox.fill(text);
  await page.keyboard.press('Enter');
  await expect.poll(async () => page.locator('.grok-comments-message-text').count(),
    {timeout: 15_000}).toBeGreaterThan(before);
  await expect(commentByText(page, text)).toHaveCount(1, {timeout: 15_000});
}

/** The rendered comment body, without the "(edited)" marker that shares its container. */
const COMMENT_BODY = '.grok-comments-message-text > :not(.grok-comments-message-edited)';

/** Locator for a comment by its exact rendered text. */
export function commentByText(page: Page, text: string) {
  return page.locator(COMMENT_BODY).filter({hasText: new RegExp(`^${escapeRe(text)}$`)});
}

/** Texts of all comments currently rendered in the chat feed, in order. */
export async function commentTexts(page: Page): Promise<string[]> {
  return page.locator(COMMENT_BODY).allTextContents();
}

/**
 * Whether the edit/delete action block is hidden for the comment at `index`.
 * The widget hides it for everyone except the comment's author
 * (`..hidden = c.user.id != Auth.current.id`, chat.dart) — the per-author
 * permission signal asserted in the collaboration scenario. Each comment also has a
 * reply block before it that everyone sees, so take the last block of each comment.
 */
export async function isCommentActionsHidden(page: Page, index: number): Promise<boolean> {
  return evalJs(page, `(() => {
    const blocks = Array.from(document.querySelectorAll('.grok-comments-message-block'))
      .map((m) => Array.from(m.querySelectorAll(':scope > .grok-comments-message-edit')).pop());
    return blocks[${index}] ? blocks[${index}].hidden : true;
  })()`);
}

/** REST: chat(s) linked to a project entity (uses the page's own auth token). */
export async function getEntityChats(page: Page, projectId: string): Promise<any[]> {
  return evalJs(page, `(async () => {
    const tok = window.localStorage.getItem('auth');
    return (await fetch('/api/chats?entityId=${projectId}', {headers: {Authorization: tok}})).json();
  })()`);
}

/** REST: comment texts of a chat (server-truth persistence check). */
export async function getChatCommentTexts(page: Page, chatId: string): Promise<string[]> {
  return evalJs(page, `(async () => {
    const tok = window.localStorage.getItem('auth');
    const comments = await (await fetch('/api/chats/comments/${chatId}', {headers: {Authorization: tok}})).json();
    return comments.map(c => c.text);
  })()`);
}

function escapeRe(s: string): string {
  return s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
}
