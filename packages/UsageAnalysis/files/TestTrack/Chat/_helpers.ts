import {Page, expect} from '@playwright/test';
import {loginToDatagrok} from '../spec-login';

export async function evalJs<T = any>(page: Page, script: string): Promise<T> {
  return page.evaluate(script) as Promise<T>;
}

export async function gotoApp(page: Page) {
  await loginToDatagrok(page);
}

export async function setupSession(page: Page) {
  await evalJs(page, `(() => {
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = true;
  })()`);
}

export async function createProject(page: Page, name: string): Promise<{id: string; name: string}> {
  return evalJs(page, `(async () => {
    const p = DG.Project.create();
    p.name = '${name}';
    await grok.dapi.projects.save(p);
    return {id: p.id, name: p.name};
  })()`);
}

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

export async function selectProjectAsCurrentObject(page: Page, projectId: string): Promise<boolean> {
  return evalJs(page, `(async () => {
    const p = await grok.dapi.projects.find('${projectId}');
    if (!p) return false;
    grok.shell.windows.showContextPanel = true;
    grok.shell.o = p;
    return true;
  })()`);
}

export async function openChatsPane(page: Page): Promise<void> {
  await evalJs(page, 'grok.shell.windows.showContextPanel = true');
  const header = page.locator('.d4-accordion-pane-header').filter({hasText: /^Chats$/}).first();
  await header.waitFor({state: 'visible', timeout: 20_000});
  const postBox = page.locator('.grok-comments-post > .grok-comments-post-input').first();
  for (let attempt = 0; attempt < 4; attempt++) {
    if (await postBox.isVisible().catch(() => false))
      return;
    await header.click();
    await postBox.waitFor({state: 'visible', timeout: 6_000}).catch(() => {});
  }
  await postBox.waitFor({state: 'visible', timeout: 8_000});
}

export async function postComment(page: Page, text: string): Promise<void> {
  const before = await page.locator('.grok-comments-message-text').count();
  const postBox = page.locator('.grok-comments-post > .grok-comments-post-input').first();
  await postBox.click();
  await postBox.fill(text);
  await page.keyboard.press('Enter');
  await expect.poll(async () => page.locator('.grok-comments-message-text').count(),
    {timeout: 15_000}).toBeGreaterThan(before);
  await expect(commentByText(page, text)).toHaveCount(1, {timeout: 15_000});
}

export function commentByText(page: Page, text: string) {
  return page.locator('.grok-comments-message-text').filter({hasText: new RegExp(`^${escapeRe(text)}$`)});
}

export async function commentTexts(page: Page): Promise<string[]> {
  return page.locator('.grok-comments-message-text').allTextContents();
}

export async function isCommentActionsHidden(page: Page, index: number): Promise<boolean> {
  return evalJs(page, `(() => {
    const blocks = Array.from(document.querySelectorAll('.grok-comments-message-edit'));
    return blocks[${index}] ? blocks[${index}].hidden : true;
  })()`);
}

export async function getEntityChats(page: Page, projectId: string): Promise<any[]> {
  return evalJs(page, `(async () => {
    const tok = window.localStorage.getItem('auth');
    return (await fetch('/api/chats?entityId=${projectId}', {headers: {Authorization: tok}})).json();
  })()`);
}

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
