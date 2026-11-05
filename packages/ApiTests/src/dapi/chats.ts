import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {after, category, expect, expectExceptionAsync, test} from '@datagrok-libraries/test/src/test';

// grok.dapi.chats and DG.ChatWidget over a chat on a throwaway project.
category('Dapi: chats', () => {
  let project: _DG.Project | undefined;
  const chatIds: string[] = [];

  async function newProject(): Promise<_DG.Project> {
    const p = DG.Project.create();
    p.name = `apitests-chats-${DG.Utils.randomString(6)}`;
    project = await grok.dapi.projects.save(p);
    return project;
  }

  function newChatId(): string {
    const id = crypto.randomUUID();
    chatIds.push(id);
    return id;
  }

  // A chat on no entity is a forum topic, whose first post emails the administrators.
  async function onProject(): Promise<{chat: {entityIds: string[]}}> {
    return {chat: {entityIds: [(project ?? await newProject()).id]}};
  }

  test('post, reply, edit, find by entity and by reference', async () => {
    const p = await newProject();
    const chatId = newChatId();
    const row = DG.ChatsDataSource.domainRowRef('assays', 'results', crypto.randomUUID());

    const question = await grok.dapi.chats.post(chatId, 'Which batch?', {chat: {friendlyName: 'ApiTests chat', entityIds: [p.id]}});
    expect(question.id != null, true, 'a posted comment has an id');
    const answer = await grok.dapi.chats.post(chatId, 'Batch 12', {parentId: question.id, refs: [{ref: row}]});
    expect(answer.parentId, question.id);

    const edited = await grok.dapi.chats.edit(chatId, question.id, 'Which batch failed?');
    expect(edited.edited != null, true, 'an edit is marked');

    const comments = await grok.dapi.chats.comments(chatId);
    expect(comments.map((c) => c.text).join('|'), 'Which batch failed?|Batch 12');

    expect((await grok.dapi.chats.forEntity(p.id)).some((c) => c.id === chatId), true, 'found by entity');
    expect((await grok.dapi.chats.forRef(row)).map((c) => c.id).join(), chatId, 'found by reference');
    expect((await grok.dapi.chats.find(chatId))?.friendlyName, 'ApiTests chat');
  });

  test('watch, unwatch, delete a comment', async () => {
    const chatId = newChatId();
    const c = await grok.dapi.chats.post(chatId, 'hello', await onProject());
    expect(await grok.dapi.chats.isWatched(chatId), true, 'the poster watches the chat');
    await grok.dapi.chats.unwatch(chatId);
    expect(await grok.dapi.chats.isWatched(chatId), false);
    await grok.dapi.chats.watch(chatId);
    expect(await grok.dapi.chats.isWatched(chatId), true);
    await grok.dapi.chats.deleteComment(c.id);
    expect((await grok.dapi.chats.comments(chatId)).length, 0);
  });

  test('import keeps the original date', async () => {
    const chatId = newChatId();
    await grok.dapi.chats.post(chatId, 'opened', await onProject());
    const c = await grok.dapi.chats.import(chatId, 'from the old system', {posted: new Date(Date.UTC(2020, 0, 1))});
    expect(new Date(c.posted!).getUTCFullYear(), 2020);
  });

  test('a comment needs text', async () => {
    await expectExceptionAsync(async () => {
      await grok.dapi.chats.post(newChatId(), '', await onProject());
    });
  });

  test('notify yourself', async () => {
    const n = await grok.dapi.users.notifications.send(DG.User.current().id, 'ApiTests chats notification');
    expect(n.sender, DG.User.current().id);
  });

  test('ChatWidget', async () => {
    const chatId = newChatId();
    await grok.dapi.chats.post(chatId, 'shown in a widget', await onProject());
    const w = await DG.ChatWidget.create(chatId);
    expect(w.root instanceof HTMLElement, true);
    w.detach();
    const p = project ?? await newProject();
    const e = DG.ChatWidget.forEntity(p);
    expect(e.root instanceof HTMLElement, true);
    e.detach();
  });

  after(async () => {
    for (const id of chatIds) {
      for (const c of await grok.dapi.chats.comments(id))
        try { await grok.dapi.chats.deleteComment(c.id); } catch (_) {}
    }
    if (project)
      await grok.dapi.projects.delete(project);
  });
}, {owner: 'aparamonov@datagrok.ai'});
