// Discussion on an entity: post to a chat attached to a project, reply, mention a user,
// reference a domain row, then show all chats on the project in a widget.

(async () => {
  const project = DG.Project.create();
  project.name = 'Chats sample';
  await grok.dapi.projects.save(project);

  // The first post creates the chat; it is visible to whoever can view the project.
  const chatId = crypto.randomUUID();
  const question = await grok.dapi.chats.post(chatId, 'Which batch failed QC?',
    {chat: {friendlyName: 'QC review', entityIds: [project.id]}});

  // A mention (entity markup) notifies the user and makes them watch the chat.
  const me = DG.User.current();
  const row = DG.ChatsDataSource.domainRowRef('assays', 'results', crypto.randomUUID());
  await grok.dapi.chats.post(chatId, `${me.toMarkup()} batch 12, see the linked result`,
    {parentId: question.id, refs: [{ref: row}]});

  for (const c of await grok.dapi.chats.comments(chatId))
    console.log(c.user?.friendlyName, c.parentId ? '(reply)' : '', c.text);

  // Chats are found from the entity they are on, and from what their comments reference.
  console.log((await grok.dapi.chats.forEntity(project.id)).map((c) => c.friendlyName));
  console.log((await grok.dapi.chats.forRef(row)).length);

  grok.shell.newView('Chats', [DG.ChatWidget.forEntity(project).root]);
})();
