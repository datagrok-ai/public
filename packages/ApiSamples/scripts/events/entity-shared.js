// React to an entity being shared — e.g. keep permissions of dependent entities in sync

let sub = grok.events.onEntityShared.subscribe((entity) =>
  grok.shell.info(`Shared: ${entity.name}`));
