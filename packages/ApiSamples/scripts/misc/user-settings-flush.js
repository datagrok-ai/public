// User settings storage - sending pending changes to the server without waiting for the periodic sync

const STORAGE_NAME = 'user-settings-flush-demo';

grok.userSettings.add(STORAGE_NAME, 'lastVisit', new Date().toISOString());
await grok.userSettings.flush();

const saved = await grok.dapi.userDataStorage.getValue(STORAGE_NAME, 'lastVisit');
grok.shell.info(`Saved on the server: ${saved}`);
