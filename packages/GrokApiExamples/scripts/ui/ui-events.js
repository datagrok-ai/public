// Demonstrates handling of Grok-originated events

demog = grok.data.testData('demog', 5000);
view = grok.shell.addTableView(demog);

function info(s) { grok.shell.balloon.info(s); }

grok.onProjectClosed(p => info(`${p.name}: closed`));
grok.onProjectModified(p => info(`${p.name}: modified`));
grok.onProjectOpened(p => info(`${p.name}: opened`));
grok.onProjectSaved(p => info(`${p.name}: saved`));
grok.onProjectUploaded(p => info(`${p.name}: uploaded`));

grok.onCurrentProjectChanged(p => info(`Current project changed: ${grok.project.name}`));

grok.onEvent('d4-current-view-changed', (_) => info('view changed'));
