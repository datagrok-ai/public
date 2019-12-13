// Demonstrates handling of Grok-originated events

demog = gr.testData('demog', 5000);
view = gr.addTableView(demog);

function info(s) { gr.balloon.info(s); }

gr.onProjectClosed(p => info(`${p.name}: closed`));
gr.onProjectModified(p => info(`${p.name}: modified`));
gr.onProjectOpened(p => info(`${p.name}: opened`));
gr.onProjectSaved(p => info(`${p.name}: saved`));
gr.onProjectUploaded(p => info(`${p.name}: uploaded`));

gr.onCurrentProjectChanged(p => info(`Current project changed: ${gr.project.name}`));

gr.onEvent('d4-current-view-changed', (_) => info('view changed'));
