const tagsInput = ui.tagsInput('Tags', ['Apps', 'Demo', 'Curves'], false);

tagsInput.addTag('New');
tagsInput.removeTag('Demo');

tagsInput.onTagAdded.subscribe((value) => grok.shell.info(`${value} tag was added.`));
tagsInput.onTagRemoved.subscribe((value) => grok.shell.info(`${value} tag was removed.`));

const tagsInput2 = ui.tagsInput('Tags', ['Apps2', 'Demo2', 'Curves2'], true);
tagsInput2.setTags(['Apps', 'Demo', 'Viewers', 'Tree']);

grok.shell.newView('TagsInput', [ui.form([tagsInput, tagsInput2])]);
