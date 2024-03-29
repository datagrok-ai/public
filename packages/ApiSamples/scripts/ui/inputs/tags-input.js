const tagsInput = ui.input.tags('Tags', {tags: ['Apps', 'Demo', 'Curves'], showBtn: false});

tagsInput.addTag('New');
tagsInput.removeTag('Demo');

tagsInput.onTagAdded.subscribe((value) => grok.shell.info(`${value} tag was added.`));
tagsInput.onTagRemoved.subscribe((value) => grok.shell.info(`${value} tag was removed.`));

const tagsInput2 = ui.input.tags('Tags', {tags: ['Apps2', 'Demo2', 'Curves2'], showBtn: true});
tagsInput2.setTags(['Apps', 'Demo', 'Viewers', 'Tree']);

grok.shell.newView('TagsInput', [ui.form([tagsInput, tagsInput2])]);
