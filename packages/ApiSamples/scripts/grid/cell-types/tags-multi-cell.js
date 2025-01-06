// Additional cell types supported in the PowerPack package

const t = DG.DataFrame.create(2);

// MultiChoice: each cell becomes a list of checkboxes
const dis = t.columns.addNewString('disease');
dis.set(0, 'Anxiety, Glaucoma');
dis.set(1, 'Hepatitis A, Glaucoma');
dis.meta.choices = ['Anxiety', 'Hepatitis A', 'Glaucoma'];
dis.setTag(DG.TAGS.CELL_RENDERER, 'MultiChoice');

// Tags: renders each value as a colored tag
const site = t.columns.addNewString('site');
site.set(0, 'Buffalo, Orlando');
site.set(1, 'Buffalo, Los Angeles');
site.setTag(DG.TAGS.CELL_RENDERER, 'Tags');

grok.shell.addTableView(t);