// This is a way to control values that can be set by a user.
// When TAGS_TOOLTIP is set, a combo box will be shown as a cell editor.

let view = grok.shell.addTableView(grok.data.demo.demog());
view.table.getCol('race').meta.choices = ['Asian', 'Black', 'Other'];
view.table.getCol('race').meta.autoChoices = false;
