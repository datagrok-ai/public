// This is a way to control values that can be set by a user.
// When TAGS_TOOLTIP is set, a combo box will be shown as a cell editor.

let view = grok.addTableView(grok.testData('demog', 5000));
view.table.getCol('race').setTag(TAGS_CHOICES, '["Asian", "Black"]');
view.table.getCol('race').setTag(TAGS_AUTO_CHOICES, '["Asian", "Black"]');