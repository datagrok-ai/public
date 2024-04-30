// Adds a view with a column combo box with int and float columns

const ccb = DG.ColumnComboBox.create(grok.data.demo.demog(), (c) => ['int', 'double'].includes(c.type));
grok.shell.newView('demo: column combo box', [ccb]);