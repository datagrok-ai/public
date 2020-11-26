// Adds a view with a column combo box with int columns

let ccb = DG.ColumnComboBox.create(grok.data.demo.demog(), (c) => c.type === 'int');
let v = grok.shell.newView('demo: column combo box', [ccb]);