// Different ways to edit table's content, filter, selection, and metadata.

let demog = grok.data.demo.demog();
demog.columns.remove('sex');
let foo = demog.columns.addNew('foo', 'int');
let bar = demog.columns.byName('weight');
// Column's version is used to track changes of the column
demog.rows.removeAt(1, 3);
grok.shell.info(`bar was changed ${bar.version} times`);
demog.rows.insertAt(2, 2);
grok.shell.info(`bar was changed ${bar.version} times`);
demog.rows.addNew([777, 'studyX', 'NYC', 32, 'Spider', 'Net', 180, 80, 666]);
demog.rows.addNew().site = 'NY';

// alternative ways of setting values
foo.set(1, 777);
grok.shell.info(`foo was changed ${foo.version} times`);
demog.set('age', 1, 44);
//demog.currentRow.age = 33;

// bit set (same applies to filter)
demog.selection.invert();
demog.selection.set(5, false);
demog.selection.findNext(0, false);

// tags
demog.columns.byName('height').meta.units = 'm';
demog.columns.byName('weight').meta.units = 'kg';
foo = demog.columns.byName('height').meta.units;

grok.shell.addTableView(demog);
