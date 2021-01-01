// Different ways to edit table's content, filter, selection, and metadata.

let demog = grok.data.demo.demog();
demog.columns.remove('sex');
foo = demog.columns.addNew('foo', 'int');
demog.rows.removeAt(1, 3);
demog.rows.insertAt(2, 2);
demog.rows.addNew([777, 'studyX', 'NYC', 32, 'Spider', 'Net', 180, 80, 666]);
demog.rows.addNew().site = 'NY';

// alternative ways of setting values
foo.set(1, 777);
demog.set('age', 1, 44);
//demog.currentRow.age = 33;

// bit set (same applies to filter)
demog.selection.invert();
demog.selection.set(5, false);
demog.selection.findNext(0, false);

// tags
demog.columns.byName('height').setTag('units', 'm');
demog.columns.byName('weight').setTag('units', 'kg');
foo = demog.columns.byName('height').getTag('units');

grok.shell.addTableView(demog);
