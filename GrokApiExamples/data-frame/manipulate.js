// Different ways to edit table's content, filter, selection, and metadata.

var demog = gr.testData('demog', 5000);
demog.cols.remove('sex');
foo = demog.cols.addNew('foo', 'int');
demog.rows.removeAt(1, 3);
demog.rows.insertAt(2, 2);
demog.rows.addNew(['Spiderman', 'studyX', 'NYC', 32, 'Spider', 'Net', new Date(2020), 180, 80, 666]);
demog.rows.addNew().subj = 'Iron Man';

// alternative ways of setting values
foo.set(1, 777);
demog.set('age', 1, 44);
//demog.currentRow.age = 33;

// bit set (same applies to filter)
demog.selection.invert();
demog.selection.set(5, false);
demog.selection.findNext(0, false);

// tags
demog.cols.byName('height').setTag('units', 'm');
demog.cols.byName('weight').setTag('units', 'kg');
foo = demog.cols.byName('height').getTag('units');

gr.addTableView(demog);