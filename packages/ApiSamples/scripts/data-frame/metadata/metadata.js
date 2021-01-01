// Different ways of associating metadata with data frame and columns

demog = grok.data.demo.demog();

// tags are objects with string keys and string values. They are serialized together with the data frame.
// To see tags associtated with the data frame, click on it and open the "Details" tab.
demog.tags.foo = 'bar';
demog.setTag('foo1', 'bar1');

// A special "temp" variable can hold values of any types.
// Use it as a temporary store of the auxiliary data; it won't get serialized.
demog.temp.zoo = {a: 1, b: 2};

// The same concepts applies to columns:
let age = demog.columns.byName('age');
age.temp.foo = {a: 'column foo'};
age.tags.bar = 'bar';
demog.columns.byName('weight').tags.bar = 'bar';

// Let's make sure we can get everything back
grok.shell.info(`${demog.tags.foo} \n ${demog.temp.zoo} \n ${age.tags.bar} \n ${age.temp.foo} `);

// We also provide a bunch of convenient methods for finding columns
// var barColumns = demog.columns.