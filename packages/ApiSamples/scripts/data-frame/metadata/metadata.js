// Different ways of associating metadata with dataframe and columns and manipulating it

demog = grok.data.demo.demog();

// Tags are objects with string keys and string values. They are serialized together with the dataframe.
// To see tags associated with the dataframe, click on it and open the "Details" tab.

// Let's add some tags to the dataframe

demog.tags.foo = 'bar';
demog.tags['quux'] = 'quuz';
demog.setTag('baz', 'qux');

// There are multiple ways to iterate through the tags, just like with the JS Map

for (const key of demog.tags.keys())
  console.log(key, ':', demog.tags[key]);
for (const value of demog.tags.values())
  console.log(value);
for (const [key, value] of demog.tags)
  console.log(key, ':', value);
for (const [key, value] of demog.tags.entries())
  console.log(key, ':', value);

// There is a `description` tag coming with the synthetic "demog" dataframe.
// We can check if a certain tag is present, and delete this tag as well.

console.log('Description is present:', 'description' in demog.tags);
delete demog.tags.description;
console.log('Description is present:', 'description' in demog.tags);

// A special `.temp` variable can hold values of any types.
// Use it as a temporary store of the auxiliary data; it won't get serialized.
// Operate with the `.temp` variable just as we did with `.tags`.

demog.temp.corge = {a: 1, b: 2};

// The same concepts apply to columns

let age = demog.columns.byName('age');
age.temp.quuz = {a: 'column xuzzy'};
age.tags.garply = 'waldo';
demog.columns.byName('weight').tags.bar = 'plugh';

// Let's make sure we can get everything back

grok.shell.info(
  `${demog.tags.foo} <br/>` +
  `${demog.tags.quux} <br/>` +
  `${demog.tags.baz} <br/>` +
  `${JSON.stringify(demog.temp.corge)} <br/>` +
  `${JSON.stringify(age.temp.quuz)} <br/>` +
  `${age.tags.garply} <br/>` +
  `${demog.columns.byName('weight').tags.bar}`
);

// Hint: this JS iteration syntax is also available, but less preferred as it materializes tags to arrays first

for (const key of Object.keys(demog.tags))
  console.log(key, ':', demog.tags[key]);
for (const value of Object.values(demog.tags))
  console.log(value);
for (const [key, value] of Object.entries(demog.tags))
  console.log(key, ':', value);
for (const key in demog.tags)
  console.log(key, ':', demog.tags[key]);