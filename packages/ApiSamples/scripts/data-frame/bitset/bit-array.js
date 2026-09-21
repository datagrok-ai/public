// Pure-JS bit masks (DG.BitArray) and applying one to the filter in a single copy

let df = grok.data.demo.demog();
grok.shell.addTableView(df);

let age = df.col('age');
let sex = df.col('sex');
let over40 = DG.BitArray.create(df.rowCount, (i) => age.get(i) > 40);
let female = DG.BitArray.create(df.rowCount, (i) => sex.get(i) === 'F');

// bitwise operations stay in JS: no Dart round trip
let mask = over40.and(female);

// one copy into the Dart-side filter; DG.BitSet.fromBitArray(mask) builds a standalone BitSet the same way
df.filter.copyFrom(mask);
grok.shell.info(`${df.filter.trueCount} of ${df.rowCount} rows match: ${mask.toString().substring(0, 16)}...`);
