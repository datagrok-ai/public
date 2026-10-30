const test = require('node:test');
const assert = require('node:assert/strict');
const DG = require('./dg.cjs');

/** A table as a row reads it: a column it lacks throws, as DataFrame.getCol does. */
function table(columns) {
  return {
    get(name, idx) {
      if (!(name in columns))
        throw new Error(`No such column: ${String(name)}`);
      return columns[name][idx];
    },
  };
}

test('a row reads its columns as properties', () => {
  const row = new DG.Row(table({age: [30, 40]}), 1);
  assert.equal(row.age, 40);
  assert.equal(row.idx, 1);
  assert.throws(() => row.weight, /No such column: weight/);
});

test('a row is no promise: then and symbol keys answer undefined, so it can be awaited', async () => {
  const row = new DG.Row(table({age: [30]}), 0);
  assert.equal(row.then, undefined);
  assert.equal(row[Symbol.toPrimitive], undefined);
  assert.equal(row[Symbol.iterator], undefined);
  assert.equal(await Promise.resolve(row), row);
  assert.equal(await row, row);
});
