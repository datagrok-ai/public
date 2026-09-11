const test = require('node:test');
const assert = require('node:assert/strict');
const DG = require('./dg.cjs');

test('Qnum keeps the value and encodes the qualifier in the low mantissa bits', () => {
  const exact = DG.Qnum.exact(1.5);
  const less = DG.Qnum.less(1.5);
  const greater = DG.Qnum.greater(1.5);
  assert.equal(DG.Qnum.getValue(exact), 1.5);
  assert.equal(DG.Qnum.getValue(less), 1.5);
  assert.equal(DG.Qnum.getValue(greater), 1.5);
  assert.notEqual(DG.Qnum.getQ(less), DG.Qnum.getQ(exact));
  assert.notEqual(DG.Qnum.getQ(greater), DG.Qnum.getQ(exact));
  assert.notEqual(DG.Qnum.getQ(less), DG.Qnum.getQ(greater));
  assert.equal(DG.Qnum.getQ(DG.Qnum.create(2, DG.Qnum.getQ(less))), DG.Qnum.getQ(less));
});
