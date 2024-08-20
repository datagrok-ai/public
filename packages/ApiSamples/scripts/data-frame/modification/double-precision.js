//tags: DataFrame, Column, modification

// Smallest number > 1 representable as 64 bit floating point number precisely
// See https://en.wikipedia.org/wiki/Double-precision_floating-point_format for more details
let preciseNumber = 1.0000000000000002;

let floats32 = new Float32Array(1);
let floats64 = new Float64Array(1);
floats32[0] = preciseNumber;
floats64[0] = preciseNumber;

let table = DG.DataFrame.fromColumns([
  DG.Column.fromFloat32Array('compact_floats', floats32),
  DG.Column.fromFloat64Array('precise_floats', floats64)
]);

console.log(table.get('compact_floats', 0)); // 1
console.log(table.get('precise_floats', 0)); // 1.0000000000000002

// Bump float column to double precision 
table.columns.byName('compact_floats').doublePrecision = true;
table.set('compact_floats', 0, preciseNumber);

console.log(table.get('compact_floats', 0)); // 1.0000000000000002
console.log(table.get('precise_floats', 0)); // 1.0000000000000002

