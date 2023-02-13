// Different ways for creating and populating columns with qualified numbers.

// Method 1: populate each value manually
let manual = DG.Column.qnum('manual', 3);
manual.set(0, DG.Qnum.greater(5));
manual.set(1, DG.Qnum.exact(5));
manual.set(2, DG.Qnum.less(5));

// Method 2: initialize with exact values (last parameter set to true)
// When last parameter is set to true, the last two bits of each number get set
// to 2 (DG.Qnum.EQUALS)
let exact = DG.Column.qnum('exact', 3, [1.123, 2.2345, 3.4567], true);

// Method 3: likely incorrect initialization with qualified numbers values (last parameter set to true)
// This method expects the qualifiers already present with the data, so be careful in order
// to not pass them by accident like in the code below.
// In the example below, the values passed will be interpreted as [>1.12, >2.23, <3.46]
// with the qualifiers that depend on the last two bits of the binary representation of that number.
let qAccidental = DG.Column.qnum('q-accidental', 3, [1.123, 2.2345, 3.4567], false);

// Method 4: correct initialization with qualified numbers values (last parameter set to true)
let qIntentional = DG.Column.qnum('q-intentional', 3, [
  DG.Qnum.exact(1.123),
  DG.Qnum.greater(2.2345),
  DG.Qnum.less(3.4567)], false);

let t = DG.DataFrame.fromColumns([manual, exact, qAccidental, qIntentional]);
grok.shell.addTableView(t);