//api: DG.DataFrame.appendAsync, DG.DataFrame.appendMergeAsync, DG.DataFrame.recalculateFormulaColumns, grok.data.appendTables
// append() moves rows only: a formula column the other table lacks stays empty for the appended rows.
// The async counterparts also calculate it, since evaluating a formula means running a function.

let t1 = DG.DataFrame.fromCsv(`make,volume,price
Honda,1.4,15000
Tesla,1.6,120000`);
let t2 = DG.DataFrame.fromCsv(`make,volume,price
BMW,1.7,60000
BMW,1.5,35000`);
await t1.columns.addNewCalculated('pricePerLiter', '${price} / ${volume}');

grok.shell.addTableView(await t1.appendAsync(t2));                // a new table
grok.shell.addTableView(await grok.data.appendTables([t1, t2]));  // any number of tables, as Data | Append Tables... does

let t3 = t1.clone();
await t3.appendMergeAsync(t2);                                    // in place, also adding the columns t3 lacks

t1.append(t2, true);                                              // rows only...
await t1.recalculateFormulaColumns();                             // ...so recalculate explicitly
grok.shell.addTableView(t1);
