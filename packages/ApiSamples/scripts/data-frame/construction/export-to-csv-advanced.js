const nums = ['>5.13257', '6.016544', '<7.02'];
const t = DG.DataFrame.fromColumns([
  DG.Column.qnum('qnum', 3).init((i) => DG.Qnum.parse(nums[i]))
]);

// NUMBER FORMATS (https://datagrok.ai/help/discover/tags#numbers)
// The usual formats, such as "0.000" or "#.0000", are supported.
const columnFormats = {'qnum': '#0.000'};
const csv = t.toCsv({filteredRowsOnly: true, qualifierAsColumn: true, columnFormats: columnFormats});
const df1 = DG.DataFrame.fromCsv(csv);
const tv = grok.shell.addTableView(df1);

const grid = tv.grid;
grid.sort(['qnum'], [false]);
grid.col('qual(qnum)').visible = false;
// if grid is specified, row sorting will be applied as well
console.log(df1.toCsv({visibleColumnsOnly: true}, grid));