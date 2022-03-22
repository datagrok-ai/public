let nums = ['>5', '6', '<7'];
let t = DG.DataFrame.fromColumns([
  DG.Column.qnum('qnum', 3).init((i) => DG.Qnum.parse(nums[i]))
]);

let csv = t.toCsv({qualifierAsColumn: true});
grok.shell.addTableView(DG.DataFrame.fromCsv(csv));