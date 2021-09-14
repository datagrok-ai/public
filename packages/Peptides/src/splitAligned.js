//name: Split Sequence
//input: column peptideColumn {semType: alignedSequence}
//tags: panel
//output: widget result
export function splitAlignedPeptides(peptideColumn) {
  let splitPeptidesArray = [];
  let flag = true;
  let splitted;

  for (let peptideStr of peptideColumn.toList()) {
    splitted = peptideStr.split("-").slice(1, -1);
    if (flag) {
      for (let i = 0; i < splitted.length; i++) {
        splitPeptidesArray.push([]);
      }
      flag = false;
    }
    splitted.forEach((value, index) => {
      splitPeptidesArray[index].push(value);
    })
  }

  let columnsArray = splitPeptidesArray.map((v, i) => { return DG.Column.fromList('string', `a${i+3}`, v) });

  //FIXME: Show it in a different way/save it somewhere? 
  grok.shell.addTableView(DG.DataFrame.fromColumns(columnsArray));
}