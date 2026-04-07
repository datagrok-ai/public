// collaborative filtering - we are filtering out everything but sex=M

demog = grok.data.demo.demog();
let view = grok.shell.addTableView(demog);
view.dataFrame.onRowsFiltering.subscribe((_) => {
  let bs = view.dataFrame.filter;
  let sex = view.dataFrame.col('sex');
  bs.init((i) => sex.get(i) === 'M');
});