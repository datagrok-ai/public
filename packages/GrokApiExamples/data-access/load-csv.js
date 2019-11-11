// Loading a csv file from the specified URL

gr.loadDataFrame('/demo/demog.csv')
  .then(t => gr.addTableView(t));
