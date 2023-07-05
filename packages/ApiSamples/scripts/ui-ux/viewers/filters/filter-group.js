let tv = grok.shell.addTableView(grok.data.demo.demog());
tv.filters({filters: [
    {type: DG.FILTER_TYPE.HISTOGRAM, column: 'height', min: 120, max: 150},
    {type: DG.FILTER_TYPE.FREE_TEXT},
    {type: DG.FILTER_TYPE.MULTI_VALUE, column: 'sex', selected: 'F'},
    {type: DG.FILTER_TYPE.CATEGORICAL, column: 'disease'},
  ]});