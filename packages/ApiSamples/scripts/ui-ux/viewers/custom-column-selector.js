// Filters column by type, if filter null, then all columns match

let view = grok.shell.addTableView(grok.data.demo.demog());
sp = view.barChart();
sp.onEvent('d4-column-combo-box-popup-show').subscribe((args) => {
  args.preventDefault();

  grok.shell.info('Column filter: ' + args.args.comboBox.property.columnFilter);

  let filtered = Array.from(view.dataFrame.columns).filter(c => c.matches(args.args.comboBox.property.columnFilter));
  sp.setOptions({[args.args.comboBox.property.name]: filtered[0].name});
});