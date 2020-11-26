// Prevents the column combo box popup event

let view = grok.shell.addTableView(grok.data.demo.demog());

bc = view.barChart();
bc.onEvent('d4-column-combo-box-popup-show').subscribe((args) => {

    args.preventDefault();

    console.log(args.args.viewer)
    grok.shell.info(args.args.selectorName)
    console.log(args.args.comboBox)
});