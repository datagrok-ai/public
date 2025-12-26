// Prevents the column combo box popup event and shows custom popup

let view = grok.shell.addTableView(grok.data.demo.demog());

sp = view.scatterPlot();
sp.onEvent('d4-column-combo-box-popup-show').subscribe((args) => {

  args.preventDefault();

  let customPopup = ui.divText('My popup');
  let host = ui.showPopup(customPopup, args.args.comboBox.root, args.args.comboBox.vertical);

  host.addEventListener('click', function(event) {
    host.remove();
    grok.shell.info('Hi! You clicked on your popup!');
  });
});