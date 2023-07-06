let table = grok.data.demo.demog();
let view = grok.shell.addTableView(table);

let div = ui.divText('');
view.dockManager.dock(div);
view.dockManager.dock(div, 'left'); //dock to another side at any time

let viewer = DG.Viewer.fromType('Scatter Plot', table);
view.addViewer(viewer);

view.dockManager.dock(viewer, 'right');

div.append(ui.bigButton('close', (_) => {
  view.dockManager.close(div);
}));