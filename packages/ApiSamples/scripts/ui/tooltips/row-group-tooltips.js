// Show a special tooltip for elements that represent multiple rows
// Note that other viewers also reflect the mouse-over state.

let t = grok.data.demo.demog();
let view = grok.shell.addTableView(t);
view.addViewer(DG.VIEWER.SCATTER_PLOT);
view.addViewer(DG.VIEWER.PIE_CHART);

let div = ui.divText('Hover to highlight odd rows');
$(div).on('mouseenter', () => ui.tooltip.showRowGroup(t, (i) => i % 2 === 0, 10, 10));
$(div).on('mouseleave', () => ui.tooltip.hide());
view.dockManager.dock(div);