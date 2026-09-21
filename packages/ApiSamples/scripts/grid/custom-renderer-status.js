//api: DG.Widget.addStatusProvider, DG.Widget.removeStatusProvider, DG.Viewer.getWidgetStatus
// Add live readings alongside the grid's native cells and headers.
const view = grok.shell.addTableView(grok.data.demo.demog(100));
view.grid.addStatusProvider('selection-summary', () => ({
  values: {'custom selected rows': view.dataFrame.selection.trueCount},
}));
grok.shell.info(view.grid.getWidgetStatus().values['custom selected rows']);
// Remove the contribution when the custom renderer is removed.
view.grid.removeStatusProvider('selection-summary');
