// Context menu lifecycle: a viewer's onContextMenu fires while the menu is being built (add items
// there), grok.events.onContextMenuShown once the popup is in the DOM, onContextMenuClosed when it
// goes away.

const tv = grok.shell.addTableView(grok.data.demo.demog());
const barChart = tv.barChart();

barChart.onContextMenu.subscribe((menu) => menu.item('Custom item', () => grok.shell.info('Clicked')));
grok.events.onContextMenuShown.subscribe(() => grok.shell.info(`Menu shown: ${document.querySelectorAll('.d4-menu-popup .d4-menu-item').length} items`));
grok.events.onContextMenuClosed.subscribe(() => grok.shell.info('Menu closed'));
