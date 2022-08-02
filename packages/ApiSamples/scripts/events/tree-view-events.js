/* Demonstrates TreeView events */

const tree = ui.tree();

const g1 = tree.group('group 1', {val: 1}, false);
g1.enableCheckBox();

g1.onNodeAdded.subscribe(_ => grok.shell.info(`Node added g1`));
tree.onChildNodeExpanding.subscribe(_ => grok.shell.info(`Child node expanding tree`));
tree.onChildNodeExpandedChanged.subscribe(_ => grok.shell.info(`Child node expanded changed tree`));
g1.onNodeCheckBoxToggled.subscribe(_ => grok.shell.info(`Node checkbox toggled g1`));
g1.onNodeContextMenu.subscribe(_ => grok.shell.info(`Node context menu g1`));
g1.onNodeEnter.subscribe(_ => grok.shell.info(`Node entered g1`));
tree.onNodeMouseEnter.subscribe(_ => grok.shell.info(`Node mouse enter tree`));
tree.onNodeMouseLeave.subscribe(_ => grok.shell.info(`Node mouse leave tree`));
tree.onSelectedNodeChanged.subscribe(_ => grok.shell.info(`Selected node changed tree`));

const i101 = g1.item('item 1.0.1');
i101.onSelected.subscribe(_ => grok.shell.info(`Node selected i101`));

const i102 = g1.item('item 1.0.2');
i102.enableCheckBox();
const g11 = g1.group('group 1.1', {val: 1.1}, false);
g11.item('item 1.1.1');
g11.item('item 1.1.2');

g1.onNodeExpanding.subscribe(_ => grok.shell.info(`Node expanding g1`));

grok.shell.newView('TreeView Demo', [tree.root]);
