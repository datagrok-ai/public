/** Demonstrates TreeViewNode functionality */

let tree = ui.tree();

let g1 = tree.group('group 1', 1);
g1.enableCheckBox(true);
g1.item('item 1.1');
g1.item('item 1.2');

let g11 = g1.group('group 1.1', 1.1, false);
g11.item('item 1.1.1');
g11.item('item 1.1.2');

let g2 = tree.group('group 2', 2);
g2.enableCheckBox();
g2.item('item 2.1');
g2.item('item 2.2').enableCheckBox();

let groups = [g1, g11, g2];
let items = g1.items.slice()
    .concat(g11.items)
    .concat(g2.items);

let button = ui.button('Push', () => {
    grok.shell.info('Selected groups:\n' + groups.filter(g => g.checked).map(g => `${g.value}`).join(', '));
    grok.shell.info('Selected items:\n' + items.filter(i => i.checked).map(i => i.text).join(', '));
});

grok.shell.newView('TreeView Demo', [tree.root, button]);
