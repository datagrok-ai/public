const dropDown = ui.dropDown(ui.iconFA('plus'), () => ui.div('Element'));
const dropDown2 = ui.dropDown(ui.iconFA('minus'), () => ui.list(['Element 1', 'Element 2', 'Element 3']));

dropDown.onElementClick.subscribe(() => grok.shell.info('Element click'));
dropDown.onExpand.subscribe(() => grok.shell.info('expand event'));

grok.shell.newView('DropDown', [ui.ribbonPanel([dropDown.root, dropDown2.root])]);
