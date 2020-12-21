// Adding new menu items for all context menus for viewers
grok.events
  .onContextMenu
  .subscribe((args) => {
    if (args.args.context instanceof DG.Viewer) {
      grok.shell.info(args.args.context.table.name);
      args.args.menu.item('Yo!', () => args.args.context.setOptions({title: 'Handled!'}));
    }
  });

// Let's add a viewer to make it more apparent
grok.shell
  .addTableView(grok.data.demo.demog())
  .scatterPlot({title: 'Right-click this plot'});