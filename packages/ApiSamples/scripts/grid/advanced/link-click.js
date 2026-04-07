// Demonstrates different link click modes in grid

const df = DG.DataFrame.fromColumns([
  DG.Column.fromList(DG.TYPE.STRING, 'open in new tab', ['To learn more about Datagrok open https://datagrok.ai/help/home']),
  DG.Column.fromList(DG.TYPE.STRING, 'open in context panel', ['To learn more about Datagrok open https://datagrok.ai/help/home']),
  DG.Column.fromList(DG.TYPE.STRING, 'custom', ['To learn more about Datagrok open https://datagrok.ai/help/home'])
]);
grok.shell.addTableView(df);

// Setting the column tags

df.col('open in new tab').meta.linkClickBehavior = DG.LINK_CLICK_BEHAVIOR.OPEN_IN_NEW_TAB;
df.col('open in context panel').meta.linkClickBehavior = DG.LINK_CLICK_BEHAVIOR.OPEN_IN_CONTEXT_PANEL;
df.col('custom').meta.linkClickBehavior = DG.LINK_CLICK_BEHAVIOR.CUSTOM;

// Displaying the clicked link - global and local grid events

grok.events.onGridCellLinkClicked.subscribe(eventData => grok.shell.info(`Global event - ${eventData.args.link}`)); // global
grok.shell.tv.grid.onGridCellLinkClicked.subscribe(eventData => grok.shell.info(`Local grid event - ${eventData.args.link}`)); // local