// Demonstrates different link click modes in grid

const df = DG.DataFrame.fromColumns([
	DG.Column.fromList(DG.TYPE.STRING, 'open in new tab', ['To learn more about Datagrok open https://datagrok.ai/help/home']),
	DG.Column.fromList(DG.TYPE.STRING, 'open in context panel', ['To learn more about Datagrok open https://datagrok.ai/help/home']),
	DG.Column.fromList(DG.TYPE.STRING, 'custom', ['To learn more about Datagrok open https://datagrok.ai/help/home'])
]);
grok.shell.addTableView(df);

// Setting the column tags

df.col('open in new tab').tags[DG.TAGS.LINK_CLICK_BEHAVIOR] = 'Open in new tab';
df.col('open in context panel').tags[DG.TAGS.LINK_CLICK_BEHAVIOR] = 'Open in context panel';
df.col('custom').tags[DG.TAGS.LINK_CLICK_BEHAVIOR] = 'Custom';

// Displaying the clicked link - global and local grid events

grok.events.onGridCellLinkClicked.subscribe(eventData => grok.shell.info(`Global event - ${eventData.args.link}`)); // global
grok.shell.tv.grid.onGridCellLinkClicked.subscribe(eventData => grok.shell.info(`Local grid event - ${eventData.args.link}`)); // local