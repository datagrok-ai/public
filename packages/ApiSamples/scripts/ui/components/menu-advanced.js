// Menu advanced item customization. For top menu, see menu.js

const dataFrame = grok.data.testData('demog', 50);
const columns = dataFrame.columns.toList().map(c => c.name);
const numericalColumns = columns.filter(c => dataFrame.columns.byName(c).isNumerical);
const menu = DG.Menu.popup();
menu.closeOnClick = false;

menu
  .header('Default header')
  .colorPalette(DG.Color.continuousSchemes, {
		onSelect: color => grok.shell.info(`Selected palette: ${DG.Color.continuousSchemes.indexOf(color) + 1}`),
	})
  .colorPalette(DG.Color.categoricalSchemes, {
		categorical: true,
		asGroup: 'Categorical Colors',
		getInitial: () => DG.Color.categoricalSchemes[0],
		onSelect: color => grok.shell.info(`Selected palette: ${DG.Color.categoricalSchemes.indexOf(color) + 1}`),
		onPreview: color => grok.shell.info(`Preview palette: ${DG.Color.categoricalSchemes.indexOf(color) + 1}`),
	})
	.separator()
	.header('Clickable header', {
		hasHoverEffect: true,
		onClick: () => grok.shell.info('Clicked on header!'),
		getDescription: () => 'Click me!',
	})
	.multiColumnSelector(dataFrame, {
		asGroup: 'multiColumSelector',
		initialValue: columns.slice(0, 5),
		onChange: grid => grok.shell.info(`Columns selected: ${grid}`),
	})
	.singleColumnSelector(dataFrame, {
		initialValue: numericalColumns.at(-2),
		onChange (grid, column, currentRowChanged) {
			if (currentRowChanged)
				grok.shell.info(`Column selected"`);
		},
		columnFilter: (c) => c.isNumerical,
	})
  .show();
