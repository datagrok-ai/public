// Table views
grok.shell.addTableView(grok.data.demo.demog());
grok.shell.addTableView(grok.data.demo.wells());

// Views
grok.shell.newView('Test', [ui.h1('Regular View')]);
grok.shell.addView(DG.View.createByType(DG.View.FUNCTIONS));

// Iterate over table views and all views
for (const tv of grok.shell.tableViews)
  grok.shell.info(`Table view: ${tv.name}`);

for (const v of grok.shell.views)
  grok.shell.info(`View: ${v.name}`);