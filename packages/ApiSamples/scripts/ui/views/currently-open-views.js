//tags: View
// Table views
grok.shell.add(grok.data.demo.demog());
grok.shell.add(grok.data.demo.wells());

// Views
grok.shell.newView('Test', [ui.h1('Regular View')]);
grok.shell.addView(DG.View.createByType(DG.View.FUNCTIONS));

// Iterate over table views and all views
for (let tv of grok.shell.tableViews)
  grok.shell.info(`Table view: ${tv.name}`);

for (let v of grok.shell.views)
  grok.shell.info(`View: ${v.name}`);
