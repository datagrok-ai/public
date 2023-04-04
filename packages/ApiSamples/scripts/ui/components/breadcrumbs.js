const breadcrumbs = ui.breadcrumbs(['dev', 'datagrok', 'ai', 'apps', 'Curves']);
breadcrumbs.onPathClick.subscribe((value) => console.log(value));

grok.shell.newView('Breadcrumbs', [breadcrumbs.root]);
