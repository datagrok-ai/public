const breadcrumbs = ui.breadcrumbs(['dev', 'datagrok', 'ai', 'apps', 'Curves']);
breadcrumbs.onPathClick.subscribe((value) => grok.shell.info(value));

grok.shell.newView('Breadcrumbs', [breadcrumbs.root]);
