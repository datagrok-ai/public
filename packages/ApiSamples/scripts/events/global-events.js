grok.events.onViewChanged.subscribe((v) => grok.shell.info(`${v.name}: changed`));
grok.events.onTableAdded.subscribe((args) => { grok.shell.info(args.args.dataFrame.name); });
