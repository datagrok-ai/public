grok.events.onViewChanged.subscribe((v) => grok.shell.balloon.info(`${v.name}: changed`));
grok.events.onTableAdded.subscribe((args) => { grok.shell.balloon.info(args.args.dataFrame.name); });