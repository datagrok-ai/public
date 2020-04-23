grok.events.onViewChanged.subscribe((v) => grok.balloon.info(`${v.name}: changed`));
grok.events.onTableAdded.subscribe((args) => { grok.balloon.info(args.args.dataFrame.name); });