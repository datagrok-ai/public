grok.events.onViewChanged.subscribe((v) => grok.shell.info(v.message));
grok.events.onTableAdded.subscribe((args) => {
  grok.shell.info(args.args.dataFrame.name);
});
