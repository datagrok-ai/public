// open a table, then do View | Layout | Clone

function info(action, layout) {
  grok.shell.info(`${action}: ${layout.viewState}`);
}

grok.events.onViewLayoutGenerated.subscribe((layout) => info('generated', layout));

grok.events.onViewLayoutApplying.subscribe((layout) => info('applying', layout));

grok.events.onViewLayoutApplied.subscribe((layout) => info('applied', layout));
