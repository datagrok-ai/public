//tags: ViewLayout
// open a table, then do View | Layout | Clone

function info(action, info) {
  grok.shell.info(`${action}: ${info.viewState}`);
}

grok.events.onViewLayoutGenerated.subscribe((info) => info('generated', info));

grok.events.onViewLayoutApplying.subscribe((info) => info('applying', info));

grok.events.onViewLayoutApplied.subscribe((info) => info('applied', info));
