// open a table, then do View | Layout | Clone

function showInfo(action, info) {
  grok.shell.info(`${action}: ${info.viewState}`);
}

grok.events.onViewLayoutGenerated.subscribe((info) => showInfo('generated', info));

grok.events.onViewLayoutApplying.subscribe((info) => showInfo('applying', info));

grok.events.onViewLayoutApplied.subscribe((info) => showInfo('applied', info));
