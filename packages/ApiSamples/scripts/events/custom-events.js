grok.events.onCustomEvent('test').subscribe((v) => grok.shell.info(v));
grok.events.fireCustomEvent('test', '42');