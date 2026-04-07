grok.events.onServerMessage.subscribe((m) => console.log(m));

grok.dapi.admin.pushMessage('grok-foo', 'here we are');