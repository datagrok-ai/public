// An example of using Dart's future as JavaScript's promises
grok.data.openTable("e1792340-bd30-11e8-ed33-916dfb17fcaa")
    .then(t => grok.shell.info(t.name));
