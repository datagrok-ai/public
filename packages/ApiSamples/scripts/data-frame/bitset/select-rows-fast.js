// Selecting or filtering rows by predicates

let demog = grok.data.demo.demog();
grok.shell.add(demog);

demog.filter.init((i) => i % 2 === 0);
