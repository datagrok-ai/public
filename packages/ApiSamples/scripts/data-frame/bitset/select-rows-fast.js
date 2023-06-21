// Selecting or filtering rows by predicates

let demog = grok.data.demo.demog();
let selection = demog.selection;
grok.shell.addTableView(demog);

// method 1: init(predicate)
selection.init((i) => i % 2 === 0);

// method 2: passing notify = false (last parameter).
// You will need to call fireChanged() in order for UI to redraw.
for (let i = 0; i < 10; i++)
  selection.set(i, true, false);
selection.fireChanged();