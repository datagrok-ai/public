// A convenient way to get reference datasets.
// Note that this is synchronous code.

function add(t) { grok.shell.addTableView(t); }

add(grok.data.demo.biosensor());

// other out-of-the-box datasets:
// grok.data.demo.demog(rows)
// grok.data.demo.wells(rows)
// grok.data.demo.randomWalk(rows, cols)