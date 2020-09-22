// Scroll to pixels.

let grid = grok.shell.addTableView(grok.data.demo.randomWalk(100, 100)).grid;
rxjs.interval(1000).pipe(rxjs.operators.startWith(0)).subscribe((i) => {
    grid.scrollToPixels(i, i);
});