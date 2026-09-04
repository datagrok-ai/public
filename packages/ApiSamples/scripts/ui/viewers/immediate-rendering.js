// Immediate rendering: the viewer repaints on a zero-delay timer instead of waiting for the next
// animation frame or the render delay. Meant for automated tests — a property set is painted by the
// time the test yields one macrotask; sets within one task still coalesce into a single repaint.

const tv = grok.shell.addTableView(grok.data.demo.demog());
const plot = tv.scatterPlot();
plot.immediateRendering = true;

plot.onViewerRendered.subscribe(() => console.log('scatter plot painted'));
plot.props.markerDefaultSize = 10;
await new Promise((r) => setTimeout(r, 0));   // painted by now
