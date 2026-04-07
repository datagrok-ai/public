// Demonstrates programmatic use of getWidgetStatus() and descriptor.events:
// - outline each widget part with a random color
// - draw semi-transparent rectangles over hit areas on every render
// - list available events in a side panel
// - subscribe to 'd4-scatterplot-point-click' and show a balloon

let t = grok.data.demo.demog();
let view = grok.shell.addTableView(t);
let sp = view.scatterPlot({x: 'height', y: 'weight', color: 'sex', size: 'age'});

setTimeout(() => {
  let s = sp.getWidgetStatus();

  // Outline each part with a distinct random color
  for (let el of Object.values(s.parts))
    el.style.outline = `2px solid hsl(${Math.random() * 360 | 0}, 80%, 55%)`;

// Draw semi-transparent rectangles over hit areas on every render
  let overlay = s.parts['overlay'];
  sp.onEvent('d4-after-draw-scene').subscribe(() => {
    let ctx = overlay.getContext('2d');
    for (let r of Object.values(s.hitAreas)) {
      ctx.fillStyle = `hsla(${Math.random() * 360 | 0}, 80%, 55%, 0.25)`;
      ctx.fillRect(r.x, r.y, r.width, r.height);
    }
  });

// List available events in a side panel
  let events = sp.descriptor.events;
  let eventsPanel = ui.divV([
    ui.h2('Events'),
    ui.tableFromMap(Object.fromEntries(events.map((e) => [e.eventName, e.name]))),
  ]);
  view.dockManager.dock(eventsPanel, 'right', null, 'Widget Events', 0.25);

// Attach to the point-click event and show a balloon
  sp.onEvent('d4-scatterplot-point-click').subscribe((e) =>
    grok.shell.info(`Clicked row ${e.args.rowId}`)
  );
}, 1000)