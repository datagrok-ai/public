// Demonstrates programmatic use of getWidgetStatus():
// - set color and size columns
// - outline each widget part with a random color
// - draw semi-transparent rectangles over hit areas on every render

let t = grok.data.demo.demog();
let view = grok.shell.addTableView(t);
let sp = view.scatterPlot({x: 'height', y: 'weight', color: 'sex', size: 'age'});

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
