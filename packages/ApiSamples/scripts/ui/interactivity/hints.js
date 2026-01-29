const df = grok.data.demo.demog();
const view = grok.shell.addTableView(df);
let scatterPlot = view.scatterPlot();
let histogram = view.histogram();

// Place a hint on a visual component in your application. 
// Remove it by clicking on the root element event or by time-out.
const icon1 = $('div.d4-ribbon-item').has('i.svg-remove-selected-rows')[0];

let indicator1 = ui.hints.addHintIndicator(icon1, false, 4000);
let indicator2 = ui.hints.addHintIndicator(scatterPlot.root, true);

// Place a hint popup with a custom HTMLElement.
let msg = ui.divV([
  ui.h1('Title'),
  ui.divV([
    ui.divText('Click on Scatter plot to close hint indicator'),
    ui.link('Remove hint-indicator and hint-popup', ()=>{
      indicator2.click();
      popup.remove();
    }, '', '')
  ])
]);

let popup = ui.hints.addHint(scatterPlot.root, msg, 'left');