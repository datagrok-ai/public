// Custom initialization script
// Use it to provide custom behavior for viewers that would survive saving/loading in a layout or project
// "v" parameter refers to the viewer. The following function gets executed:
// (v) => { <script> }
// For function-based initialization, see init-function.js


let view = grok.shell.addTableView(grok.data.demo.demog());

let plot = view.scatterPlot({
  onInitializedScript:
    'v.onAfterDrawScene.subscribe((_) => { ' +
    '  v.canvas.getContext(\'2d\').setFillStyle(\'red\').fillRect(100, 100, 50, 50); ' +
    '});'
});