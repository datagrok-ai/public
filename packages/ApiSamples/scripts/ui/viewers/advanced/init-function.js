// Custom viewer initialization function.
// Нou can attach bevahior to functions that could be saved in the layout, etc.
// For JS-based initialization, see init-script.js


// This function is dynamically registered for the sake of simplicity.
// Normally, you would use a package function
grok.functions.register({
  signature: 'void initScatterSquare(viewer v)',
  run: (v) => {
    v.onAfterDrawScene.subscribe((_) => {
      v.canvas.getContext('2d')
        .setFillStyle('red')
        .fillRect(100, 100, 50, 50);
    });
  }
});

grok.shell
  .addTableView(grok.data.demo.demog())
  .scatterPlot({initializationFunction: 'initScatterSquare'});