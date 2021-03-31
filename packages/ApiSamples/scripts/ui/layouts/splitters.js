// Creating custom dialogs

let t = grok.data.demo.demog();
let block =
  $(ui.splitV([
    ui.textInput('', 'this \n is a \n multi-line \n script').root,
    ui.splitH([
      DG.Viewer.scatterPlot(t),
      DG.Viewer.histogram(t)])
  ])).css('flex-grow', '1');

ui.dialog('Windows')
  .add(block)
  .showModal(true);