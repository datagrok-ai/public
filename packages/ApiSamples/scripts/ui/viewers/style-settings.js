// Style-related context commands - https://datagrok.ai/help/visualize/viewers/#common-actions

const df = grok.data.demo.demog();
const view = grok.shell.addTableView(df);

const sp = view.scatterPlot({
  colorColumnName: 'sex',
  backColor: DG.Color.lightBlue,
});

// Sets the current state of viewer properties as the default configuration
// used to create new viewer instances of this type (works both for native and custom viewers).
// Corresponds to the "Pick Up / Apply | Set as Default"
// Takes optional arguments to specify whether [data] and [style] settings should be copied.
// By default, copies only style settings: setDefault() is the same as setDefault(false, true)
sp.props.setDefault(true, true);

// Clears the previously remembered default settings.
// Note that the reset method can be invoked on any instance of this viewer type.
// Corresponds to "Pick Up / Apply | Reset Default"
const sp1 = view.scatterPlot();
sp1.props.resetDefault();
