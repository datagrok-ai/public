// controlling the visibility of the UI parts

let windows = grok.shell.windows;
ui.dialog('Windows')
  .add(ui.input.bool('help', {value: windows.showHelp, onValueChanged: (inp, x) => windows.showHelp = x}))
  .add(ui.input.bool('sidebar', {value: windows.showSidebar, onValueChanged: (inp, x) => windows.showSidebar = x}))
  .add(ui.input.bool('toolbox', {value: windows.showToolbox, onValueChanged: (inp, x) => windows.showToolbox = x}))
  .add(ui.input.bool('status bar', {value: windows.showStatusBar, onValueChanged: (inp, x) => windows.showStatusBar = x}))
  .add(ui.input.bool('console', {value: windows.showConsole, onValueChanged: (inp, x) => windows.showConsole = x}))
  .add(ui.input.bool('properties', {value: windows.showProperties, onValueChanged: (inp, x) => windows.showProperties = x}))
  .add(ui.input.bool('tables', {value: windows.showTables, onValueChanged: (inp, x) => windows.showTables = x}))
  .add(ui.input.bool('columns', {value: windows.showColumns, onValueChanged: (inp, x) => windows.showColumns = x}))
  .add(ui.input.bool('variables', {value: windows.showVariables, onValueChanged: (inp, x) => windows.showVariables = x}))
  .show();