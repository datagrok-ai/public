// controlling the visibility of the UI parts

let windows = grok.shell.windows;
ui.dialog('Windows')
  .add(ui.input.bool('help', {value: windows.showHelp, onValueChanged: (x) => windows.showHelp = x}))
  .add(ui.input.bool('sidebar', {value: windows.showSidebar, onValueChanged: (x) => windows.showSidebar = x}))
  .add(ui.input.bool('toolbox', {value: windows.showToolbox, onValueChanged: (x) => windows.showToolbox = x}))
  .add(ui.input.bool('status bar', {value: windows.showStatusBar, onValueChanged: (x) => windows.showStatusBar = x}))
  .add(ui.input.bool('console', {value: windows.showConsole, onValueChanged: (x) => windows.showConsole = x}))
  .add(ui.input.bool('properties', {value: windows.showProperties, onValueChanged: (x) => windows.showProperties = x}))
  .add(ui.input.bool('tables', {value: windows.showTables, onValueChanged: (x) => windows.showTables = x}))
  .add(ui.input.bool('columns', {value: windows.showColumns, onValueChanged: (x) => windows.showColumns = x}))
  .add(ui.input.bool('variables', {value: windows.showVariables, onValueChanged: (x) => windows.showVariables = x}))
  .show();