// controlling the visibility of the UI parts

let windows = grok.shell.windows;
ui.dialog('Windows')
    .add(ui.boolInput('help', windows.showHelp, (x) => windows.showSidebar = x))
    .add(ui.boolInput('sidebar', windows.showSidebar, (x) => windows.showSidebar = x))
    .add(ui.boolInput('toolbox', windows.showToolbox, (x) => windows.showToolbox = x))
    .add(ui.boolInput('console', windows.showConsole, (x) => windows.showConsole = x))
    .add(ui.boolInput('properties', windows.showProperties, (x) => windows.showProperties = x))
    .add(ui.boolInput('tables', windows.showTables, (x) => windows.showTables = x))
    .add(ui.boolInput('columns', windows.showColumns, (x) => windows.showColumns = x))
    .add(ui.boolInput('variables', windows.showVariables, (x) => windows.showVariables = x))
    .show();