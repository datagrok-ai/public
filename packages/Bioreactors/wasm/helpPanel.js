export function showHelpPanel() {
  grok.shell.windows.help.visible = true;

  const info = `# Try
  Vary inputs and check results
  # No-code
  Construction of complex phenomena simulators is provided by Datagrok WebAutosolver tool.
  # Model
  Only declarative equations description is required.
  # Essence
  Simulation of controlled fab-arm exchange kinetic mechanism is performed here.
  # Performance
  1000 times faster than the previous version.
  # Complexity
  Each time you change inputs, a system of 13 non-linear ordinary differential equations is solved.`;

  grok.shell.windows.help.showHelp(ui.markdown(info));
  
  grok.shell.windows.context.visible = true;

  grok.shell.windows.showContextPanel = false;
  grok.shell.windows.showProperties = false;
}