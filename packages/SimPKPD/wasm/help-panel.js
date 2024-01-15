// Help panels

export function showHelpPanel() {

  grok.shell.windows.help.visible = true;

  const info = `# Try
  Interactive results when changing inputs
  # No-code
  Complex phenomena simulators are provided by Datagrok WebAutosolver tool.
  # Model
  Only declarative equations description is required.
  # Essence
  Two-compartment pharmacokinetic-pharmacodynamic (PK-PD) modeling is performed.
  # Performance
  Nonlinear system of differential equations within a few milliseconds
  # Learn more
  * [Compute](https://datagrok.ai/help/compute/)
  * [PK-PD](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7348046/)`;

  grok.shell.windows.help.showHelp(ui.markdown(info));
  
  grok.shell.windows.context.visible = true;  
  grok.shell.windows.showContextPanel = false;
  grok.shell.windows.showProperties = false;
  grok.shell.windows.help.visible = true;
}