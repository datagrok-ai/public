import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export function newCompoundDialog(register: (name: string, smiles: string, mw: number | null) => Promise<void>): DG.Dialog {
  const name = ui.input.string('Name', {nullable: false});
  const smiles = ui.input.string('SMILES', {nullable: false});
  const mw = ui.input.float('MW', {min: 0});
  mw.addPostfix('g/mol');
  const inputs: DG.InputBase[] = [name, smiles, mw];
  const dialog = ui.dialog('New compound').add(ui.form(inputs));
  let busy = false;
  dialog.addButton('Register', async () => {
    if (busy || inputs.some((i) => i.validity !== null))
      return;
    busy = true;
    ui.setUpdateIndicator(dialog.root, true, 'Registering...');
    try {
      await register(name.value, smiles.value, mw.value);
      dialog.close();
    }
    catch (e) {
      grok.shell.error(e instanceof Error ? e.message : String(e));
    }
    finally {
      busy = false;
      ui.setUpdateIndicator(dialog.root, false);
    }
  });
  const registerButton = dialog.getButton('Register');
  const sync = () => {
    const problem = inputs.find((i) => i.validity !== null);
    registerButton.classList.toggle('d4-disabled', problem !== undefined);
    ui.tooltip.bind(registerButton, problem === undefined ? 'Register the compound' : `${problem.caption}: ${problem.validity}`);
  };
  for (const input of inputs)
    input.onChanged.subscribe(sync);
  sync();
  return dialog.show();
}
