import {Control, Dialog, Form, NumberInput, TextInput, button, component, notify} from '@datagrok-libraries/u2';

export function newCompoundDialog(register: (name: string, smiles: string, mw: number | null) => Promise<void>): Dialog {
  const name = new TextInput({label: 'Name', required: true});
  const smiles = new TextInput({label: 'SMILES', required: true});
  const mw = new NumberInput({label: 'MW', min: 0, postfix: 'g/mol'});
  return new Dialog('New compound')
    .add(new Form().addAll([name, smiles, mw]))
    .onOK(() => register(name.value, smiles.value, mw.value), {text: 'Register', primary: true})
    .show();
}

// Live demo
async function registerCompound(name: string): Promise<void> {
  await new Promise((resolve) => setTimeout(resolve, 800));
  if (name.trim().toLowerCase() === 'aspirin')
    throw new Error(`${name} is already registered as CHEMBL25`);
  notify.info(`Registered ${name}`);
}

export function render(host: HTMLElement): Control {
  const demo = component(() => button('New compound...', () => newCompoundDialog(registerCompound)));
  host.append(demo.root);
  return demo;
}
