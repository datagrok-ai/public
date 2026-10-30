import * as grok from 'datagrok-api/grok';
import {DynamicInput, Form, TextInput} from '@datagrok-libraries/u2';

export function compoundSummary(compound: {id: string, name: string, smiles: string}): Form {
  const form = new Form();
  form.build(() => form.addAll([
    new TextInput({label: 'ID', value: compound.id, readOnly: true}),
    new TextInput({label: 'Name', value: compound.name, readOnly: true}),
    new DynamicInput<string>({label: 'Structure', value: compound.smiles,
      render: (smiles) => smiles ? grok.chem.drawMolecule(smiles, 160, 100) : null}),
  ]));
  return form;
}
