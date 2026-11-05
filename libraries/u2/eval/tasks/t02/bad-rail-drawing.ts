import * as grok from 'datagrok-api/grok';
import {Form, TextInput} from '@datagrok-libraries/u2';

export function compoundSummary(compound: {id: string, name: string, smiles: string}): Form {
  const form = new Form();
  form.build(() => {
    const structure = new TextInput({label: 'Structure', value: compound.smiles, readOnly: true});
    structure.addOptions(grok.chem.drawMolecule(compound.smiles, 160, 100));
    form.addAll([
      new TextInput({label: 'ID', value: compound.id, readOnly: true}),
      new TextInput({label: 'Name', value: compound.name, readOnly: true}),
      structure,
    ]);
  });
  return form;
}
