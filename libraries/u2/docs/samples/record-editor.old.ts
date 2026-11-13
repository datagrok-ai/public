import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

interface Compound {
  name: string;
  smiles: string;
  mw: number | null;
  stage: string | null;
}

export function compoundEditor(compound: Compound, save: (c: Compound) => void): HTMLElement {
  let saved = compound;
  const name = ui.input.string('Name', {nullable: false});
  const smiles = ui.input.string('SMILES', {nullable: false});
  const mw = ui.input.float('MW', {min: 0});
  mw.addPostfix('g/mol');
  const stage = ui.input.choice('Stage', {items: ['Hit', 'Lead', 'Candidate']});
  const inputs: DG.InputBase[] = [name, smiles, mw, stage];
  const load = (c: Compound) => {
    name.value = c.name;
    smiles.value = c.smiles;
    mw.value = c.mw;
    stage.value = c.stage;
  };
  const edited = (): Compound => ({name: name.value, smiles: smiles.value, mw: mw.value, stage: stage.value});
  const changed = () => (Object.keys(compound) as (keyof Compound)[]).some((k) => edited()[k] !== saved[k]);
  const problem = () => {
    for (const input of inputs) {
      if (input.validity !== null)
        return `${input.caption}: ${input.validity}`;
    }
    return null;
  };
  let reason: string | null = null;
  const status = ui.divText('');
  const saveButton = ui.bigButton('Save', () => {
    if (reason !== null)
      return;
    save(edited());
    saved = edited();
    refresh();
  });
  ui.tooltip.bind(saveButton, () => reason);
  const revertButton = ui.button('Revert', () => {
    load(saved);
    refresh();
  });
  const refresh = () => {
    reason = problem() ?? (changed() ? null : 'No changes to save');
    saveButton.classList.toggle('d4-disabled', reason !== null);
    revertButton.classList.toggle('d4-disabled', !changed());
    status.textContent = problem() ?? (changed() ? 'Unsaved changes' : 'Saved');
  };
  load(compound);
  for (const input of inputs)
    input.onChanged.subscribe(() => refresh());
  refresh();
  return ui.divV([ui.form(inputs), ui.buttonsInput([saveButton, revertButton]), status]);
}
