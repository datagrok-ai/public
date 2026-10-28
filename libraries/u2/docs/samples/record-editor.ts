import {ChoiceInput, Control, Form, NumberInput, TextInput, button, component, computed, notify, signal, statusLine}
  from '@datagrok-libraries/u2';

interface Compound {
  name: string;
  smiles: string;
  mw: number | null;
  stage: string | null;
}

export function compoundEditor(compound: Compound, save: (c: Compound) => void): Control {
  return component(() => {
    const saved = signal(compound);
    const name = new TextInput({label: 'Name', required: true});
    const smiles = new TextInput({label: 'SMILES', required: true});
    const mw = new NumberInput({label: 'MW', min: 0, postfix: 'g/mol'});
    const stage = new ChoiceInput({label: 'Stage', items: ['Hit', 'Lead', 'Candidate']});
    const load = (c: Compound) => {
      name.reset(c.name);
      smiles.reset(c.smiles);
      mw.reset(c.mw);
      stage.reset(c.stage);
    };
    const edited = computed((): Compound => ({name: name.value, smiles: smiles.value, mw: mw.value, stage: stage.value}));
    const changed = computed(() => (Object.keys(compound) as (keyof Compound)[]).some((k) => edited.value[k] !== saved.value[k]));
    const form = new Form().addAll([name, smiles, mw, stage]);
    const blocked = computed(() => form.problem.value ?? (changed.value ? null : 'No changes to save'));
    load(compound);
    return form.addButtons((row) => row.append(
      button('Save', () => {
        save(edited.value);
        saved.value = edited.value;
      }, {primary: true, disabled: blocked}),
      button('Revert', () => load(saved.value), {disabled: computed(() => changed.value ? null : 'Nothing to revert')}),
      statusLine(form.shownProblem, computed(() => changed.value ? 'Unsaved changes' : 'Saved')),
    ));
  }, {name: 'compound-editor'});
}

// Live demo
export function render(host: HTMLElement): Control {
  const demo = compoundEditor({name: 'Aspirin', smiles: 'CC(=O)OC1=CC=CC=C1C(=O)O', mw: 180.16, stage: 'Lead'},
    (c) => notify.info(`Saved ${c.name}`));
  host.append(demo.root);
  return demo;
}
