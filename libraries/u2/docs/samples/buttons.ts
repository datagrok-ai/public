import {Control, Form, MultiChoiceInput, RadioInput, button, component, computed, dropDownButton, iconButton, notify}
  from '@datagrok-libraries/u2';

export function compoundActions(ids: string[]): Control {
  return component(() => {
    const compounds = new MultiChoiceInput({label: 'Compounds', items: ids});
    const view = new RadioInput({label: 'Show as', items: ['Grid', 'Cards'], value: 'Grid', buttons: true});
    const none = computed(() => compounds.value.length === 0 ? 'Pick at least one compound' : null);
    const exportAs = (format: string) => notify.info(`Exported ${compounds.value.length} compounds as ${format}`);
    return new Form().addAll([compounds, view]).addButtons((row) => row.append(
      button('Submit', () => notify.info(`Submitted ${compounds.value.join(', ')}`), {primary: true, disabled: none}),
      dropDownButton('Export', (menu) => menu.item('As CSV', () => exportAs('CSV')).item('As SDF', () => exportAs('SDF')),
        {split: () => exportAs('CSV')}),
      iconButton('sync', () => notify.info('Refreshed'), {tooltip: 'Refresh the list'}),
    ));
  }, {name: 'compound-actions'});
}

// Live demo
export function render(host: HTMLElement): Control {
  const demo = compoundActions(['CHEMBL25', 'CHEMBL112', 'CHEMBL521']);
  host.append(demo.root);
  return demo;
}
