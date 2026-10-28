import {Control, Toolbar, component, notify} from '@datagrok-libraries/u2';

export function importActions(): Control {
  return component(() => new Toolbar({ariaLabel: 'Import'})
    .addButton('Info', () => notify.info('Saved 12 compounds to KHL-12'), {icon: 'info-circle'})
    .addButton('Warning', () => notify.warning('2 rows had no SMILES and were skipped', {singleKey: 'skipped'}),
      {icon: 'exclamation-triangle'})
    .addButton('Error', () => notify.error('Could not reach the registry', {copyText: 'GET /api/registry: 503'}),
      {icon: 'times-circle'}), {name: 'import-actions'});
}

// Live demo
export function render(host: HTMLElement): Control {
  const demo = importActions();
  host.append(demo.root);
  return demo;
}
