import {Breadcrumbs, Control, Toolbar, component, notify, signal} from '@datagrok-libraries/u2';

export function platesPanel(): Control {
  return component(() => {
    const crumbs = new Breadcrumbs({items: ['Home', 'Projects', 'Kinase panel', 'Plates'],
      onClick: (index) => crumbs.setItems(crumbs.items.value.slice(0, index + 1))});
    const compact = signal(false);
    const bar = new Toolbar({ariaLabel: 'Plate actions'})
      .addButton('Run', () => notify.info('Running the analysis'), {icon: 'play', tooltip: 'Run the analysis'})
      .addButton('Save', () => notify.info('Saved'), {icon: 'save'})
      .addSeparator()
      .addToggle('Compact', compact, {tooltip: 'Denser rows'})
      .addMenu('Export', (m) => m.item('As CSV', () => notify.info('CSV')).item('As PDF', () => notify.info('PDF')));
    return [crumbs, bar];
  }, {name: 'plates-panel'});
}

// Live demo
export function render(host: HTMLElement): Control {
  const demo = platesPanel();
  host.append(demo.root);
  return demo;
}
