import {Control, Menu, MenuBar, component, notify, panel, signal, span} from '@datagrok-libraries/u2';

export function campaignsPanel(): Control {
  return component(() => {
    const sort = signal('Name');
    const bar = new MenuBar()
      .menu('File', (m) => m
        .item('New campaign...', () => notify.info('New campaign'), {shortcut: 'Ctrl+N', icon: 'plus'})
        .item('Open...', () => notify.info('Open'), {shortcut: 'Ctrl+O'})
        .separator()
        .group('Export').item('As CSV', () => notify.info('CSV')).item('As SDF', () => notify.info('SDF')).endGroup())
      .menu('View', (m) => {
        m.group('Sort by');
        for (const by of ['Name', 'Created', 'Molecules'])
          m.item(by, () => sort.value = by, {check: sort.value === by});
        m.endGroup().item('Show archived', () => notify.info('Archived shown'), {detail: '3 campaigns'});
      });
    const area = panel([span('Right-click here for the context menu')]);
    new Menu()
      .item('Copy', () => notify.info('Copied'), {shortcut: 'Ctrl+C'})
      .item('Paste', () => notify.info('Pasted'), {shortcut: 'Ctrl+V', disabled: 'Copy something first'})
      .separator()
      .item('Select all', () => notify.info('Selected all'))
      .bindContext(area);
    return [bar, area];
  }, {name: 'campaigns-panel'});
}

// Live demo
export function render(host: HTMLElement): Control {
  const demo = campaignsPanel();
  host.append(demo.root);
  return demo;
}
