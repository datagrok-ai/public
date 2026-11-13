import {AsyncValue, Control, StatCard, Toolbar, component, div} from '@datagrok-libraries/u2';
import './stat-cards.css';

type Count = (abort: AbortSignal) => Promise<number>;

export function inventoryCards(compounds: Count, assays: Count, plates: Count): Control {
  return component(() => {
    const values = [new AsyncValue(compounds), new AsyncValue(assays), new AsyncValue(plates)];
    const cards = div([
      new StatCard({label: 'Compounds', source: values[0], icon: 'flask', delta: 0.12}),
      new StatCard({label: 'Assays this month', source: values[1], icon: 'vial', delta: -0.04}),
      new StatCard({label: 'Plates in stock', source: values[2], icon: 'th'}),
    ], {gap: 'm', cls: 'kpis'});
    const refresh = async () => {
      await Promise.all(values.map((v) => v.refresh()));
    };
    return [new Toolbar({ariaLabel: 'Inventory'}).addButton('Refresh', refresh, {icon: 'sync'}), cards];
  }, {name: 'inventory-cards'});
}

// Live demo
const after = (value: number, ms: number): Count => (abort) => new Promise((resolve, reject) => {
  const timer = setTimeout(() => resolve(value), ms);
  abort.addEventListener('abort', () => {
    clearTimeout(timer);
    reject(abort.reason);
  });
});

export function render(host: HTMLElement): Control {
  const demo = inventoryCards(after(48210, 600), after(1375, 900), () => Promise.reject(new Error('Inventory is offline')));
  host.append(demo.root);
  return demo;
}
