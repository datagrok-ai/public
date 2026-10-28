import {Control, Splitter, VirtualList, component, computed, divV, emptyState, show, span, tableFromMap}
  from '@datagrok-libraries/u2';
import './master-detail.css';

interface Compound {
  id: string;
  name: string;
  mw: number;
  stage: string;
}

export function compoundBrowser(compounds: Compound[]): Control {
  return component(() => {
    const list = new VirtualList<Compound>({itemHeight: 28, keyOf: (c) => c.id, render: (c) => span(`${c.id} · ${c.name}`)});
    list.setItems(compounds);
    const current = computed(() => compounds[list.selectedIndex.value] ?? null);
    const field = (read: (c: Compound) => string) => span(computed(() => current.value === null ? '' : read(current.value)));
    const details = divV([], 'compound-details');
    show(details, computed(() => current.value !== null),
      () => tableFromMap({'Name': field((c) => c.name), 'MW': field((c) => `${c.mw.toFixed(2)} g/mol`), 'Stage': field((c) => c.stage)}),
      () => emptyState('Select a compound to see its details'));
    const split = new Splitter([list, details], {direction: 'horizontal', sizes: [45, 55]});
    split.root.classList.add('compound-browser');
    return split;
  }, {name: 'compound-browser'});
}

// Live demo
const STAGES = ['Hit', 'Lead', 'Candidate'];

export function render(host: HTMLElement): Control {
  const browser = compoundBrowser(Array.from({length: 2000}, (_, i) =>
    ({id: `CPD-${String(i + 1).padStart(4, '0')}`, name: `Compound ${i + 1}`, mw: 180 + (i * 37) % 420, stage: STAGES[i % 3]})));
  host.append(browser.root);
  return browser;
}
