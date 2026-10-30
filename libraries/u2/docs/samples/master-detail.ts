import {Control, MasterDetail, component, divV, tableFromMap} from '@datagrok-libraries/u2';
import './master-detail.css';

interface Compound {
  id: string;
  name: string;
  mw: number;
  stage: string;
}

export function compoundBrowser(compounds: Compound[]): Control {
  return component(() => divV([new MasterDetail<Compound>({
    items: compounds,
    keyOf: (c) => c.id,
    row: (c) => ({label: `${c.id} · ${c.name}`, hint: c.stage}),
    detail: (c) => tableFromMap({'Name': c.name, 'MW': `${c.mw.toFixed(2)} g/mol`, 'Stage': c.stage}),
    emptyText: 'No compounds to show',
    label: 'Compounds',
  })], 'compound-browser'), {name: 'compound-browser'});
}

// Live demo
const STAGES = ['Hit', 'Lead', 'Candidate'];

export function render(host: HTMLElement): Control {
  const browser = compoundBrowser(Array.from({length: 40}, (_, i) =>
    ({id: `CPD-${String(i + 1).padStart(4, '0')}`, name: `Compound ${i + 1}`, mw: 180 + (i * 37) % 420, stage: STAGES[i % 3]})));
  host.append(browser.root);
  return browser;
}
