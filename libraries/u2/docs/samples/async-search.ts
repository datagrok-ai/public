import {AsyncView, BasicTable, Control, Form, TextInput, component} from '@datagrok-libraries/u2';

interface Compound {
  id: string;
  name: string;
  mw: number;
}

export function compoundSearch(search: (query: string, abort: AbortSignal) => Promise<Compound[]>): Control {
  return component(() => {
    const results = AsyncView.owned(search, (items) => new BasicTable<Compound>({items, columns: [
      {header: 'Id', render: (c) => c.id},
      {header: 'Name', render: (c) => c.name},
      {header: 'MW', render: (c) => c.mw.toFixed(1), align: 'right'},
    ]}).root, {empty: 'No compounds match', skeleton: true});
    const query = new TextInput({label: 'Search', search: true, placeholder: 'Name or id', onChanged: (q) => results.refresh(q)});
    results.refresh('');
    return [new Form().add(query), results];
  }, {name: 'compound-search'});
}

// Live demo
const NAMES = ['Aspirin', 'Caffeine', 'Ibuprofen', 'Naproxen', 'Paracetamol', 'Imatinib', 'Gefitinib', 'Erlotinib'];

async function searchCompounds(query: string, abort: AbortSignal): Promise<Compound[]> {
  await new Promise((resolve) => setTimeout(resolve, 500));
  abort.throwIfAborted();
  if (query.trim().toLowerCase() === 'error')
    throw new Error('The registry did not answer (503)');
  return NAMES.map((name, i) => ({id: `CHEMBL${(i + 1) * 25}`, name, mw: 150 + i * 47.3}))
    .filter((c) => `${c.id} ${c.name}`.toLowerCase().includes(query.trim().toLowerCase()));
}

export function render(host: HTMLElement): Control {
  const demo = compoundSearch(searchCompounds);
  host.append(demo.root);
  return demo;
}
