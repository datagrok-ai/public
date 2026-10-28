import {Control, DataTable, Format, component, computed, notify, signal} from '@datagrok-libraries/u2';
import './data-table.css';

interface Measurement {
  id: string;
  compound: string;
  assay: string;
  ic50: number;
  measured: Date;
}

export function measurementsTable(rows: Measurement[]): DataTable<Measurement> {
  const sort = signal('!ic50');
  const sorted = computed(() => {
    const key = sort.value.replace(/^!/, '') as keyof Measurement;
    const direction = sort.value.startsWith('!') ? -1 : 1;
    return sort.value === '' ? rows : [...rows].sort((a, b) => (a[key] < b[key] ? -1 : a[key] > b[key] ? 1 : 0) * direction);
  });
  return new DataTable<Measurement>({items: sorted, sort, keyOf: (m) => m.id,
    columns: [
      {name: 'id', header: 'Id', width: 'minmax(70px, 1fr)'},
      {name: 'compound', header: 'Compound', width: 'minmax(90px, 1fr)'},
      {name: 'assay', header: 'Assay', width: 'minmax(110px, 2fr)'},
      {name: 'ic50', header: 'IC50, nM', align: 'right', width: '90px', render: (m) => m.ic50.toFixed(1)},
      {name: 'measured', header: 'Measured', width: '100px', render: (m) => Format.date(m.measured)},
    ],
    onActivate: (m) => notify.info(`${m.compound}: ${m.ic50.toFixed(1)} nM in ${m.assay}`),
  });
}

// Live demo
export function render(host: HTMLElement): Control {
  const demo = component(() => {
    const table = measurementsTable(Array.from({length: 10000}, (_, i) => ({id: `M-${i + 1}`,
      compound: `CPD-${String(1 + (i * 7919) % 2500).padStart(4, '0')}`, assay: ['Kinase', 'Cell viability', 'hERG'][i % 3],
      ic50: ((i * 2654435761) % 100000) / 10, measured: new Date(2026, i % 10, 1 + i % 28)})));
    table.root.classList.add('measurements');
    return table;
  });
  host.append(demo.root);
  return demo;
}
