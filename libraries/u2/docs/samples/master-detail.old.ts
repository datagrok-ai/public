import * as ui from 'datagrok-api/ui';

interface Compound {
  id: string;
  name: string;
  mw: number;
  stage: string;
}

export function compoundBrowser(compounds: Compound[]): HTMLElement {
  const details = ui.divV([ui.divText('Select a compound to see its details')], 'compound-details');
  let selected: HTMLElement | null = null;
  const rows = compounds.map((c) => {
    const row = ui.divText(`${c.id} · ${c.name}`, 'compound-row');
    row.addEventListener('click', () => {
      selected?.classList.remove('compound-row-selected');
      row.classList.add('compound-row-selected');
      selected = row;
      ui.empty(details);
      details.append(ui.tableFromMap({'Name': c.name, 'MW': `${c.mw.toFixed(2)} g/mol`, 'Stage': c.stage}));
    });
    return row;
  });
  return ui.splitH([ui.divV(rows, 'compound-list'), details], null, true);
}
