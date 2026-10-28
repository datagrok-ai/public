import * as ui from 'datagrok-api/ui';

interface Campaign {
  code: string;
  name: string;
  template: string;
  status: 'In progress' | 'Submitted' | 'Archived';
  molecules: number;
  created: Date;
  author: string;
}

const COLUMNS: [string, keyof Campaign][] = [['Name', 'name'], ['Code', 'code'], ['Created', 'created'],
  ['Molecules', 'molecules'], ['Status', 'status']];

export function campaignCatalog(campaigns: Campaign[], open: (c: Campaign) => void): HTMLElement {
  let items = campaigns;
  let sortBy: keyof Campaign | null = null;
  let ascending = true;
  const host = ui.div();
  const remove = (c: Campaign) => ui.dialog('Delete campaign')
    .add(ui.divText(`Delete campaign ${c.name} (${c.code})?`))
    .onOK(() => {
      items = items.filter((x) => x !== c);
      rebuild();
    })
    .show();
  const compare = (a: Campaign, b: Campaign) => sortBy === null ? b.created.getTime() - a.created.getTime() :
    (a[sortBy] < b[sortBy] ? -1 : a[sortBy] > b[sortBy] ? 1 : 0) * (ascending ? 1 : -1);
  const rebuild = () => {
    const table = ui.table([...items].sort(compare), (c) => [ui.link(c.name, () => open(c), 'Continue the campaign'),
      c.code, c.created.toLocaleDateString(), c.molecules, c.status,
      ui.iconFA('trash-alt', () => remove(c), 'Delete campaign')], [...COLUMNS.map(([header]) => header), '']);
    table.querySelectorAll('th').forEach((th, i) => {
      const column = COLUMNS[i];
      if (column === undefined)
        return;
      th.textContent = column[0] + (sortBy === column[1] ? (ascending ? ' ▲' : ' ▼') : '');
      th.style.cursor = 'pointer';
      th.addEventListener('click', () => {
        ascending = sortBy === column[1] ? !ascending : true;
        sortBy = column[1];
        rebuild();
      });
    });
    ui.empty(host);
    host.append(table, ui.divText(`${items.length} campaigns`));
  };
  rebuild();
  return host;
}
