import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

const host = 'https://api.chem-space.com';
let token: string | null = null;

enum SEARCH_MODE {
  SIMILAR = 'Similar',
  SUBSTRUCTURE = 'Substructure',
}

enum CATALOG_TYPE {
  BB = 'BB',
  SCR = 'SCR',
  REAL = 'REAL',
}

//TODO: Replace with acuatl type
type Price = any;

const searchModeToCommand = {'Similar': 'smarts/bb/sim', 'Substructure': 'smarts/bb/sub'};
const catalogToParam = {'BB': 'smarts/bb', 'SCR': 'smarts/sc', 'REAL': 'advanced/1'};
const modeToParam = {[SEARCH_MODE.SIMILAR]: 'sim', [SEARCH_MODE.SUBSTRUCTURE]: 'sub'};

//tags: app
//name: Chemspace
export async function app(): Promise<void> {
  const token = await getApiToken();

  const molecule = ui.moleculeInput('', 'c1ccccc1O');
  const mode = ui.choiceInput('Mode', SEARCH_MODE.SIMILAR, Object.values(SEARCH_MODE)) as DG.InputBase<SEARCH_MODE>;
  const similarity = ui.choiceInput('Similarity', '0.6', ['0.2', '0.4', '0.6', '0.8']) as DG.InputBase<string>;
  const catalog =
    ui.choiceInput('Catalog', CATALOG_TYPE.SCR, Object.values(CATALOG_TYPE)) as DG.InputBase<CATALOG_TYPE>;
  const filterForm = ui.form([molecule, mode, similarity, catalog]);
  const filtersHost = ui.div([filterForm], 'chemspace-controls,pure-form');

  const emptyTable = DG.DataFrame.create();
  const view = grok.shell.addTableView(emptyTable);
  view.name = 'Chemspace';
  view.basePath = '';
  view.description = 'Chemspace search viewer';
  view.root.className = 'grok-view grok-table-view chemspace';

  function update(): void {
    ui.setUpdateIndicator(view.root, true);

    function setDataFrame(t: DG.DataFrame): void {
      view.dataFrame = t;
      const order = ['smiles', 'CS-id', 'id', 'similarity', 'iupac_name', 'molformula', 'mw', 'molweight', 'cas', 'hac',
        'logp', 'rotb', 'hba', 'hbd', 'ring_count', 'fsp3', 'tpsa', 'mfcd', 'price_category', 'vendor_id', 'link'];
      if (t.rowCount > 0) {
        const names = t.columns.names();
        view.grid.columns.setOrder(order.filter((o) => names.includes(o)));
      }
      ui.setUpdateIndicator(view.root, false);
    }

    queryMultipart(`search/${catalogToParam[catalog.value]}/${modeToParam[mode.value]}`, molecule.value,
      mode.value === SEARCH_MODE.SIMILAR ? {'simThreshold': parseFloat(similarity.value) * 100} : null, token)
      .then((t) => setDataFrame(t))
      .catch((_) => setDataFrame(emptyTable));
  }

  update();

  molecule.onChanged(() => update());
  mode.onChanged(() => {
    similarity.enabled = mode.value === SEARCH_MODE.SIMILAR;
    update();
  });
  similarity.onChanged(() => update());
  catalog.onChanged(() => update());

  const acc = view.toolboxPage.accordion;
  acc.addPane('Chemspace', () => filtersHost, true, acc.panes[0]);
}

//name: Chemspace Samples
//description: Chemspace Samples
//tags: panel, widgets
//input: string smiles {semType: Molecule}
//output: widget result
//condition: true
export function samplesPanel(smiles: string): DG.Widget {
  const panels = ui.div();
  getApiToken().then(() => {
    panels.appendChild(createSearchPanel(SEARCH_MODE.SIMILAR, smiles));
    panels.appendChild(createSearchPanel(SEARCH_MODE.SUBSTRUCTURE, smiles));
  });
  return new DG.Widget(panels);
}

//description: Creates search panel
export function createSearchPanel(panelName: SEARCH_MODE, smiles: string): HTMLDivElement {
  const headerHost = ui.divH([ui.h2(panelName)], 'chemspace-panel-header');
  const compsHost = ui.div([ui.loader()], 'd4-flex-wrap');
  const panel = ui.divV([headerHost, compsHost], 'chemspace-panel');

  queryMultipart(`search/${searchModeToCommand[panelName]}`, smiles,
    panelName === SEARCH_MODE.SIMILAR ? {'simThreshold': 40} : null, token!)
    .then((t) => {
      compsHost.firstChild?.remove();
      if (t.rowCount === 0) {
        compsHost.appendChild(ui.divText('No matches'));
        return;
      }

      function getTooltip(idx: number): HTMLDivElement {
        const props = {
          'ID': t.get('CS-id', idx),
          'IUPAC': t.get('iupac_name', idx),
          'Formula': t.get('molformula', idx),
          'MW': t.get('molweight', idx),
        };
        return ui.divV([ui.tableFromMap(props), ui.divText('Click to open in the store.')]);
      }

      for (let n = 0; n < Math.min(t.rowCount, 20); n++) {
        const smiles = t.get('smiles', n);
        const mol = grok.chem.svgMol(smiles, 150, 75);
        ui.tooltip.bind(mol, () => getTooltip(n));
        mol.addEventListener('click', function() {
          window.open(t.get('link', n), '_blank');
        });
        compsHost.appendChild(mol);
      }
      headerHost.appendChild(ui.iconFA('arrow-square-down', () => {
        t.name = `Chemspace ${panelName}`;
        grok.shell.addTableView(t);
      }, 'Open compounds as table'));
      compsHost.style.overflowY = 'auto';
    })
    .catch((err) => {
      compsHost.firstChild?.remove();
      const div = ui.divText('No matches');
      ui.tooltip.bind(div, `${err}`);
      compsHost.appendChild(div);
    });

  return panel;
}

//name: Chemspace Prices
//description: Chemspace Prices
//tags: panel, widgets
//input: string id {semType: chemspace-id}
//output: widget result
//condition: true
export function pricesPanel(id: string): DG.Widget {
  const panel = ui.div([ui.loader()]);
  getApiToken().then(() => {
    const xhr = new XMLHttpRequest();
    xhr.open('GET', `${host}/v2/cs-id/${id}/prices`);
    xhr.setRequestHeader('Authorization', `Bearer ${token}`);
    xhr.setRequestHeader('Accept', 'application/json; version=2.6');

    function onError(): void {
      while (panel.firstChild)
        panel.firstChild.remove();
      panel.appendChild(ui.divText('Not found'));
    }

    xhr.onload = function(): void {
      if (this.status >= 200 && this.status < 300) {
        const map = JSON.parse(xhr.responseText);
        const t = pricesDataToTable(map['items']);
        const grid = DG.Grid.create(t);
        grid.root.style.width = '400px';
        grid.root.style.height = '300px';
        while (panel.firstChild)
          panel.removeChild(panel.firstChild);
        panel.appendChild(ui.div([grid.root]));
        const button = ui.bigButton('ORDER', () => window.open(map['link'], '_blank'));
        button.style.marginTop = '6px';
        panel.appendChild(button);
      } else
        onError();
    };
    xhr.onerror = (): void => onError();
    xhr.send();
  });

  return new DG.Widget(panel);
}

// description: Converts prices JSON items into DataFrame
function pricesDataToTable(items: Price[]): DG.DataFrame {
  const table = DG.DataFrame.fromJson(JSON.stringify(items));
  table.columns.remove('vendor_code');
  const packsArrays = new Map();
  for (let n = 0; n < items.length; n++) {
    const packs = items[n]['prices'];
    for (let m = 0; m < packs.length; m++) {
      const pack = packs[m];
      const name = `${pack['pack_g']} g`;
      if (!packsArrays.has(name))
        packsArrays.set(name, new Array(items.length));
      packsArrays.get(name)[n] = pack['price_usd'];
    }
  }
  for (const name of Array.from(packsArrays.keys()).sort()) {
    const column = DG.Column.fromList(DG.TYPE.FLOAT, name, packsArrays.get(name));
    column.semType = 'Money';
    column.setTag('format', 'money($)');
    table.columns.add(column);
  }
  return table;
}

//description: Gets access token
async function getApiToken(): Promise<string> {
  if (token === null) {
    const t = await grok.data.query('Chemspace:AuthToken', null, true, 100);
    token = t.get('access_token', 0) as string;
  }
  return token;
}

//description: Perform query with multipart form data
function queryMultipart(path: string, smiles: string, params: {[key: string]: string | number} | null, token: string,
): Promise<DG.DataFrame> {
  // TODO: Deprecate after WebQuery 'multipart/form-data' support
  return new Promise(function(resolve, reject) {
    const xhr = new XMLHttpRequest();
    const formData = new FormData();
    formData.append('SMILES', smiles);
    const queryParams =
      params !== null ? `?${Object.keys(params).map((key) => `${key}=${params[key]}`).join('&')}` : '';
    xhr.open('POST', `${host}/v2/${path}${queryParams}`);
    xhr.setRequestHeader('Authorization', `Bearer ${token}`);
    xhr.setRequestHeader('Accept', 'application/json; version=2.6');
    xhr.onload = function(): void {
      if (this.status >= 200 && this.status < 300) {
        const list = JSON.parse(xhr.responseText)['items'];
        if (list.length > 0)
          resolve(DG.DataFrame.fromJson(JSON.stringify(list)));
        else
          reject(new Error('No matches'));
      } else
        reject(new Error(xhr.statusText));
    };
    xhr.onerror = (): void => reject(new Error(xhr.statusText));
    xhr.send(formData);
  });
}
