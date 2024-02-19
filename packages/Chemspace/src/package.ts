import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

const host = 'https://api.chem-space.com';
let token: string | null = null;

const WIDTH = 150;
const HEIGHT = 75;

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
  const filterForm = ui.form([mode, similarity, catalog]);
  const filtersHost = ui.div([molecule, filterForm], 'chemspace-controls,pure-form');

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

//name: Databases | Chemspace
//description: Chemspace Samples
//tags: panel, widgets
//input: string smiles {semType: Molecule}
//output: widget result
//condition: true
export async function samplesPanel(smiles: string): Promise<DG.Widget> {
  if (DG.chem.isMolBlock(smiles))
    smiles = DG.chem.convert(smiles, DG.chem.Notation.MolBlock, DG.chem.Notation.Smiles)
  await getApiToken();
  const acc = ui.accordion();
  const catalogToData: {[catalogType in CATALOG_TYPE]?: {[searchMode in SEARCH_MODE]?: HTMLDivElement}} = {};
  const catalog = ui.choiceInput('Catalog', CATALOG_TYPE.SCR, Object.values(CATALOG_TYPE), () => {
    const similarPanel = acc.getPane(SEARCH_MODE.SIMILAR);
    const substructurePanel = acc.getPane(SEARCH_MODE.SUBSTRUCTURE);
    const similarExpanded = similarPanel?.expanded ?? false;
    const substructureExpanded = substructurePanel?.expanded ?? false;
    for (const pane of acc.panes)
      acc.removePane(pane);

    acc.addPane(SEARCH_MODE.SIMILAR, () => {
      catalogToData[catalog.value] ??= {};
      catalogToData[catalog.value]![SEARCH_MODE.SIMILAR] ??=
        createSearchPanel(SEARCH_MODE.SIMILAR, smiles, catalog.value);
      return catalogToData[catalog.value]![SEARCH_MODE.SIMILAR]!;
    }, similarExpanded);
    acc.addPane(SEARCH_MODE.SUBSTRUCTURE, () => {
      catalogToData[catalog.value] ??= {};
      catalogToData[catalog.value]![SEARCH_MODE.SUBSTRUCTURE] ??=
        createSearchPanel(SEARCH_MODE.SUBSTRUCTURE, smiles, catalog.value);
      return catalogToData[catalog.value]![SEARCH_MODE.SUBSTRUCTURE]!;
    }, substructureExpanded);
  }) as DG.InputBase<CATALOG_TYPE>;
  catalog.fireChanged();

  const form = ui.form([catalog]);
  const panels = ui.divV([form, acc.root]);
  return new DG.Widget(panels);
}

//description: Creates search panel
export function createSearchPanel(searchMode: SEARCH_MODE, smiles: string, catalog: CATALOG_TYPE = CATALOG_TYPE.BB,
): HTMLDivElement {
  const headerHost = ui.divH([/*ui.h2(panelName)*/], 'chemspace-panel-header');
  const compsHost = ui.div([ui.loader()], 'd4-flex-wrap chem-viewer-grid');
  const panel = ui.divV([headerHost, compsHost], 'chemspace-panel');

  queryMultipart(`search/${catalogToParam[catalog]}/${modeToParam[searchMode]}`, smiles,
    searchMode === SEARCH_MODE.SIMILAR ? {'simThreshold': 40} : null, token!)
    .then(async (table) => {
      compsHost.firstChild?.remove();
      if (table.rowCount === 0) {
        compsHost.appendChild(ui.divText('No matches'));
        return;
      }

      function getTooltip(idx: number): HTMLDivElement {
        const props = {
          'ID': table.get('CS-id', idx),
          'IUPAC': table.get('iupac_name', idx),
          'Formula': table.get('molformula', idx),
          'MW': table.get('molweight', idx),
        };
        return ui.divV([ui.tableFromMap(props), ui.divText('Click to open in the store.')]);
      }

      const smilesCol: DG.Column<string> = table.getCol('smiles');

      let similarityResult: DG.DataFrame | null = null;
      if (searchMode === SEARCH_MODE.SIMILAR) {
        similarityResult = await grok.chem.findSimilar(smilesCol, smiles, { limit: 20, cutoff: 0.8 });
      }

      for (let i = 0; i < Math.min((similarityResult ?? table).rowCount, 20); i++) {
        const idx = searchMode === SEARCH_MODE.SIMILAR ? similarityResult!.get('index', i) : i;
        const smiles = smilesCol.get(idx);
        const molHost = ui.div();
        grok.functions.call('Chem:drawMolecule', {'molStr': smiles, 'w': WIDTH, 'h': HEIGHT, 'popupMenu': true})
          .then((res: HTMLElement) => {
            molHost.append(res);
            if (searchMode === SEARCH_MODE.SIMILAR)
              molHost.appendChild(ui.divText(`Score: ${similarityResult?.get('score', i).toFixed(2)}`));
          });
        ui.tooltip.bind(molHost, () => getTooltip(idx));
        molHost.addEventListener('click', () => window.open(table.get('link', idx), '_blank'));
        compsHost.appendChild(molHost);
      }
      headerHost.appendChild(ui.iconFA('arrow-square-down', () => {
        table.name = `Chemspace ${searchMode}`;
        grok.shell.addTableView(table);
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
