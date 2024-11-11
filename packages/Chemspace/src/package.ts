import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { COUNTRY_CODES } from './country-codes';
import { ChemspacePricesTableItem, ChemspaceResult } from './model';

const host = 'https://api.chem-space.com';
let token: string | null = null;
export const _package = new DG.Package();

const WIDTH = 150;
const HEIGHT = 75;

enum SEARCH_MODE {
  SIMILAR = 'Similar',
  SUBSTRUCTURE = 'Substructure',
  TEXT = 'Text'
}

enum CATEGORY {
  CSMS = 'CSMS',
  CSMB = 'CSMB', 
  CSCS = 'CSCS', 
  CSSB = 'CSSB', 
  CSSS = 'CSSS',
}

const modeToParam = {[SEARCH_MODE.SIMILAR]: 'sim', [SEARCH_MODE.SUBSTRUCTURE]: 'sub', [SEARCH_MODE.TEXT]: 'text'};

//tags: app
//name: Chemspace
export async function app(): Promise<void> {
  await getApiToken();

  const molecule = ui.input.molecule('', {value: 'c1ccccc1O'});
  const mode = ui.input.choice('Mode', {
    value: SEARCH_MODE.SIMILAR, 
    items: Object.values(SEARCH_MODE),
    onValueChanged: () => update(),
  }) as DG.InputBase<SEARCH_MODE>;
  const shipToCountry = ui.input.choice('Ship to country', {
    value: 'United States',
    items: Object.keys(COUNTRY_CODES),
    onValueChanged: () => update(),
  }) as DG.InputBase<string>;
  const category = ui.input.choice('Category', {
    value: CATEGORY.CSCS, items: Object.values(CATEGORY),
    onValueChanged: () => update(),
  }) as DG.InputBase<CATEGORY>;
  const filterForm = ui.form([mode, shipToCountry, category]);
  const filtersHost = ui.div([molecule, filterForm], 'chemspace-controls,ui-form');

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
      ui.setUpdateIndicator(view.root, false);
    }
    grok.functions.call(`${_package.name}:queryMultipart`, {
      path: `search/${modeToParam[mode.value]}`,
      formParamsStr: JSON.stringify({'SMILES': molecule.value}),
      paramsStr: JSON.stringify({category: category.value, shipToCountry: COUNTRY_CODES[shipToCountry.value! as keyof typeof COUNTRY_CODES]})
    }).then((res) => setDataFrame(DG.DataFrame.fromJson(res)))
      .catch((_) => setDataFrame(emptyTable));
  }

  update();

  molecule.onChanged.subscribe(() => update());
  mode.onChanged.subscribe(() => {
    update();
  });

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

  const updateSearchResults = (acc: DG.Accordion, categoryToData: {[key: string]: {[searchMode in SEARCH_MODE]?: HTMLDivElement}},
    cat: CATEGORY, shipToCountry: COUNTRY_CODES) => {
    const similarPanel = acc.getPane(SEARCH_MODE.SIMILAR);
    const substructurePanel = acc.getPane(SEARCH_MODE.SUBSTRUCTURE);
    const similarExpanded = similarPanel?.expanded ?? false;
    const substructureExpanded = substructurePanel?.expanded ?? false;
    const cacheKey = getCategoryCacheKey(cat, shipToCountry);
    for (const pane of acc.panes)
      acc.removePane(pane);
    acc.addPane(SEARCH_MODE.SIMILAR, () => {
      categoryToData[cacheKey] ??= {};
      categoryToData[cacheKey]![SEARCH_MODE.SIMILAR] ??= createSearchPanel(SEARCH_MODE.SIMILAR, smiles, cat, shipToCountry);
      return categoryToData[cacheKey]![SEARCH_MODE.SIMILAR]!;
    }, similarExpanded);
    acc.addPane(SEARCH_MODE.SUBSTRUCTURE, () => {
      categoryToData[cacheKey] ??= {};
      categoryToData[cacheKey]![SEARCH_MODE.SUBSTRUCTURE] ??= createSearchPanel(SEARCH_MODE.SUBSTRUCTURE, smiles, cat, shipToCountry);
      return categoryToData[cacheKey]![SEARCH_MODE.SUBSTRUCTURE]!;
    }, substructureExpanded);
  };

  let panels: HTMLDivElement | null = null;
  try {
    if (DG.chem.isMolBlock(smiles))
      smiles = DG.chem.convert(smiles, DG.chem.Notation.MolBlock, DG.chem.Notation.Smiles);
    await getApiToken();
    const acc = ui.accordion();
    const categoryToData: { [key: string]: { [searchMode in SEARCH_MODE]?: HTMLDivElement } } = {};
    const shipToCountry = ui.input.choice('Ship to country', {
      value: 'United States',
      items: Object.keys(COUNTRY_CODES),
      onValueChanged: (value) => updateSearchResults(acc, categoryToData, category.value, COUNTRY_CODES[value as keyof typeof COUNTRY_CODES]),
    }) as DG.InputBase<string>;
    const category = ui.input.choice('Category', {
      value: CATEGORY.CSCS, items: Object.values(CATEGORY),
      onValueChanged: (value) => updateSearchResults(acc, categoryToData, value, COUNTRY_CODES[shipToCountry.value! as keyof typeof COUNTRY_CODES]),
    }) as DG.InputBase<CATEGORY>;
    category.fireChanged();
    panels = ui.divV([ui.form([shipToCountry, category]), acc.root]);
  } catch (e: any) {
    panels = ui.divText(e.message);
  }
  return new DG.Widget(panels);
}

//description: Creates search panel
export function createSearchPanel(searchMode: SEARCH_MODE, smiles: string, category: CATEGORY = CATEGORY.CSCS,
shipToCountry: COUNTRY_CODES = COUNTRY_CODES['United States']): HTMLDivElement {
  const headerHost = ui.divH([/*ui.h2(panelName)*/], 'chemspace-panel-header');
  const compsHost = ui.div([ui.loader()], 'd4-flex-wrap chem-viewer-grid');
  const panel = ui.divV([headerHost, compsHost], 'chemspace-panel');

  const queryParams: {[key: string]: any} = {
    'shipToCountry': shipToCountry,
    'categories': category
  };
  grok.functions.call(`${_package.name}:queryMultipart`, {
    path: `search/${modeToParam[searchMode]}`,
    formParamsStr: JSON.stringify({'SMILES': smiles}),
    paramsStr: JSON.stringify(queryParams)
  }).then(async (resStr: string) => {
    const res: ChemspaceResult[] = JSON.parse(resStr);
      compsHost.firstChild?.remove();
      if (res.length === 0) {
        compsHost.appendChild(ui.divText('No matches'));
        return;
      }

      function getTooltip(idx: number): HTMLDivElement {
        const props: {[key: string]: any} = {
          'ID': res[idx].csId,
          'Formula': res[idx].molFormula,
        }
        Object.keys(res[idx].properties).forEach((prop) => props[prop] = (res[idx].properties as any)[prop])
        return ui.divV([ui.tableFromMap(props), ui.divText('Click to open in the store.')]);
      }

      const smilesCol: DG.Column<string> = DG.Column.string('smiles', res.length).init((i) => res[i].smiles);

      let similarityResult: DG.DataFrame | null = null;
      if (searchMode === SEARCH_MODE.SIMILAR)
        similarityResult = await grok.chem.findSimilar(smilesCol, smiles, {limit: 20, cutoff: 0.1});

      for (let i = 0; i < Math.min((similarityResult ? similarityResult.rowCount : res.length), 20); i++) {
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
        molHost.addEventListener('click', () => window.open(res[idx].link, '_blank'));
        compsHost.appendChild(molHost);
      }
      headerHost.appendChild(ui.iconFA('arrow-square-down', () => {
        const df = DG.DataFrame.fromJson(JSON.stringify(res))
        df.name = `Chemspace ${searchMode}`;
        grok.shell.addTableView(df);
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
export async function pricesPanel(id: string): Promise<DG.Widget> {
  let prices: HTMLDivElement | null = null;
  const resData = ui.div([ui.loader()]);
  const updatePrices = (categoryToData: { [key: string]: HTMLDivElement }, cat: CATEGORY, shipToCountry: COUNTRY_CODES) => {
    ui.empty(resData);
    const cacheKey = getCategoryCacheKey(cat, shipToCountry);
    if (!categoryToData[cacheKey]) {
      resData.append(ui.loader());
      grok.functions.call(`${_package.name}:queryMultipart`, {
        path: `search/${modeToParam[SEARCH_MODE.TEXT]}`,
        formParamsStr: JSON.stringify({'query': id}),
        paramsStr: JSON.stringify({'shipToCountry': shipToCountry, 'categories': cat})
      }).then(async (resStr: string) => {
          const res: ChemspaceResult[] = JSON.parse(resStr);
          ui.empty(resData);
          if (res.length === 0) {
            resData.appendChild(ui.divText('No matches'));
            return;
          };
          if (res.length > 1) {
            resData.appendChild(ui.divText(`Ambigous results for id ${id}`));
            return;
          };
          const offers = res[0].offers;
          const chemspacePricesArray: ChemspacePricesTableItem[] = [];
          offers.forEach((offer) => {
            for (let i = 0; i < offer.prices.length; i++)
              chemspacePricesArray.push({
                packMg: offer.prices[i].packMg,
                priceUsd: offer.prices[i].priceUsd,
                priceEur: offer.prices[i].priceEur,
                vendorName: offer.vendorName,
                vendorCode: offer.vendorCode,
                leadTimeDays: offer.leadTimeDays,
                purity: offer.purity
            });
          });
          const table = DG.DataFrame.fromJson(JSON.stringify(chemspacePricesArray));
          //setting format for money columns
          ['priceUsd', 'priceEur'].forEach((name: string) => {
            table.col(name)!.semType = 'Money';
            if(name === 'priceUsd')
              table.col(name)!.meta.format = 'money($)';
          })
          const grid = DG.Grid.create(table);
          grid.root.classList.add('chemspace-prices-grid');
          const button = ui.bigButton('ORDER', () => window.open(res[0].link, '_blank'));
          button.classList.add('chemspace-order-button');
          resData.append(grid.root);
          resData.append(button);
      })
      .catch((err: any) => {
        ui.empty(resData);
        resData.append(ui.divText(err.message));
      });
    }
    else
      resData.append(categoryToData[cacheKey]);
  }
  try {
    await getApiToken();
    const categoryToData: { [key: string]: HTMLDivElement } = {};
    const shipToCountry = ui.input.choice('Ship to country', {
      value: 'United States',
      items: Object.keys(COUNTRY_CODES),
      onValueChanged: (value) => updatePrices(categoryToData, category.value, COUNTRY_CODES[value as keyof typeof COUNTRY_CODES]),
    }) as DG.InputBase<string>;
    const category = ui.input.choice('Category', {
      value: CATEGORY.CSCS, items: Object.values(CATEGORY),
      onValueChanged: (value) => updatePrices(categoryToData, value, COUNTRY_CODES[shipToCountry.value! as keyof typeof COUNTRY_CODES]),
    }) as DG.InputBase<CATEGORY>;
    category.fireChanged();
    prices = ui.divV([ui.form([shipToCountry, category]), resData]);
  } catch (e: any) {
    prices = ui.divText(e.message);
  }
  return new DG.Widget(prices);
}

//description: Gets access token
async function getApiToken(): Promise<void> {
  if (token === null) {
    const t = await grok.data.query('Chemspace:AuthToken', null, true);
    token = t.get('access_token', 0) as string;
  }
}

//description: Perform query with multipart form data
//meta.cache: all
//meta.cache.invalidateOn: 0 0 1 * *
//input: string path
//input: string formParamsStr
//input: string paramsStr {optional: true}
//output: string result
export function queryMultipart(path: string, formParamsStr: string, paramsStr?: string): Promise<string> {
  const formParams: {[key: string]: string} = JSON.parse(formParamsStr);
  const params: {[key: string]: string | number} | undefined = paramsStr ? JSON.parse(paramsStr) : undefined;
  // TODO: Deprecate after WebQuery 'multipart/form-data' support
  return new Promise(function(resolve, reject) {
    const xhr = new XMLHttpRequest();
    const formData = new FormData();
    Object.keys(formParams).forEach((key) => formData.append(key, formParams[key]));
   // formData.append('SMILES', smiles);
    const queryParams = params ? `?${Object.keys(params).map((key) => `${key}=${params[key]}`).join('&')}` : '';
    xhr.open('POST', `${host}/v4/${path}${queryParams}`);
    xhr.setRequestHeader('Authorization', `Bearer ${token}`);
    xhr.setRequestHeader('Accept', 'application/json; version=2.6');
    xhr.onload = function(): void {
      if (this.status >= 200 && this.status < 300) {
        const list: ChemspaceResult[] = JSON.parse(xhr.responseText)['items'];
        if (list.length > 0)
          resolve(JSON.stringify(list));
        else
          reject(new Error('No matches'));
      } else
        reject(new Error(xhr.statusText));
    };
    xhr.onerror = (): void => reject(new Error(xhr.statusText));
    xhr.send(formData);
  });
}

function getCategoryCacheKey(cat: CATEGORY, shipToCountry: COUNTRY_CODES): string {
  return `${cat}|${shipToCountry}`;
}
