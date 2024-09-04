import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

export enum SEARCH_MODE {
  EXACT = 'Exact',
  SIMILAR = 'Similar',
  SUBSTRUCTURE = 'Substructure'
}

enum CURRENCY {
  USD = 'USD',
  EUR = 'EUR',
}

export type EnamineStorePack = {
  'amount':  number,
  'measure': string,
  'price': number,
  'priceOld': number,
  'stockAvailability': number,
  'usCurK': number,
};

export type EnamineStoreSearchResult = {
  "code": string,
  "stock": number,
  "productType": string,
  "id": string,
  "cdId": number,
  "restriction": string,
  "restrictionErrorMsg": string,
  "formula": string,
  "smile": string,
  "name": string,
  "cas": string,
  "clogp": number,
  "mw": number,
  "tpsa": number,
  "hbd": number,
  "hba": number,
  "rotb": number,
  "fsp3": number,
  "purity": number,
  "MFCD": number,
  "transportConditions": number,
  "storageConditions": number,
  "solubilityLogarithm": number,
  "polarSurfaceArea": number,
  "hydrogenAcceptors": number,
  "hydrogenDonors": number,
  "collection": string,
  "mfcd": string,
};

export const searchModeToCommandMap = {
  [SEARCH_MODE.EXACT]: 'EXACT',
  [SEARCH_MODE.SIMILAR]: 'SIM',
  [SEARCH_MODE.SUBSTRUCTURE]: 'SUB',
};

type EnamineMolProperties =
  {'ID': string, 'Code': string, 'Formula': string, 'MW': number,/*  'Availability': number, 'Delivery': string */};

enum CATALOG_TYPE {
  SCR = 'SCR',
  REAL = 'REAL',
  BB = 'BB',
}

const WIDTH = 150;
const HEIGHT = 75;

//tags: app
//name: Enamine Store
export function enamineStoreApp(): void {
  const molecule = ui.input.molecule('', {value: 'c1ccccc1O'});
  const searchMode =
    ui.input.choice('Mode', {value: SEARCH_MODE.SIMILAR, items: Object.keys(searchModeToCommandMap)}) as DG.InputBase<SEARCH_MODE>;
  const currency = ui.input.choice('Currency', {value: CURRENCY.USD, items: Object.values(CURRENCY)}) as DG.InputBase<string>;
  const similarity = ui.input.choice('Similarity', {value: '0.8', items: ['0.2', '0.4', '0.6', '0.8']}) as DG.InputBase<string>;
  const catalog =
    ui.input.choice('Catalog', {value: CATALOG_TYPE.BB, items: Object.values(CATALOG_TYPE)}) as DG.InputBase<CATALOG_TYPE>;
  const filterForm = ui.form([molecule, searchMode, currency, similarity, catalog]);
  const filtersHost = ui.div([filterForm], 'enamine-store-controls,pure-form');

  const emptyTable = DG.DataFrame.create();
  const view = grok.shell.addTableView(emptyTable);
  view.name = 'Enamine Store';
  view.basePath = '';
  view.description = 'Enamine Store search viewer';
  view.root.className = 'grok-view grok-table-view enamine-store';

  function update(): void {
    ui.setUpdateIndicator(view.root, true);
    grok.data.callQuery('EnamineStore:Search', {
      'q': molecule.value,
      'cat': catalog.value,
      'currency': currency.value,
      'type': 'SMARTS',
      'sstype': searchModeToCommandMap[searchMode.value],
      'sim': parseFloat(similarity.value)
    }, true).then(async (fc) => {
      const data = JSON.parse(fc.getParamValue('stringResult'))['searchResults'] as EnamineStoreSearchResult[];
      //get packs information for each id
      const promises: Promise<DG.FuncCall>[] = [];
      data.forEach((it) => {
        const opts: { [key: string]: any } = {
          'id': it['code'],
          'cat': catalog.value,
          'cur': currency.value,
        };
        promises.push(grok.data.callQuery('EnamineStore:Price', opts, true))
      });
      const packs: { [key: string]: EnamineStorePack[] } = {};
      const packsResFuncCalls = await Promise.all(promises);
      packsResFuncCalls.forEach((fc) => {
        const res = JSON.parse(fc.getParamValue('stringResult'));
        packs[res['code']] = res['samples'] ?? [];
      });


      view.dataFrame = data !== null ? dataToTable(data, packs, 'enaminestore') : emptyTable;
      ui.setUpdateIndicator(view.root, false);
    });
  }

  update();

  molecule.onChanged.subscribe(() => update());
  searchMode.onChanged.subscribe((value) => {
    similarity.enabled = value === SEARCH_MODE.SIMILAR;
    update();
  });
  currency.onChanged.subscribe(() => update());
  similarity.onChanged.subscribe(() => update());
  catalog.onChanged.subscribe(() => update());

  const acc = view.toolboxPage.accordion;
  acc.addPane('Enamine Store', () => filtersHost, true, acc.panes[0]);
}

//name: Databases | Enamine Store
//description: Enamine Store Samples
//tags: panel, widgets
//input: string smiles {semType: Molecule}
//output: widget result
//condition: true
export function enamineStorePanel(smiles: string): DG.Widget {
  if (DG.chem.isMolBlock(smiles))
    smiles = DG.chem.convert(smiles, DG.chem.Notation.MolBlock, DG.chem.Notation.Smiles);
  const acc = ui.accordion();
  const catalogToData: {[catalogType in CATALOG_TYPE]?: {[searchMode in SEARCH_MODE]?: HTMLDivElement}} = {};
  const catalog = ui.input.choice('Catalog', {value: CATALOG_TYPE.SCR, items: Object.values(CATALOG_TYPE),
    onValueChanged: (inp, value) => {
      const exactPanel = acc.getPane(SEARCH_MODE.EXACT);
      const similarPanel = acc.getPane(SEARCH_MODE.SIMILAR);
      const substructurePanel = acc.getPane(SEARCH_MODE.SUBSTRUCTURE);
      const exactExpanded = exactPanel?.expanded ?? false;
      const similarExpanded = similarPanel?.expanded ?? false;
      const substructureExpanded = substructurePanel?.expanded ?? false;
      for (const pane of acc.panes)
        acc.removePane(pane);

      acc.addPane(SEARCH_MODE.EXACT, () => {
        catalogToData[value] ??= {};
        catalogToData[value]![SEARCH_MODE.EXACT] ??= createSearchPanel(SEARCH_MODE.EXACT, smiles, value);
        return catalogToData[value]![SEARCH_MODE.EXACT]!;
      }, exactExpanded);
      acc.addPane(SEARCH_MODE.SIMILAR, () => {
        catalogToData[value] ??= {};
        catalogToData[value]![SEARCH_MODE.SIMILAR] ??=
          createSearchPanel(SEARCH_MODE.SIMILAR, smiles, value);
        return catalogToData[value]![SEARCH_MODE.SIMILAR]!;
      }, similarExpanded);
      acc.addPane(SEARCH_MODE.SUBSTRUCTURE, () => {
        catalogToData[value] ??= {};
        catalogToData[value]![SEARCH_MODE.SUBSTRUCTURE] ??=
          createSearchPanel(SEARCH_MODE.SUBSTRUCTURE, smiles, value);
        return catalogToData[value]![SEARCH_MODE.SUBSTRUCTURE]!;
      }, substructureExpanded);
    }}) as DG.InputBase<CATALOG_TYPE>;
  catalog.fireChanged();

  const form = ui.form([catalog]);
  const panels = ui.divV([form, acc.root]);

  return DG.Widget.fromRoot(panels);
}

//description: Creates search panel
export function createSearchPanel(searchMode: SEARCH_MODE, smiles: string, catalog: CATALOG_TYPE = CATALOG_TYPE.BB,
): HTMLDivElement {
  const currency = CURRENCY.USD;
  const headerHost = ui.divH([/*ui.h2(searchMode)*/], 'enamine-store-panel-header');
  const compsHost = ui.div([ui.loader()], 'd4-flex-wrap chem-viewer-grid');
  const panel = ui.divV([headerHost, compsHost], 'enamine-store-panel');
  const options: {[key: string]: any} = {
    'q': smiles,
    'cat': catalog,
    'currency': currency,
    'type': 'SMARTS',
    'sstype': searchModeToCommandMap[searchMode]
  };
  grok.data.callQuery('EnamineStore:Search', options, true).then(async (fc) => {
    compsHost.firstChild?.remove();
    const data = JSON.parse(fc.getParamValue('stringResult'))['searchResults'] as EnamineStoreSearchResult[];
    if (data === null) {
      compsHost.appendChild(ui.divText('No matches'));
      return;
    }
    //get packs information for each id
    const promises: Promise<DG.FuncCall>[] = [];
    data.forEach((it)=> {
      const opts: {[key: string]: any} = {
        'id': it['code'],
        'cat': catalog,
        'cur': currency,
      };
      promises.push(grok.data.callQuery('EnamineStore:Price', opts, true))
    });
    const packs: {[key: string]: EnamineStorePack[]} = {};
    const packsResFuncCalls = await Promise.all(promises);
    packsResFuncCalls.forEach((fc) => {
      const res = JSON.parse(fc.getParamValue('stringResult'));
      packs[res['code']] = res['samples'] ?? [];
    });


    let similarityResult: DG.DataFrame | null = null;
    if (searchMode === SEARCH_MODE.SIMILAR) {
      const df = DG.DataFrame.create(data.length);
      const smCol = df.columns.addNewString('smiles').init((i) => data[i]['smile']);
      similarityResult = await grok.chem.findSimilar(smCol, smiles, {limit: 20, cutoff: 0.01});
    }

    if (similarityResult?.rowCount === 0) {
      compsHost.appendChild(ui.divText('No matches'));
      return;
    }

    for (let i = 0; i < Math.min(similarityResult?.rowCount ?? data.length, 20); i++) {
      const idx = searchMode === SEARCH_MODE.SIMILAR ? similarityResult?.get('index', i) : i;
      const comp = data[idx];
      const currentSmiles = comp['smile'];
      const molHost = ui.divV([]);
      grok.functions.call('Chem:drawMolecule', {'molStr': currentSmiles, 'w': WIDTH, 'h': HEIGHT, 'popupMenu': true})
        .then((res: HTMLElement) => {
          molHost.append(res);
          if (searchMode === SEARCH_MODE.SIMILAR)
            molHost.appendChild(ui.divText(`Score: ${similarityResult?.get('score', i).toFixed(2)}`));
        });

      const id = comp['id'];
      const code = comp['code'];
      const props: EnamineMolProperties = {
        'ID': id,
        'Code': code,
        'Formula': comp['formula'],
        'MW': comp['mw'],
      };

      for (const pack of packs[code]) {
        //@ts-ignore: idk how to properly define type of props so that it has both required fields and optional any field
        props[`${pack['amount']} ${pack['measure']}`] = `${pack['price']} ${currency}`;
      }
      ui.tooltip.bind(molHost, () => ui.divV([ui.tableFromMap(props), ui.divText('Click to open in the store.')]));
      
    //product URL is unavailable in new Enamine API version
     // molHost.addEventListener('click', () => window.open(comp['productUrl'], '_blank'));
      compsHost.appendChild(molHost);
    }
    headerHost.appendChild(ui.iconFA('arrow-square-down', () =>
      grok.shell.addTableView(dataToTable(data, packs, `EnamineStore ${searchMode}`)), 'Open compounds as table'));
    compsHost.style.overflowY = 'auto';
  }).catch((err) => {
    compsHost.firstChild?.remove();
    const div = ui.divText('No matches');
    ui.tooltip.bind(div, `${err}`);
    compsHost.appendChild(div);
  });
  return panel;
}

// description: Converts JSON data into DataFrame
function dataToTable(data: EnamineStoreSearchResult[], packsList: {[key: string]: EnamineStorePack[]},name: string): DG.DataFrame {
  const columns = [
    DG.Column.fromStrings('smiles', data.map((comp) => comp['smile'])),
    DG.Column.fromStrings('ID', data.map((comp) => comp['id'])),
    DG.Column.fromStrings('Formula', data.map((comp) => comp['formula'])),
    DG.Column.fromFloat32Array('MW', new Float32Array(data.map((comp) => comp['mw']))),
  // availability and deliveryDays are not available in new Enamine API
  //  DG.Column.fromInt32Array('Availability', new Int32Array(data.map((comp) => comp['availability']))),
  //  DG.Column.fromStrings('Delivery', data.map((comp) => comp['deliveryDays'])),
  ];
  let currency = null;
  const packsArrays = new Map();
  for (let searchResultIdx = 0; searchResultIdx < data.length; searchResultIdx++) {
    const packs = packsList[data[searchResultIdx]['code']];
    for (let packIdx = 0; packIdx < packs.length; packIdx++) {
      const pack = packs[packIdx];
      const packName = `${pack['amount']} ${pack['measure']}`;
      if (!packsArrays.has(packName))
        packsArrays.set(packName, new Float32Array(data.length));
      packsArrays.get(packName)[searchResultIdx] = pack['price'];
/*       if (currency === null && pack['currencyName'] !== null)
        currency = pack['currencyName']; */
    }
  }
  for (const packKey of packsArrays.keys()) {
    const column = DG.Column.fromFloat32Array(packKey, packsArrays.get(packKey));
    column.semType = 'Money';
    column.meta.format = `money(${currency === CURRENCY.USD ? '$' : '€'})`;
    columns.push(column);
  }
  const table = DG.DataFrame.fromColumns(columns);
  table.name = name;
  return table;
}
