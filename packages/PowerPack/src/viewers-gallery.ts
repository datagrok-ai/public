import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import $ from 'cash-dom';

let view: DG.TableView;
let table: DG.DataFrame;

const groupComparisons = [
  'Bar chart',
  'Radar',
  'Chord viewer',
  'Matrix plot',
  'Word cloud',
  'Word cloud viewer',
  'Pie chart',
  'Tree map',
];
const groupTrends = [
  'Line chart',
  'Histogram',
  'Box plot',
  'Sunburst',
];
const groupCorrelations = [
  'Scatter plot',
  'Density plot',
  '3d scatter plot',
  'PC plot',
  'Correlation plot',
  'Heat map',
];
const groupRelationships = [
  'Sankey',
  'Network diagram',
  'Tree',
  'Phylo tree',
  'Scaffold Tree',
  'PhylocanvasGL',
];
const groupMaps = [
  'Shape map',
  'Leaflet',
  'Globe',
  'Map',
];

interface ViewersTestData {
  checkFunc: ((table: DG.DataFrame) => boolean),
  tooltip: string,
}

const CHARTS_VIEWERS_TEST_DATA: {[key: string] : ViewersTestData} = {
  'Chord': {checkFunc: (table) => [...table.columns.categorical].length >= 2 && [...table.columns.numerical].length >= 1,
    tooltip: 'Chord viewer needs at least 2 categorical columns and 1 numerical column'},
  'Globe': {checkFunc: (table) => table.columns.toList().filter((col) => ['double', 'int']
    .includes(col.type)).length >= 1, tooltip: 'Globe viewer needs at least 1 numerical column'},
  'Group Analysis': {checkFunc: (table) => table.columns.length >= 1,
    tooltip: 'Group Analysis viewer needs at least 1 column'},
  'Radar': {checkFunc: (table) => table.columns.toList().filter((col) => ['double', 'int']
    .includes(col.type)).length >= 1, tooltip: 'Radar viewer needs at least 1 numerical column'},
  'Sankey': {checkFunc: (table) => table.columns.toList().filter((col) => col.type === 'string' &&
    col.categories.length <= 50).length >= 2 && table.columns.toList().filter((col) => ['double', 'int'].includes(col.type)).length >= 1,
    tooltip: 'Sankey viewer needs at least 2 string columns with less than 50 categories and 1 numerical column'},
  'Sunburst': {checkFunc: (table) => table.columns.toList().length >= 1,
    tooltip: 'Sunburst viewer needs at least 1 column'},
  'Surface plot': {checkFunc: (table) => table.columns.toList().length >= 3,
    tooltip: 'Surface plot viewer needs at least 3 columns'},
  'Timelines': {checkFunc: (table) => table.columns.toList().filter((col) => col.type === DG.COLUMN_TYPE.STRING)
    .length >= 1 && [...table.columns.numerical].length >= 1, tooltip: 'Timelines viewer needs at least 1 string column and 1 numerical column'},
  'Tree': {checkFunc: (table) => table.columns.toList().length >= 1, tooltip: 'Tree viewer needs at least 1 column'},
  'Word cloud': {checkFunc: (table) => table.columns.toList().filter((col) => col.type === DG.TYPE.STRING)
    .length >= 1, tooltip: 'Word cloud viewer needs at least 1 string column'},
};

const BIOCHEM_VIEWERS_TEST_DATA: {[key: string] : ViewersTestData} = {
  'Chem Similarity Search': {checkFunc: (table) => table.columns.bySemType(DG.SEMTYPE.MOLECULE) !== null,
    tooltip: 'Chem Similarity Search viewer needs at least 1 molecular column'},
  'Chem Diversity Search': {checkFunc: (table) => table.columns.bySemType(DG.SEMTYPE.MOLECULE) !== null,
    tooltip: 'Chem Diversity Search viewer needs at least 1 molecular column'},
  'Sequence Similarity Search': {checkFunc: (table) => table.columns.bySemType(DG.SEMTYPE.MACROMOLECULE) !== null,
    tooltip: 'Chem Similarity Search viewer needs at least 1 sequence column'},
  'Sequence Diversity Search': {checkFunc: (table) => table.columns.bySemType(DG.SEMTYPE.MACROMOLECULE) !== null,
    tooltip: 'Chem Diversity Search viewer needs at least 1 sequence column'},
};

const viewers: any = {};
const jsViewers: any = {};
const rootViewers = ui.divH([], 'viewer-gallery');
let dlg: DG.Dialog;
let search: DG.InputBase;
const viewersCount = ui.div([], 'vg-counter-label');
const recentViewersRoot = ui.divH([], 'viewer-gallery');
const recentLabel = ui.div(['Recently used'], 'vg-counter-label');
const recentBlock = ui.divV([recentLabel, recentViewersRoot], {style: {marginBottom: '10px'}});

export function viewersDialog(currentView: DG.TableView, currentTable: DG.DataFrame) {
  getViewers(viewers, currentTable);
  getJsViewers(jsViewers, currentTable);

  view = currentView;
  table = currentTable;

  dlg = ui.dialog('Add Viewer');
  $(dlg.root.lastChild).hide();

  dlg.root.addEventListener('keypress',(event)=>{
    if (event.key === 'Escape'){
      dlg.close();
    }
  })
  search = ui.input.search('', {value: '', onValueChanged: (value) => findViewer(value)});
  search.input.setAttribute('tabindex', '-1');
  search.input.setAttribute('placeholder', 'Search by name, keywords, description, tag, or package');

  var delta = 500;
  var lastKeypressTime = 0;

  search.input.onkeyup = (event) => {
    if (event.key === 'Escape'){
      search.fireChanged();
      var thisKeypressTime = Number(new Date());
      if ( thisKeypressTime - lastKeypressTime <= delta ) {
        if (search.value === ''){
          dlg.close();
        }
        thisKeypressTime = 0;
      }
      lastKeypressTime = thisKeypressTime;
    }
  };

  const searchIcon = ui.iconFA('search');
  searchIcon.classList.add('vg-search-icon');

  recentViewersRoot.innerHTML = '';
  const recentNames = getRecentViewersList();

  if (recentNames.length > 0) {
    const allViewersMap: {[name: string]: {idx: string, list: any}} = {};
    for (const i in viewers)
      allViewersMap[viewers[i].name] = {idx: i, list: viewers};
    for (const i in jsViewers)
      allViewersMap[jsViewers[i].name] = {idx: i, list: jsViewers};
    let hasValidRecents = false;
    for (const name of recentNames) {
      if (allViewersMap[name]) {
        recentViewersRoot.append(renderCard(allViewersMap[name].list, allViewersMap[name].idx));
        hasValidRecents = true;
      }
    }
    recentBlock.style.display = hasValidRecents ? 'flex' : 'none';
  }
  else
    recentBlock.style.display = 'none';

  rootViewers.innerHTML = '';

  for (const i in viewers)
    rootViewers.append(renderCard(viewers, i));


  for (const i in jsViewers)
    rootViewers.append(renderCard(jsViewers, i));

  const tags = ui.divV([
    ui.label('Relative tags:'),
    generateTags(),
  ], 'vg-tags');

  const root = ui.divH([
    ui.block([recentBlock, viewersCount, rootViewers], 'viewer-gallery-root'),
    tags,
  ]);

  dlg.add(ui.block([ui.div([searchIcon, search.input], 'd4-search-ba')], 'vg-controls grok-gallery-search-bar'));
  dlg.add(root);
  dlg.showModal(true);

  getTotalViewer();
  setTabIndex(rootViewers);
}


function getViewers(viewers: { [v: string]: { [k: string]: any } }, table: DG.DataFrame) {
  let viewerList = [];

  for (const value of Object.values(DG.CORE_VIEWER)) {
    if (value !== 'Shape Map')
      viewerList.push(value);
  }
  viewerList.push('Column viewer');
  viewerList.push('Web viewer');
  viewerList.push('Scripting viewer');

  viewerList = [...new Set(viewerList)];
  for (const i in viewerList) {
    const isViewerEnabledMsg = DG.Viewer.canVisualize(viewerList[i], table);
    Object.assign(viewers, {
      [i]: {
        name: viewerList[i],
        icon: 'grok-icon svg-icon svg-' + viewerList[i].toLowerCase().replace(/(\s)/g, '-'),
        enabled: isViewerEnabledMsg == null,
        tooltip: isViewerEnabledMsg == null ? '' : isViewerEnabledMsg,
        group: '',
        type: 'viewer',
      },
    });
    if (groupComparisons.includes(viewerList[i])) viewers[i]['group'] = 'Comparison';
    else if (groupCorrelations.includes(viewerList[i])) viewers[i]['group'] = 'Correlation';
    else if (groupRelationships.includes(viewerList[i])) viewers[i]['group'] = 'Relationship';
    else if (groupTrends.includes(viewerList[i])) viewers[i]['group'] = 'Trend';
    else if (groupMaps.includes(viewerList[i])) viewers[i]['group'] = 'Map';
    else viewers[i]['group'] = 'Misc';
  }
}

function getJsViewers(jsViewers: { [v: string]: { [k: string]: any } }, table: DG.DataFrame) {
  const skip = ['TestViewerForProperties', 'OutliersSelectionViewer'];
  const list = DG.Func.find({meta: {role: DG.FUNC_TYPES.VIEWER}}).filter((v) => !skip.includes(v.friendlyName));
  let i = 0;
  for (const v of list) {
    let isViewerEnabled = true;
    if (v.package.name === 'Charts') {
      const viewerTestData = CHARTS_VIEWERS_TEST_DATA[v.friendlyName];
      if (viewerTestData !== undefined) {
        if (!viewerTestData.checkFunc(table))
          isViewerEnabled = false;
      }
      else
        isViewerEnabled = false;
    } else if ((v.package.name === 'Chem' || v.package.name === 'Bio') && BIOCHEM_VIEWERS_TEST_DATA[v.friendlyName])
      isViewerEnabled = BIOCHEM_VIEWERS_TEST_DATA[v.friendlyName].checkFunc(table);
    else {
      if (v.options['showInGallery'] === 'false')
        continue;
    }
    Object.assign(jsViewers, {
      [i]: {
        name: v.friendlyName,
        enabled: isViewerEnabled,
        tooltip: isViewerEnabled ? '' : CHARTS_VIEWERS_TEST_DATA[v.friendlyName] ?
          CHARTS_VIEWERS_TEST_DATA[v.friendlyName].tooltip : BIOCHEM_VIEWERS_TEST_DATA[v.friendlyName] ?
          BIOCHEM_VIEWERS_TEST_DATA[v.friendlyName].tooltip : 'Viewer cannot be created from viewer gallery',
        icon: (v.options['icon'] != undefined) ? `${v.package.webRoot.endsWith('/') ?
          v.package.webRoot : v.package.webRoot + '/'}${v.options['icon']}` : 'svg-project',
        description: v.description,
        keywords: (v.options['keywords'] != null) ? v.options['keywords'] : '',
        group: '',
        package: v.package.name,
        type: 'js-viewer',
      },
    });
    if (groupComparisons.includes(v.friendlyName)) jsViewers[i]['group'] = 'Comparison';
    else if (groupCorrelations.includes(v.friendlyName)) jsViewers[i]['group'] = 'Correlation';
    else if (groupRelationships.includes(v.friendlyName)) jsViewers[i]['group'] = 'Relationship';
    else if (groupTrends.includes(v.friendlyName)) jsViewers[i]['group'] = 'Trend';
    else if (groupMaps.includes(v.friendlyName)) jsViewers[i]['group'] = 'Map';
    else jsViewers[i]['group'] = 'Misc';
    i++;
  }
}

function findViewer(value: string) {
  rootViewers.innerHTML = '';
  if (value != '') {
    recentBlock.style.display = 'none';
    for (const i in viewers) {
      if (viewers[i].name.toLowerCase().includes(value.toLowerCase()) ||
        viewers[i].group.toLowerCase().includes(value.toLowerCase()))
        rootViewers.append(renderCard(viewers, i));
    }
    for (const i in jsViewers) {
      if (jsViewers[i].name.toLowerCase().includes(value.toLowerCase()))
        rootViewers.append(renderCard(jsViewers, i));
      else if (jsViewers[i].group.toLowerCase().includes(value.toLowerCase()))
        rootViewers.append(renderCard(jsViewers, i));
      else if (jsViewers[i].description.toLowerCase().includes(value.toLowerCase()))
        rootViewers.append(renderCard(jsViewers, i));
      else if (jsViewers[i].keywords.toLowerCase().includes(value.toLowerCase()))
        rootViewers.append(renderCard(jsViewers, i));
      else if (jsViewers[i].package.toLowerCase().includes(value.toLowerCase()))
        rootViewers.append(renderCard(jsViewers, i));
    }
  } else {
    if (recentViewersRoot.childElementCount > 0)
      recentBlock.style.display = 'flex';
    for (const i in viewers)
      rootViewers.append(renderCard(viewers, i));

    for (const i in jsViewers)
      rootViewers.append(renderCard(jsViewers, i));
  }
  getTotalViewer();
  setTabIndex(rootViewers);
}

function setTabIndex(root: HTMLDivElement) {
  for (let i = 0; i < root.childElementCount; i++)
    root.children[i].setAttribute('tabindex', String(i + 1));
}

function renderCard(viewers: { [v: string]: { [k: string]: any } }, index: string) {
  let icon: HTMLElement;
  const viewer = viewers[index];

  if (viewer.type == 'viewer') {
    icon = ui.iconFA('');
    icon.className = 'grok-icon svg-icon ' + viewer.icon;
  }
  if (viewer.type == 'js-viewer') {
    if (viewer.icon != 'svg-project') {
      icon = ui.iconImage('', viewer.icon);
      icon.classList.add('svg-icon');
    } else {
      icon = ui.iconFA('');
      icon.className = 'grok-icon svg-icon ' + viewer.icon;
      icon.style.filter = 'grayscale(.5) invert(.5) contrast(120%)';
    }
  }
  const name = ui.label(viewer.name);
  name.classList.add('card-label');
  const card = ui.div([
    ui.divH([icon!, name]),
  ], `d4-item-card viewer-gallery vg-card-small${viewer.enabled ? '' : ' disabled'}`);

  if (viewer.enabled) {
    card.addEventListener('click', () => {
      dockViewers(viewer.name, view, table);
      dlg.close();
    });

    card.addEventListener('keydown', (e) => {
      if (e.key === 'Enter') {
        e.preventDefault();
        card.click();
      }
    });
  }
  else
    ui.tooltip.bind(card, viewer.tooltip);

  if (viewer.name.length > 18) {
    name.addEventListener('mouseenter', () => {
      if (name.offsetWidth < name.scrollWidth)
        name.setAttribute('title', viewer.name);
    }, {once: true});
  }

  return card;
}

function generateTags() {
  const root = ui.div([]);
  const list = [];
  for (const i in viewers)
    list.push(viewers[i].group);

  for (const i in jsViewers)
    list.push(jsViewers[i].group);

  for (const i in jsViewers)
    list.push(jsViewers[i].package);

  const tags = [...new Set(list)];

  for (const i in tags) {
    const tag = ui.div(tags[i], 'd4-tag');
    tag.addEventListener('click', () => {
      search.value = tags[i];
      search.fireChanged();
      search.input.focus();
    });
    root.append(tag);
  }

  return root;
}

function getTotalViewer() {
  viewersCount.innerHTML = '';
  viewersCount.append(`Viewers: ${rootViewers.childElementCount}`);
}

function dockViewers(viewer: string, view: DG.TableView, table: DG.DataFrame) {
  try {
    view.addViewer(DG.Viewer.fromType(viewer, table));
    logRecentViewer(viewer);
  } catch (e) {
    grok.shell.error(`Cannot add ${viewer} for current table. Check browser console for specific errors`);
    console.error(e);
  }
}

const STORAGE_NAME = 'recentViewerSettings';
const RECENT_KEY = 'recentViewers';

function getRecentViewersList(): string[] {
  try {
    const jsonStr = grok.userSettings.getValue(STORAGE_NAME, RECENT_KEY);
    if (!jsonStr)
      return [];
    return JSON.parse(jsonStr);
  } catch (e) {
    console.error('Failed to load recent viewers', e);
    return [];
  }
}

function logRecentViewer(viewerName: string) {
  let recent = getRecentViewersList();
  recent = [viewerName, ...recent.filter((n) => n !== viewerName)].slice(0, 8);
  grok.userSettings.add(STORAGE_NAME, RECENT_KEY, JSON.stringify(recent));
}
