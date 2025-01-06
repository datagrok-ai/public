import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import $ from 'cash-dom';

let view: DG.TableView;
let table: DG.DataFrame;

const groupСomparisons = [
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

const viewers: any = {};
const jsViewers: any = {};
const rootViewers = ui.divH([], 'viewer-gallery');
let dlg: DG.Dialog;
let search: DG.InputBase;
const viewersCount = ui.div([], 'vg-counter-label');

export function viewersDialog(currentView: DG.TableView, currentTable: DG.DataFrame) {
  getViewers(viewers);
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
    ui.block([viewersCount, rootViewers], 'viewer-gallery-root'),
    tags,
  ]);

  dlg.add(ui.block([ui.div([searchIcon, search.input], 'd4-search-ba')], 'vg-controls grok-gallery-search-bar'));
  dlg.add(root);
  dlg.showModal(true);

  getTotalViewer();
  setTabIndex(rootViewers);
};

function getViewers(viewers: { [v: string]: { [k: string]: any } }) {
  let viewerList = [];

  for (const value of Object.values(DG.VIEWER)) {
    switch (String(value)) {
    case 'Globe': break;
    case 'Google map': break;
    case 'RadarViewer': break;
    case 'SurfacePlot': break;
    case 'TimelinesViewer': break;
    case 'Word cloud': break;
    case 'Scaffold Tree': break;
    default:
      if (value !== 'Shape Map') // return Shape Map back when it is reincarnated
        viewerList.push(value);
      break;
    }
  }
  viewerList.push('Pivot table');
  viewerList.push('Column viewer');
  viewerList.push('Web viewer');
  viewerList.push('Card');
  viewerList.push('Viewer host');
  viewerList.push('Info panel');
  viewerList.push('Scripting viewer');

  viewerList = [...new Set(viewerList)];
  for (const i in viewerList) {
    Object.assign(viewers, {
      [i]: {
        name: viewerList[i],
        icon: 'grok-icon svg-icon svg-' + viewerList[i].toLowerCase().replace(/(\s)/g, '-'),
        enabled: true,
        group: '',
        type: 'viewer',
      },
    });
    if (groupСomparisons.includes(viewerList[i])) viewers[i]['group'] = 'Comparison';
    else if (groupCorrelations.includes(viewerList[i])) viewers[i]['group'] = 'Correlation';
    else if (groupRelationships.includes(viewerList[i])) viewers[i]['group'] = 'Relationship';
    else if (groupTrends.includes(viewerList[i])) viewers[i]['group'] = 'Trend';
    else if (groupMaps.includes(viewerList[i])) viewers[i]['group'] = 'Map';
    else viewers[i]['group'] = 'Misc';
  }
}

function getJsViewers(jsViewers: { [v: string]: { [k: string]: any } }, table: DG.DataFrame) {
  const skip = ['TestViewerForProperties', 'OutliersSelectionViewer'];
  const list = DG.Func.find({tags: ['viewer']}).filter((v) => !skip.includes(v.friendlyName));
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
    }
    Object.assign(jsViewers, {
      [i]: {
        name: v.friendlyName,
        enabled: isViewerEnabled,
        tooltip: isViewerEnabled ? '' : CHARTS_VIEWERS_TEST_DATA[v.friendlyName].tooltip,
        icon: (v.options['icon'] != undefined) ? `${v.package.webRoot.endsWith('/') ?
          v.package.webRoot : v.package.webRoot + '/'}${v.options['icon']}` : 'svg-project',
        description: v.description,
        keywords: (v.options['keywords'] != null) ? v.options['keywords'] : '',
        group: '',
        package: v.package.name,
        type: 'js-viewer',
      },
    });
    if (groupСomparisons.includes(v.friendlyName)) jsViewers[i]['group'] = 'Comparison';
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
  } catch (e) {
    grok.shell.error(`Cannot add ${viewer} for current table. Check browser console for specific errors`);
    console.error(e);
  }
}
