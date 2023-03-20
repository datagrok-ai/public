import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import { DataFrame, InputBase } from 'datagrok-api/dg';

let view: DG.TableView;
let table: DataFrame;

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

const viewers: any = {};
const jsViewers: any = {};
const rootViewers = ui.divH([], 'viewer-gallery');
let dlg: DG.Dialog;
let search: InputBase;
const viewersCount = ui.div([], 'vg-counter-label');

export function viewersDialog(currentView: DG.TableView, currentTable: DataFrame) {
    getViewers(viewers);
    getJsViewers(jsViewers);

    view = currentView;
    table = currentTable;

    dlg = ui.dialog('Add Viewer');
    $(dlg.root.lastChild).hide();

    search = ui.searchInput('', '', (v: string) => findViewer(v));
    search.input.setAttribute('tabindex', '-1');
    search.input.setAttribute('placeholder', 'Search by name, keywords, description, tag, or package');

    search.input.addEventListener('keydown', (e) => {
        if (e.key === 'Escape') {
            e.preventDefault();
            search.value = '';
            search.fireChanged();
        }
    });

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
    const viewerList = [];

    for (const [key, value] of Object.entries(DG.VIEWER)) {
        switch (String(value)) {
            case 'Globe': break;
            case 'Google map': break;
            case 'RadarViewer': break;
            case 'SurfacePlot': break;
            case 'TimelinesViewer': break;
            case 'Word cloud': break;
            default:
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

    for (const i in viewerList) {
        Object.assign(viewers, {
            [i]: {
                name: viewerList[i],
                icon: 'grok-icon svg-icon svg-' + viewerList[i].toLowerCase().replace(/(\s)/g, '-'),
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

function getJsViewers(jsViewers: { [v: string]: { [k: string]: any } }) {
    const list = DG.Func.find({ returnType: 'viewer' });
    let i = 0;
    for (const v of list) {
        if (v.package.name != 'ApiTests') {
            Object.assign(jsViewers, {
                [i]: {
                    name: v.friendlyName,
                    icon: (v.options['icon'] != undefined) ? `${v.package.webRoot}${v.options['icon']}` : 'svg-project',
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
        }
        i++;
    }
}

function findViewer(value: string) {
    rootViewers.innerHTML = '';
    if (value != '') {
        for (const i in viewers) {
            if (viewers[i].name.toLowerCase().includes(value.toLowerCase()) || viewers[i].group.toLowerCase().includes(value.toLowerCase()))
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
    if (viewers[index].type == 'viewer') {
        icon = ui.iconFA('');
        icon.className = 'grok-icon svg-icon ' + viewers[index].icon;
    }
    if (viewers[index].type == 'js-viewer') {
        if (viewers[index].icon != 'svg-project') {
            icon = ui.iconImage('', viewers[index].icon);
            icon.classList.add('svg-icon');
        } else {
            icon = ui.iconFA('');
            icon.className = 'grok-icon svg-icon ' + viewers[index].icon;
            icon.style.filter = 'grayscale(.5) invert(.5) contrast(120%)';
        }
    }
    const name = ui.label(viewers[index].name);
    name.classList.add('card-label');
    const card = ui.div([
        ui.divH([icon!, name]),
    ], 'd4-item-card viewer-gallery vg-card-small');

    card.addEventListener('click', () => {
        dockViewers(viewers[index].name, view, table);
        dlg.close();
    });

    card.addEventListener('keydown', (e) => {
        if (e.key === 'Enter') {
            e.preventDefault();
            card.click();
        }
    });

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
    view.addViewer(DG.Viewer.fromType(viewer, table));
}