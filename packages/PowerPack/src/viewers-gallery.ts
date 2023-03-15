import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import { DataFrame, InputBase, TableView, View } from 'datagrok-api/dg';

let view: DG.TableView;
let table: DataFrame;

const group_comparisons = [
    'Bar chart',
    'Radar',
    'Chord viewer',
    'Matrix plot',
    'Word cloud',
    'Word cloud viewer',
    'Pie chart',
    'Tree map',
];
const group_trends = [
    'Line chart',
    'Histogram',
    'Box plot',
    'Sunburst',
];
const group_correlations = [
    'Scatter plot',
    'Density plot',
    '3d scatter plot',
    'PC plot',
    'Correlation plot',
    'Heat map',
];
const group_relationships = [
    'Sankey',
    'Network diagram',
    'Tree',
    'Phylo tree',
    'Scaffold Tree',
    'PhylocanvasGL'
];
const group_maps = [
    'Shape map',
    'Leaflet',
    'Globe',
    'Map',
];

const viewer_list = [];

let viewers: any = {};
let jsViewers: any = {};
let rootViewers = ui.divH([], 'viewer-gallery');
let dlg: DG.Dialog;
let search: InputBase;
let viewersCount = ui.div([], 'vg-counter-label');

export async function viewersDialog(current_view: DG.TableView, current_table: DataFrame) {
    getViewers(viewers);
    getJsViewers(jsViewers);

    view = current_view;
    table = current_table;

    dlg = ui.dialog('Add Viewer');
    $(dlg.root.lastChild).hide();

    search = ui.searchInput('', '', (v: string) => findViewer(v));
    search.input.setAttribute('tabindex', '-1');
    search.input.setAttribute('placeholder', 'Search by name, keywords, description, tag, or package');

    search.input.addEventListener('keydown', (e) => {
        if (e.key === "Escape") {
            e.preventDefault();
            search.value = '';
            search.fireChanged();
        }
    });

    let searchIcon = ui.iconFA('search');
    searchIcon.classList.add('vg-search-icon');

    rootViewers.innerHTML = '';

    for (let i in viewers) {
        rootViewers.append(renderCard(viewers, i));
    }
    for (let i in jsViewers) {
        rootViewers.append(renderCard(jsViewers, i));
    }

    let tags = ui.divV([
        ui.label('Relative tags:'),
        generateTags()
    ], 'vg-tags');

    let root = ui.divH([
        ui.block([viewersCount, rootViewers], 'viewer-gallery-root'),
        tags
    ]);

    dlg.add(ui.block([ui.div([searchIcon, search.input], 'd4-search-ba')], 'vg-controls grok-gallery-search-bar'));
    dlg.add(root);
    dlg.showModal(true);

    getTotalViewer();
    setTabIndex(rootViewers);

};

function getViewers(viewers: { [v: string]: { [k: string]: any } }) {
    let i = 0;
    let viewer_list = [];

    for (const [key, value] of Object.entries(DG.VIEWER)) {
        switch (String(value)) {
            case 'Globe': break;
            case 'Google map': break;
            case 'RadarViewer': break;
            case 'SurfacePlot': break;
            case 'TimelinesViewer': break;
            case 'Word cloud': break;
            default:
                viewer_list.push(value);
                break;
        }

    }
    viewer_list.push('Pivot table');
    viewer_list.push('Column viewer');
    viewer_list.push('Web viewer');
    viewer_list.push('Card');
    viewer_list.push('Viewer host');
    viewer_list.push('Info panel');
    viewer_list.push('Scripting viewer');

    for (let i in viewer_list) {
        Object.assign(viewers, {
            [i]: {
                name: viewer_list[i],
                icon: 'grok-icon svg-icon svg-' + viewer_list[i].toLowerCase().replace(/(\s)/g, '-'),
                group: '',
                type: 'viewer'
            }
        });
        if (group_comparisons.includes(viewer_list[i])) { viewers[i]['group'] = 'Comparison' }
        else if (group_correlations.includes(viewer_list[i])) { viewers[i]['group'] = 'Correlation' }
        else if (group_relationships.includes(viewer_list[i])) { viewers[i]['group'] = 'Relationship' }
        else if (group_trends.includes(viewer_list[i])) { viewers[i]['group'] = 'Trend' }
        else if (group_maps.includes(viewer_list[i])) { viewers[i]['group'] = 'Map' }
        else { viewers[i]['group'] = 'Misc' }
    }
}

function getJsViewers(jsViewers: { [v: string]: { [k: string]: any } }) {
    let list = DG.Func.find({ returnType: 'viewer' });
    let i = 0;
    for (let v of list) {
        if (v.package.name != 'ApiTest') {
            Object.assign(jsViewers, {
                [i]: {
                    name: v.friendlyName,
                    icon: (v.options['icon'] != undefined) ? `${v.package.webRoot}${v.options['icon']}` : 'svg-project',
                    description: v.description,
                    keywords: (v.options['keywords'] != null) ? v.options['keywords'] : '',
                    group: '',
                    package: v.package.name,
                    type: 'js-viewer'
                }
            });
            if (group_comparisons.includes(v.friendlyName)) { jsViewers[i]['group'] = 'Comparison' }
            else if (group_correlations.includes(v.friendlyName)) { jsViewers[i]['group'] = 'Correlation' }
            else if (group_relationships.includes(v.friendlyName)) { jsViewers[i]['group'] = 'Relationship' }
            else if (group_trends.includes(v.friendlyName)) { jsViewers[i]['group'] = 'Trend' }
            else if (group_maps.includes(v.friendlyName)) { jsViewers[i]['group'] = 'Map' }
            else { jsViewers[i]['group'] = 'Misc' }
        }
        i++;
    }
}

function findViewer(value: string) {
    rootViewers.innerHTML = '';
    if (value != '') {
        for (const i in viewers) {
            if (viewers[i].name.toLowerCase().includes(value.toLowerCase()) || viewers[i].group.toLowerCase().includes(value.toLowerCase())) {
                rootViewers.append(renderCard(viewers, i));
            }
        }
        for (const i in jsViewers) {
            if (jsViewers[i].name.toLowerCase().includes(value.toLowerCase())) {
                rootViewers.append(renderCard(jsViewers, i));
            } else if (jsViewers[i].group.toLowerCase().includes(value.toLowerCase())) {
                rootViewers.append(renderCard(jsViewers, i));
            } else if (jsViewers[i].description.toLowerCase().includes(value.toLowerCase())) {
                rootViewers.append(renderCard(jsViewers, i));
            } else if (jsViewers[i].keywords.toLowerCase().includes(value.toLowerCase())) {
                rootViewers.append(renderCard(jsViewers, i));
            } else if (jsViewers[i].package.toLowerCase().includes(value.toLowerCase())) {
                rootViewers.append(renderCard(jsViewers, i));
            }
        }
    } else {
        for (const i in viewers) {
            rootViewers.append(renderCard(viewers, i))
        }
        for (const i in jsViewers) {
            rootViewers.append(renderCard(jsViewers, i))
        }
    }
    getTotalViewer();
    setTabIndex(rootViewers);
}

function setTabIndex(root: HTMLDivElement) {
    for (let i = 0; i < root.childElementCount; i++) {
        root.children[i].setAttribute('tabindex', String(i + 1));
    }
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
    let name = ui.label(viewers[index].name);
    name.classList.add('card-label');
    let card = ui.div([
        ui.divH([icon!, name])
    ], 'd4-item-card viewer-gallery vg-card-small');

    card.addEventListener('click', () => {
        dockViewers(viewers[index].name, view, table);
        dlg.close();
    });

    card.addEventListener('keydown', (e) => {
        if (e.key === "Enter") {
            e.preventDefault();
            card.click();
        }
    });

    return card;
}

function generateTags() {
    let root = ui.div([]);
    let list = [];
    for (const i in viewers) {
        list.push(viewers[i].group)
    }
    for (const i in jsViewers) {
        list.push(jsViewers[i].group)
    }
    for (const i in jsViewers) {
        list.push(jsViewers[i].package)
    }
    const tags = [...new Set(list)]

    for (let i in tags) {
        let tag = ui.div(tags[i], 'd4-tag');
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
    let count = 0;
    for (let i = 0; i < rootViewers.childElementCount; i++) {
        count++;
    }
    viewersCount.innerHTML = '';
    viewersCount.append(`Viewers: ${count}`);
}

function dockViewers(viewer: string, view: DG.TableView, table: DG.DataFrame) {
    view.addViewer(DG.Viewer.fromType(viewer, table));
}
