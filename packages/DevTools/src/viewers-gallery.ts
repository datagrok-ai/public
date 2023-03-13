import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import '../css/viwes-gallery.css';
import { DataFrame, InputBase, TableView, View } from 'datagrok-api/dg';

const table = grok.shell.t;
const view = grok.shell.tableView(table.name);

const group_comparisons = [
    'Bar chart',
    'Radar viewer',
    'Chord viewer',
    'Matrix plot',
    'Word cloud',
    'Word cloud viewer',
    'Pie chart',
    'Tree map',
    'Tree map viewer'
];
const group_trends = [
    'Line chart',
    'Histogram',
    'Box plot',
    'Sunbrust viewer',
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
    'Sankey viewer',
    'Network diagram',
    'Tree viewer',
    'Phylo tree'
];
const group_maps = [
    'Shape map',
    'Google map',
    'Globe',
    'Map',
];

let viewers: any = {};
let jsViewers: any = {};
let root = ui.div([], 'viewer-gallery-root');
let dlg: DG.Dialog;
let search: InputBase;

export async function viewersDialog() {
    getViewers(viewers);
    getJsViewers(jsViewers);

    dlg = ui.dialog('Viewer Gallery');
    search = ui.searchInput('','',(v)=>findViewer(v));
    search.input.setAttribute('tabindex','-1');

    for (let i in viewers){
        root.append(renderCard(viewers, i));
    }
    for (let i in jsViewers){
        root.append(renderCard(jsViewers, i));
    }
    dlg.add(search);
    dlg.add(generateTags())
    dlg.add(root);
    dlg.showModal(true);
    setTabIndex(root);
};

function getViewers(viewers:{}){
    let i = 0;
    for(const [key, value] of Object.entries(DG.VIEWER)){
        if (key == 'SURFACE_PLOT') null
        else if (key == 'RADAR_VIEWER') null
        else if (key == 'TIMELINES') null
        else{
            Object.assign(viewers,{
                [i]:{
                    name: value,
                    icon: 'grok-icon svg-icon svg-' + value.toLowerCase().replace(/(\s)/g, '-'),
                    group: '',
                    type: 'viewer'
                }
            });
            if (group_comparisons.includes(value)){viewers[i].group='Comparisons'}
            else if (group_correlations.includes(value)){viewers[i].group='Correlations'}
            else if (group_relationships.includes(value)){viewers[i].group='Relationships'}
            else if (group_trends.includes(value)){viewers[i].group='Trends'}
            else if (group_maps.includes(value)){viewers[i].group='Map'}
            else {viewers[i].group='Other'}
        }
        i++;
    }
}

function getJsViewers(jsViewers:{}){
    let list = DG.Func.find({returnType:'viewer'});
    let i = 0;
    for (let v of list){
        Object.assign(jsViewers,{
            [i]:{
                name: v.friendlyName,
                icon: (v.options['icon']!= undefined) ? `${v.package.webRoot}${v.options['icon']}`: 'svg-project',
                description: v.description,
                keywords: (v.options['keywords']!=null)?v.options['keywords']:'',
                group: v.package.name,
                type: 'js-viewer'
            }
        });
        i++;
    }
}

function findViewer (value:string) {
    root.innerHTML = '';
    if (value != null) {
        for (const i in viewers){
            if (viewers[i].name.toLowerCase().includes(value.toLowerCase()) || viewers[i].group.toLowerCase().includes(value.toLowerCase())) {
                root.append(renderCard(viewers, i));
            }
        }
        for (const i in jsViewers){
            if (jsViewers[i].name.toLowerCase().includes(value.toLowerCase())) {
                root.append(renderCard(jsViewers, i));
            } else if (jsViewers[i].group.toLowerCase().includes(value.toLowerCase())) {
                root.append(renderCard(jsViewers, i));
            } else if (jsViewers[i].description.toLowerCase().includes(value.toLowerCase())) {
                root.append(renderCard(jsViewers, i));
            } else if (jsViewers[i].keywords.toLowerCase().includes(value.toLowerCase()))    {
                root.append(renderCard(jsViewers, i));
            }
        }
    } else {
        for (const i in viewers){
            root.append(ui.div([viewers[i].name], {style:{padding:'6px', border:'1px solid var(--grey-2)'}}))
        }
        for (const i in jsViewers){
            root.append(ui.div([viewers[i].name], {style:{padding:'6px', border:'1px solid var(--grey-2)'}}))
        }
    }
    setTabIndex(root);
}

function setTabIndex(root: HTMLDivElement) {
    for (let i = 0; i < root.childElementCount; i++){
        root.children[i].setAttribute('tabindex', String(i+1));
    }
}

function renderCard(viewers:{}, index:string) {
    let icon: HTMLElement;
    if (viewers[index].type == 'viewer'){
        icon = ui.iconFA('');
        icon.className = 'grok-icon svg-icon ' + viewers[index].icon;
    }
    if (viewers[index].type == 'js-viewer'){
        if (viewers[index].icon != 'svg-project'){
            icon = ui.iconImage('',viewers[index].icon);
        } else {
            icon = ui.iconFA('');
            icon.className = 'grok-icon svg-icon ' + viewers[index].icon;
            icon.style.filter = 'grayscale(.5) invert(.5) contrast(120%)';
        }
    }
    let name = ui.label(viewers[index].name);
    let card = ui.divH([
        icon, name  
    ]);

    card.style.margin = '5px';
    card.style.padding = '5px';
    card.style.gap = '5px';

    card.addEventListener('click', ()=>{
        dockViewers(viewers[index].name, table, view);
        dlg.close();
    });

    card.addEventListener('keydown', (e)=>{
        if (e.key === "Enter") {
            e.preventDefault();
            card.click();
        }
    });

    return card;
}

function generateTags(){
    let root = ui.div([]);
    let list = [];
    for (const i in viewers){
        list.push(viewers[i].group)
    }
    for (const i in jsViewers){
        list.push(jsViewers[i].group)
    }
    const tags = [...new Set(list)]
    
    for (let i in tags){
        let tag = ui.div(tags[i], 'd4-tag');
        tag.addEventListener('click', ()=>{
            search.value = tags[i];
            search.fireChanged();
            search.input.focus();
        });
        root.append(tag);
    }

    return root;
}

function dockViewers(viewer: string, table: DG.DataFrame, view: DG.TableView) {
    view.addViewer(DG.Viewer.fromType(viewer, table));
}
