import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import '../css/viewers-gallery.css';
import { filter } from 'rxjs';
import { DataFrame, TableView } from 'datagrok-api/dg';


export async function viewersGallery() {

    let viewers = {};
    const table = grok.shell.t;
    const view = grok.shell.tableView(table.name);
    let root = ui.div([], {style:{position:'relative', paddingTop:'40px'}});
    let recommends = ui.divH([],'viewer-gallery');
    let selectioText = ui.div(['Selected: 0'],'vg-selection-text');
    let cards = ui.divH([],'viewer-gallery');
    let search = ui.searchInput('','',(value)=>{
        clearRoot([recommends]);
        clearRoot([cards]);
        unSelectAll(viewers)
        if (value != ''){
            for (let i in viewers){
                if (viewers[i].name.toLowerCase().includes(value.toLowerCase())){
                    if (viewers[i].recommend)
                        recommends.append(render(viewers[i],table))
                    else
                        cards.append(render(viewers[i],table));
                }
            }
        } else {
            for (let i in viewers){
                if (viewers[i].recommend)
                    recommends.append(render(viewers[i],table))
                else
                    cards.append(render(viewers[i],table));
            }
        }
        $('.vg-selection-text').html('Selected: '+ $('.viewer-gallery').find('.fa-minus').length);   
    });

    //@ts-ignore
    search.input.placeholder = 'Search by name or type';
    
    let controls = ui.block([ui.div([search.input],'d4-search-ba')], 'vg-controls grok-gallery-search-bar');

    for (let i=0; i<DG.Viewer.getViewerTypes().length; i++){

        Object.assign(viewers, {[i]:{
            icon:'grok-icon svg-icon svg-'+DG.Viewer.getViewerTypes()[i].toLowerCase().replace(/(\s)/g,'-'),
            name: insertSpaces(DG.Viewer.getViewerTypes()[i]),
            type: 'Viewer',
            recommend: false,
            dock: false
        }});

        if (i<5){
            viewers[i].recommend = true;
        }
        if(viewers[i].name.includes('widget')){
            viewers[i].type = 'Widget';
            viewers[i].icon = 'grok-icon svg-icon svg-project';
        }

    }

    let dlg = ui.dialog('Add viewer')
        .add(ui.divV([
            controls,
            ui.h1('Recommends'),
            recommends,
            ui.h1('All viewers'),
            cards,
        ], 'viewer-gallery-root'))
        .onOK(() => dockViewers(viewers, table, view));

    $(dlg.root).find('.d4-dialog-contents').removeClass('ui-form');
    $(dlg.root).find('.d4-dialog-contents').removeClass('ui-panel');
    $(dlg.root).find('.d4-dialog-contents').addClass('ui-box');
    console.log($(dlg.root).find('.d4-command-bar'));
    $(dlg.root).find('.d4-command-bar').append(selectioText);


    if (grok.shell.v.type == 'TableView'){
        dlg.show({y: 50, width:940, height:620})
    }

    for (let i in viewers){
        if (viewers[i].recommend)
            recommends.append(render(viewers[i],table))
        else
            cards.append(render(viewers[i],table));    
    }

    (async () => {
        console.log($(root).find('i'));
    })();

}

function insertSpaces(string) {
    string = string.replace(/^(_|-)/g, '')
	string = string.replace(/(_|-)/g, ' ');
    string = string.replace(/([a-z]\s)/g, m => m.toLowerCase());
    string = string.replace(/([a-z])([A-Z])/g, '$1 $2');
    string = string.replace(/([A-Z])([A-Z][a-z])/g, '$1 $2');
    string = string.replace(/([A-Z])([a-z])/g, m => m.toLowerCase());
    string = string.replace(/^./g, m => m.toUpperCase())
    return string
}

function clearRoot(root: HTMLDivElement[]){
    for (const i in root){
        root[i].innerHTML = '';
    }
}

function render(viewer:any, table:DG.DataFrame){
    let root = ui.div([])
    let icon = ui.iconFA('');
    icon.className = 'grok-icon svg-icon '+viewer.icon;
    let label = ui.div([viewer.name], 'card-label');

    let addBtn = ui.button(ui.iconFA('plus'), ()=>{
        if($(addBtn).find('i').hasClass('fa-plus')){
            viewer.dock = true
            root.style.border = '1px solid var(--blue-1)';
            $(addBtn).find('i').removeClass('fa-plus');
            $(addBtn).find('i').addClass('fa-minus')
        } else {
            viewer.dock = false
            root.style.border = '1px solid var(--grey-2)';
            $(addBtn).find('i').removeClass('fa-minus');
            $(addBtn).find('i').addClass('fa-plus')
        }
        $('.vg-selection-text').html('Selected: '+ $('.viewer-gallery').find('.fa-minus').length);
    });
    
    if (viewer.recommend){
        let viewerRoot = ui.box(DG.Viewer.fromType(viewer.name, table, {"LegendVisibility":"Never"}).root);
        root.className = 'd4-item-card viewer-gallery vg-card';
        root.append(ui.block([viewerRoot,ui.divH([label,addBtn])]));
    }else{
        root.className = 'd4-item-card viewer-gallery vg-card-small';
        root.append(ui.divH([icon,label,addBtn]));
    }

    root.addEventListener('click', ()=>addBtn.click());
    return root;
}

function unSelectAll(viewers:object){
    for (let i in viewers){
        viewers[i].dock = false
    }
}

function dockViewers(viewers:object, table:DG.DataFrame, view: DG.TableView){
    for (let i in viewers){
        if(viewers[i].dock)
            view.addViewer(DG.Viewer.fromType(viewers[i].name, table));
    }
}