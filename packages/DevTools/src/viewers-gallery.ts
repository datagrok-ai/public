import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import '../css/viewers-gallery.css';
import { filter } from 'rxjs';
import { DataFrame, TableView } from 'datagrok-api/dg';let viewerOptions = {
    "showXAxis": false,
    "showYAxis": false,
    "showXSelector": false,
    "showYSelector": false,
    "showYSelectors": false,
    "showXSelectors": false,
    "showSizeSelector": false,
    "showColorSelector": false,
    "showVerticalGridLines": false,
    "showHorizontallGridLines": false,
    "showContextMenu": false,
    "showColumnSelector": false,
    "showBinSelector": false,
    "showCurrentRow": false,
    "showRangeSlider": false,
    "showRangeInputs": false,
    "showTopPanel": false,
    "showAggrSelectors": false,
    "showValueAxis": false,
    "showValueSelector": false,
    "showCategorySelector": false,
    "showStackSelector": false,
    "multiAxis": true,
    "showSplitSelector": false,
    "showColumnSelectionPanel": false,
    "legendVisibility":"Never",
    "marginLeft": 10,
    "marginRight": 10,
    "marginTop": 10,
    "marginBottom": 10,
  }
  let tempName = ''; 

  let contentBox = ui.box();
  contentBox.classList.add('vg-root-content');

  let descriptionBox = ui.box();
  descriptionBox.classList.add('vg-root-description');

  let filterBox = ui.box();
  filterBox.classList.add('vg-root-filters');
  
  filterBox.style.maxWidth = '230px';

export async function viewersGallery() {

    let viewers = {};
    const table = grok.shell.t;
    const view = grok.shell.tableView(table.name);
    let root = ui.box();
    
    let recommends = ui.divH([],'viewer-gallery');
    let selectioText = ui.div(['Selected: '],'vg-selection-text');
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
    
    let searchBlock = ui.block([ui.div([search.input],'d4-search-ba')], 'vg-controls grok-gallery-search-bar');
    
    clearRoot([filterBox,contentBox,descriptionBox]);
    
    for (let i=0; i<DG.Viewer.getViewerTypes().length; i++){

        Object.assign(viewers, {[i]:{
            icon:'grok-icon svg-icon svg-'+DG.Viewer.getViewerTypes()[i].toLowerCase().replace(/(\s)/g,'-'),
            name: insertSpaces(DG.Viewer.getViewerTypes()[i]),
            type: 'Viewer',
            recommend: false,
            dock: false
        }});

        if (i<7){
            viewers[i].recommend = true;
        }
        if(viewers[i].name.includes('widget')){
            viewers[i].type = 'Widget';
            viewers[i].icon = 'grok-icon svg-icon svg-project';
        }

    }

    filterBox.append(ui.divV([
        ui.h3('Type'),
        ui.boolInput('Viewer',false),
        ui.boolInput('Widget',false),
        ui.h3('Category'),
        ui.boolInput('Comparison',false),
        ui.boolInput('Trends',false),
        ui.boolInput('Part to whole',false),
        ui.boolInput('Correlations',false),
        ui.boolInput('Relationships',false),
        ui.boolInput('Maps',false),
    ], 'vg-filter-panel'));
    

    contentBox.append(
        ui.divV([
            searchBlock,
            ui.h1('Recommends'),
            recommends,
            ui.h1('All viewers'),
            cards,
        ], 'viewer-gallery-root')
    );

    $(descriptionBox).hide();

    //root.innerHTML = '';
    root = ui.splitH([
        filterBox,
        contentBox,
        descriptionBox
    ])
    ui.dialog('Add viewer')
    let dlg = ui.dialog('Add viewer')
        .add(root
            //viewerRoot
        )
        .onOK(() => {
            //clearRoot([root]) 
            if(tempName!='')
                view.addViewer(DG.Viewer.fromType(tempName, table))
        });

    $(dlg.root).find('.d4-dialog-contents').removeClass('ui-form');
    $(dlg.root).find('.d4-dialog-contents').removeClass('ui-panel');
    $(dlg.root).find('.d4-dialog-contents').addClass('ui-box');
    $(dlg.root).find('.d4-command-bar').append(selectioText);
    
    $(dlg.root).find('.d4-dialog-contents').css('padding','0px');
    $(dlg.root).find('.d4-dialog-header').css('border-bottom','1px solid var(--grey-2)');
    $(dlg.root).find('.d4-command-bar').css('border-top','1px solid var(--grey-2)');

    if (grok.shell.v.type == 'TableView'){
        dlg.showModal(true)
    }

    for (let i in viewers){
        if (viewers[i].recommend)
            recommends.append(render(viewers[i],table))
        else
            cards.append(render(viewers[i],table));    
    }

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

    let details = ui.button(ui.iconFA('list'), ()=>{
        //show details block
        $(descriptionBox).show();
        $(contentBox).hide();
        //$(filterBox).hide();
        clearRoot([descriptionBox]);
        //$(descriptionBox).html('');
        let markup  = ui.div();
        //@ts-ignore
        let link = DG.Viewer.fromType(viewer.name, table).helpUrl;
        grok.dapi.fetchProxy('https://raw.githubusercontent.com/datagrok-ai/public/master/'+link)
            .then(response => response.text())
            .then(data => markup.append(ui.markdown(data)));

        $(descriptionBox).append(ui.splitH([
            ui.box(DG.Viewer.fromType(viewer.name, table).root),
            ui.divV([
                ui.button('back', ()=>{
                    $(descriptionBox).hide();
                    $(contentBox).show();
                    //$(filterBox).show();
                }),
                
                markup
                //DG.Viewer.fromType('Web widget', table).root,
            ])
        ]))
    });

    if (viewer.recommend){
        let viewerRoot = ui.box(DG.Viewer.fromType(viewer.name, table, viewerOptions).root);
        root.className = 'd4-item-card viewer-gallery vg-card';
        root.append(ui.block([viewerRoot,ui.divH([label,details])]));
    }else{
        root.className = 'd4-item-card viewer-gallery vg-card-small';
        root.append(ui.divH([icon,label,details]));
    }
    root.addEventListener('click', ()=>{
        $('.vg-card').css('border','1px solid var(--grey-2)');
        $('.vg-card-small').css('border','1px solid var(--grey-2)');
        root.style.border = '1px solid var(--blue-1)';
        tempName = viewer.name;
        $('.vg-selection-text').html('Selected: '+ viewer.name);
    });

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