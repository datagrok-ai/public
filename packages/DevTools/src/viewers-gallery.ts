import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import '../css/viewers-gallery.css';
import { filter } from 'rxjs';
import { DataFrame, InputBase, TableView } from 'datagrok-api/dg';
const cat_comparisons = [
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
const cat_trends = [
    'Line chart',
    'Histogram',
    'Box plot',
    'Sunbrust viewer',
];
const cat_correlations = [
    'Scatter plot',
    'Density plot',
    '3d scatter plot',
    'PC plot',
    'Correlation plot',
    'Heat map',
];
const cat_relationships = [
    'Sankey viewer',
    'Network diagram',
    'Tree viewer',
    'Phylo tree'
];
const cat_maps = [
    'Shape map',
    'Google map',
    'Globe',
];
const viewerOptions = {
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

  let dlgTitle;

  let okBtn;
  //let okBtn = ui.bigButton('ADD',()=>{});

  let backBtn = ui.iconFA('arrow-left',()=>{
    $('.vg-selection-text').html('Selected: ');
    dlgTitle.text('Add viewer');
    $(okBtn).hide();
    $(backBtn).hide();
    $(descriptionBox).hide();
    $(contentBox).show();
    $(filterBox).show();
  })
  backBtn.style.marginRight = '10px';
  backBtn.style.fontSize = '16px';
  backBtn.classList.remove('fal');
  backBtn.classList.add('far');

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
        $(okBtn).hide();
        tempName = '';
        let filter_value = getFilter([filters_type, filters_cat], viewers);

        if (filter_value.length === 0){
            for (let i in viewers){
                filter_value.push(viewers[i].name)
            }
        }

        if (value != ''){
            for (let i in viewers){
                if (filter_value.indexOf(viewers[i].name) != -1 && viewers[i].name.toLowerCase().includes(value.toLowerCase())){
                    if (viewers[i].recommend)
                        recommends.append(render(viewers[i],table))
                    else
                        cards.append(render(viewers[i],table));
                }
            }
        } 
        else {
            for (let i in viewers){
                if (filter_value.indexOf(viewers[i].name) != -1){
                    if (viewers[i].recommend)
                        recommends.append(render(viewers[i],table))
                    else
                        cards.append(render(viewers[i],table));
                }
            }
        }
        $('.vg-selection-text').html('Selected: '+ $('.viewer-gallery').find('.fa-minus').length);   
    });

    //@ts-ignore
    let filters_type = ui.multiChoiceInput('', [''], ['Viewer','Widget'], (value)=>{
        getFilter([filters_type, filters_cat], viewers);
        search.fireChanged();
    });

    //@ts-ignore
    let filters_cat = ui.multiChoiceInput('', [''], ['Comparisons','Trends','Correlations','Relationships', 'Maps', 'Others'], (value)=>{
        getFilter([filters_type, filters_cat], viewers);
        search.fireChanged();
        /*
        search.value = '';
        clearRoot([recommends]);
        clearRoot([cards]);
        $(okBtn).hide();
        tempName = '';

        if (value.length>0){
            for (let i in viewers){
                if (viewers[i].recommend){
                    if (value.indexOf(viewers[i].category) != -1){
                        recommends.append(render(viewers[i],table));
                    }
                }
                else{
                    if (value.indexOf(viewers[i].category) != -1){
                        cards.append(render(viewers[i],table));
                    }
                }
            }
        }else{
            for (let i in viewers){
                if (viewers[i].recommend)
                    recommends.append(render(viewers[i],table))
                else
                    cards.append(render(viewers[i],table));
            }
        }
        

        $('.vg-selection-text').html('Selected: '+ $('.viewer-gallery').find('.fa-minus').length); */
    })

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
            category: 'Other'
        }});

        if (i<7){
            viewers[i].recommend = true;
        }
        
        if (cat_comparisons.indexOf(viewers[i].name) != -1){
            viewers[i].category = 'Comparisons';
        }else if (cat_trends.indexOf(viewers[i].name) != -1){
            viewers[i].category = 'Trends';
        }else if (cat_correlations.indexOf(viewers[i].name) != -1){
            viewers[i].category = 'Correlations';
        }else if (cat_relationships.indexOf(viewers[i].name) != -1){
            viewers[i].category = 'Relationships';
        }else if (cat_maps.indexOf(viewers[i].name) != -1){
            viewers[i].category = 'Maps';
        }else{
            viewers[i].category = 'Others';
        }

        if(viewers[i].name.includes('widget')){
            viewers[i].type = 'Widget';
            viewers[i].icon = 'grok-icon svg-icon svg-project';
        }

    }

    filterBox.append(ui.divV([
        ui.h3('Type'),
        filters_type.root,
        ui.h3('Category'),
        filters_cat.root
       /* ui.boolInput('Comparison',false),
        ui.boolInput('Trends',false),
        ui.boolInput('Part to whole',false),
        ui.boolInput('Correlations',false),
        ui.boolInput('Relationships',false),
        ui.boolInput('Maps',false),
       */ 
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
    $(backBtn).hide();
    $(contentBox).show();
    $(filterBox).show();

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
            if(tempName!=''){view.addViewer(DG.Viewer.fromType(tempName, table))}
            clearRoot([filterBox,contentBox,descriptionBox]);
            clearRoot([root]) 
            clearRoot([recommends]);
            clearRoot([cards]);
            tempName = '';
        });

    $(dlg.root).find('.d4-dialog-contents').removeClass('ui-form');
    $(dlg.root).find('.d4-dialog-contents').removeClass('ui-panel');
    $(dlg.root).find('.d4-dialog-contents').addClass('ui-box');
    $(dlg.root).find('.d4-command-bar').append(selectioText);
    
    $(dlg.root).find('.d4-dialog-contents').css('padding','0px');
    $(dlg.root).find('.d4-dialog-header').css('border-bottom','1px solid var(--grey-2)');
    $(dlg.root).find('.d4-dialog-header').prepend(backBtn);
    
    $(dlg.root).find('.d4-command-bar').css('border-top','1px solid var(--grey-2)');

    //$(dlg.root).find('.d4-command-bar').prepend(okBtn);
    okBtn = $(dlg.root).find('.d4-command-bar > button').first();
    okBtn.addClass('ui-btn-raised');
    okBtn.text('ADD')
    $(okBtn).hide();
    /*
    okBtn.addEventListener('click',()=>{
        if(tempName!=''){view.addViewer(DG.Viewer.fromType(tempName, table))}
        clearRoot([filterBox,contentBox,descriptionBox]);
        clearRoot([root]) 
        clearRoot([recommends]);
        clearRoot([cards]);
        tempName = '';
    })
    */

    dlgTitle = $(dlg.root).find('.d4-dialog-header>.d4-dialog-title');
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

function getFilter(inputs:InputBase[], viewers:any){
    let arr = [];
    let values = [];
    for (let i in inputs){
        arr = arr.concat(inputs[i].value);
    }

    if (arr.length>0){
        for (let i in viewers){
            if (arr.indexOf(viewers[i].category) != -1){
                values.push(viewers[i].name);
            }
            if (arr.indexOf(viewers[i].type) != -1){
                values.push(viewers[i].name);
            }
        }
    }

    return values;
};

function insertSpaces(string:string) {
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

        tempName = viewer.name;
        $('.vg-selection-text').html('Selected: '+ viewer.name);

        dlgTitle.text(viewer.name);
        //show details block
        $(okBtn).show();
        root.style.border = '1px solid var(--grey-2)';
        $('.vg-selection-text').html('Selected: '+ viewer.name);

        $(backBtn).show();
        $(descriptionBox).show();
        $(contentBox).hide();
        $(filterBox).hide();
        //clearRoot([descriptionBox]);
        
        //descriptionBox = ui.box();
        let markup = ui.panel();
        let options = ui.div();
        let viewerbox = ui.box();
        viewerbox.append(DG.Viewer.fromType(viewer.name, table).root);
        let tabs = ui.box();
        tabs.append(ui.tabControl({
            'Description': ()=> markup,
            'Options': ()=> options
        }).root);
        tabs.style.paddingLeft = '20px';
        
        clearRoot([descriptionBox])
        //descriptionBox.innerHTML = '';
        descriptionBox.append(ui.splitH([viewerbox,tabs], {style:{height:'100%'}}));   

        //@ts-ignore
        let link = DG.Viewer.fromType(viewer.name, table).helpUrl;
        grok.dapi.fetchProxy('https://raw.githubusercontent.com/datagrok-ai/public/master'+link)
            .then(response => response.text())
            .then(data => {
                let res = data.replace(new RegExp('../../', 'g'), 'https://raw.githubusercontent.com/datagrok-ai/public/master/help/')
                markup.append(ui.markdown(res));
                $(markup).find('img').css('width','300px');
            });
           
          

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
        $(okBtn).show();
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