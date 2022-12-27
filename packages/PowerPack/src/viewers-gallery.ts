import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import { DataFrame, InputBase, TableView, View } from 'datagrok-api/dg';
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
    'Map',
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
    "legendVisibility": "Never",
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

let dlgTitle: any;
let dlgFooter: any;
let dlg: any;
let view: DG.TableView;

let okBtn: any;

let backBtn = ui.iconFA('arrow-left', () => {
    dlgFooter.hide();
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

    let viewers: any = {};
    const table = grok.shell.t;
    view = grok.shell.tableView(table.name);
    let root = ui.box();

    let recommend = ui.divH([], 'viewer-gallery');
    let selectioText = ui.div(['Selected: '], 'vg-selection-text');
    let cards = ui.divH([], 'viewer-gallery');

    $(recommend).css('min-height', '360px');
    $(recommend).css('overflow', 'hidden');


    let showMore = ui.link('Show more', () => {
        if ($(showMore).hasClass('vg-link-expand')) {
            $(showMore).text('Show less');
            $(showMore).addClass('vg-link-collapse');
            $(showMore).removeClass('vg-link-expand');
            console.log($(showMore))
        } else if ($(showMore).hasClass('vg-link-collapse')) {
            $(showMore).text('Show more');
            $(showMore).addClass('vg-link-expand');
            $(showMore).removeClass('vg-link-collapse');
            console.log($(showMore))
        }

        $(recommend).toggleClass('vg-expanded');

        if ($(recommend).hasClass('vg-expanded')) {
            $(recommend).css('min-height', 'auto');
            $(recommend).css('overflow', 'initial');
        } else {
            $(recommend).css('min-height', '360px');
            $(recommend).css('overflow', 'hidden');
        }
    });
    showMore.className = 'ui-link vg-link vg-link-expand';
    $(showMore).attr('data-content', '\f078');

    let search = ui.searchInput('', '', (value: String) => {
        clearRoot([recommend]);
        clearRoot([cards]);
        $(okBtn).hide();
        tempName = '';
        let filter_value = getFilter([filters_type, filters_cat], viewers);

        if (filter_value.length === 0) {
            for (let i in viewers) {
                filter_value.push(viewers[i].name)
            }
        }

        if (value != '') {
            for (let i in viewers) {
                if (filter_value.indexOf(viewers[i].name) != -1 && viewers[i].name.toLowerCase().includes(value.toLowerCase())) {
                    if (viewers[i].recommend)
                        recommend.append(render(viewers[i], table, Number(i)+1))
                    else
                        cards.append(render(viewers[i], table, Number(i)+1));
                }
            }
        }
        else {
            for (let i in viewers) {
                if (filter_value.indexOf(viewers[i].name) != -1) {
                    if (viewers[i].recommend)
                        recommend.append(render(viewers[i], table, Number(i)+1))
                    else
                        cards.append(render(viewers[i], table, Number(i)+1));
                }
            }
        }
        if ($(recommend).children().length == 0) {
            $(recommend).css('min-height', '0px');
            $(showMore).hide();
        } else {
            $(recommend).css('min-height', '360px');
            $(showMore).show();
        }
        $('.vg-selection-text').html('Selected: ' + $('.viewer-gallery').find('.fa-minus').length);
    });
    search.input.setAttribute('tabindex','-1');

    //@ts-ignore
    let filters_type = ui.multiChoiceInput('', [''], ['Viewer', 'Widget'], (value) => {
        getFilter([filters_type, filters_cat], viewers);
        search.fireChanged();
    });

    //@ts-ignore
    let filters_cat = ui.multiChoiceInput('', [''], ['Comparisons', 'Trends', 'Correlations', 'Relationships', 'Maps', 'Others'], (value) => {
        getFilter([filters_type, filters_cat], viewers);
        search.fireChanged();
    })

    //@ts-ignore
    search.input.placeholder = 'Search by name or type';

    let searchBlock = ui.block([ui.div([search.input], 'd4-search-ba')], 'vg-controls grok-gallery-search-bar');

    clearRoot([filterBox, contentBox, descriptionBox]);

    for (let i = 0; i < DG.Viewer.getViewerTypes().length; i++) {

        Object.assign(viewers, {
            [i]: {
                icon: 'grok-icon svg-icon svg-' + DG.Viewer.getViewerTypes()[i].toLowerCase().replace(/(\s)/g, '-'),
                name: insertSpaces(DG.Viewer.getViewerTypes()[i]),
                type: 'Viewer',
                recommend: false,
                category: 'Other',
                disabled: false
            }
        });
        /*
        if (i < 7) {
            viewers[i].recommend = true;
        }
        */
        if (cat_comparisons.indexOf(viewers[i].name) != -1) {
            viewers[i].category = 'Comparisons';
        } else if (cat_trends.indexOf(viewers[i].name) != -1) {
            viewers[i].category = 'Trends';
        } else if (cat_correlations.indexOf(viewers[i].name) != -1) {
            viewers[i].category = 'Correlations';
        } else if (cat_relationships.indexOf(viewers[i].name) != -1) {
            viewers[i].category = 'Relationships';
        } else if (cat_maps.indexOf(viewers[i].name) != -1) {
            viewers[i].category = 'Maps';
        } else {
            viewers[i].category = 'Others';
        }

        if (viewers[i].name.includes('widget')) {
            viewers[i].type = 'Widget';
            viewers[i].icon = 'grok-icon svg-icon svg-project';
        }

        //disable viewers
        /*
        if (i > 40) {
            viewers[i].disabled = true;
        }
        */
    }

    filterBox.append(ui.divV([
        ui.h3('Type'),
        filters_type.root,
        ui.h3('Category'),
        filters_cat.root
    ], 'vg-filter-panel'));


    contentBox.append(
        ui.divV([
            searchBlock,
            ui.h1('Recommend'),
            recommend,
            showMore,
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

    dlg = ui.dialog('Add viewer')
        .add(root
            //viewerRoot
        )
        .onOK(() => {
            if (tempName != '') { view.addViewer(DG.Viewer.fromType(tempName, table)) }
            clearRoot([filterBox, contentBox, descriptionBox]);
            clearRoot([root])
            clearRoot([recommend]);
            clearRoot([cards]);
            tempName = '';
        });

    $(dlg.root).find('.d4-dialog-contents').removeClass('ui-form');
    $(dlg.root).find('.d4-dialog-contents').removeClass('ui-panel');
    $(dlg.root).find('.d4-dialog-contents').addClass('ui-box');
    $(dlg.root).find('.d4-command-bar').append(selectioText);

    $(dlg.root).find('.d4-dialog-contents').css('padding', '0px');
    $(dlg.root).find('.d4-dialog-header').css('border-bottom', '1px solid var(--grey-2)');
    $(dlg.root).find('.d4-dialog-header').prepend(backBtn);

    $(dlg.root).find('.d4-command-bar').css('border-top', '1px solid var(--grey-2)');

    okBtn = $(dlg.root).find('.d4-command-bar > button').first();
    okBtn.addClass('ui-btn-raised');
    okBtn.text('ADD')
    $(okBtn).hide();

    dlgTitle = $(dlg.root).find('.d4-dialog-header>.d4-dialog-title');
    dlgFooter = $(dlg.root).find('.d4-dialog-footer');
    dlgFooter.hide();

    if (grok.shell.v.type == 'TableView') {
        dlg.showModal(true);
        setTimeout(function(){
            $(search.input).trigger('focus');
        }, 200);
    }

    for (let i in viewers) {
        if (viewers[i].recommend)
            recommend.append(render(viewers[i], table, Number(i)+1))
        else {
            cards.append(render(viewers[i], table, Number(i)+1));
        }
    }

    if ($(recommend).children().length == 0) {
        $(recommend).css('min-height', '0px');
        $(showMore).hide();
    } else {
        $(recommend).css('min-height', '360px');
        $(showMore).show();
    }

}

function getFilter(inputs: InputBase[], viewers: any) {
    let arr: any[] = [];
    let values = [];
    for (let i in inputs) {
        arr = arr.concat(inputs[i].value);
    }

    if (arr.length > 0) {
        for (let i in viewers) {
            if (arr.indexOf(viewers[i].category) != -1) {
                values.push(viewers[i].name);
            }
            if (arr.indexOf(viewers[i].type) != -1) {
                values.push(viewers[i].name);
            }
        }
    }

    return values;
};

function insertSpaces(string: string) {
    string = string.replace(/^(_|-)/g, '')
    string = string.replace(/(_|-)/g, ' ');
    string = string.replace(/([a-z]\s)/g, m => m.toLowerCase());
    string = string.replace(/([a-z])([A-Z])/g, '$1 $2');
    string = string.replace(/([A-Z])([A-Z][a-z])/g, '$1 $2');
    string = string.replace(/([A-Z])([a-z])/g, m => m.toLowerCase());
    string = string.replace(/^./g, m => m.toUpperCase())
    return string
}

function clearRoot(root: HTMLDivElement[]) {
    for (const i in root) {
        root[i].innerHTML = '';
    }
}

function render(viewer: any, table: DG.DataFrame, index:number) {
    let root = ui.div([]);
    root.setAttribute('tabindex', String(index));
    let icon = ui.iconFA('');
    icon.className = 'grok-icon svg-icon ' + viewer.icon;
    let label = ui.div([viewer.name], 'card-label');
    let add = ui.button(ui.iconFA('plus'), () => {
        root.click();
    });

    let details = ui.button(ui.icons.help(() => { }, 'Click to read more about ' + viewer.name), () => {
        dlgFooter.show();

        tempName = viewer.name;
        $('.vg-selection-text').html('Selected: ' + viewer.name);

        dlgTitle.text(viewer.name);
        //show details block
        $(okBtn).show();
        $('.vg-selection-text').html('Selected: ' + viewer.name);

        $(backBtn).show();
        $(descriptionBox).show();
        $(contentBox).hide();
        $(filterBox).hide();

        let markup = ui.panel();
        let options = ui.div();
        let viewerbox = ui.box();
        viewerbox.append(DG.Viewer.fromType(viewer.name, table).root);
        let tabs = ui.box();
        tabs.append(ui.tabControl({
            'Description': () => markup,
            'Options': () => options
        }).root);
        tabs.style.paddingLeft = '20px';

        clearRoot([descriptionBox])
        descriptionBox.append(ui.splitH([viewerbox, tabs], { style: { height: '100%' } }));

        //@ts-ignore
        let link = DG.Viewer.fromType(viewer.name, table).helpUrl;
        grok.dapi.fetchProxy('https://raw.githubusercontent.com/datagrok-ai/public/master' + link)
            .then(response => response.text())
            .then(data => {
                let res = data.replace(new RegExp('../../', 'g'), 'https://raw.githubusercontent.com/datagrok-ai/public/master/help/')
                markup.append(ui.markdown(res));
                $(markup).find('img').css('width', '300px');
            });
    });

    $(add).hide();

    if (viewer.recommend) {
        let viewerRoot = ui.box(DG.Viewer.fromType(viewer.name, table, viewerOptions).root);
        root.className = 'd4-item-card viewer-gallery vg-card';
        root.append(ui.block([viewerRoot, ui.divH([label, add, details])]));
    } else {
        root.className = 'd4-item-card viewer-gallery vg-card-small';
        if (viewer.disabled)
            root.classList.add('disabled')
        root.append(ui.divH([icon, label, add, details]));
    }
    root.addEventListener('click', () => {
        //view.addViewer(DG.Viewer.fromType(tempName, table))
        if (root.classList.contains('disabled') != true) {
            dockViewers(viewer.name, table, view);
            dlg.close();
        }
    });

    root.addEventListener('keydown', function(event){
        if (event.key === "Enter") {
            event.preventDefault();
            root.click();
        }
    });

    if (viewer.disabled != true) {
        root.addEventListener('mouseover', () => {
            $(add).show();
            ui.tooltip.bind(add, 'Click to add ' + viewer.name)
        });
        root.addEventListener('mouseout', () => {
            $(add).hide();
        });
    }

    return root
}

function dockViewers(viewer: string, table: DG.DataFrame, view: DG.TableView) {
    view.addViewer(DG.Viewer.fromType(viewer, table));
}
