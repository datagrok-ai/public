/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";
import scheme from "./TPP000153303.json";
import mutcodes from "./mutcodes.json";


export let _package = new DG.Package();

//name: Repertoire Browser
//tags: app
export async function RepertoireBrowserApp() {
    let tname = grok.shell.tableNames;
    let view = grok.shell.getTableView(tname[0]);
    grok.shell.v = view;
    launchBrowser(view);
}

//name: launchBrowser
export async function launchBrowser(view) {

    // processes JSON to derive scheme names
    function extract_schemes() {
        let raw_scheme_names = Object.keys(scheme.cdr_ranges);
        let schemes_lst = ['default'];
        raw_scheme_names.forEach((str)=>{
            str = str.split('_')
            if (schemes_lst.includes(str[0]) === false) {
                schemes_lst.push(str[0]);
            }
        })
        return schemes_lst;
    }

    // palette generation
    function interpolateColor(color1, color2, factor) {
        if (arguments.length < 3) {
            factor = 0.5;
        }
        var result = color1.slice();
        for (var i = 0; i < 3; i++) {
            result[i] = Math.round(result[i] + factor * (color2[i] - color1[i]));
        }
        let hex_col = "#" + ((1 << 24) + (result[0] << 16) + (result[1] << 8) + result[2]).toString(16).slice(1);
        return hex_col;
    };

    function interpolateColors(color1, color2, steps) {
        var stepFactor = 1 / (steps - 1),
            interpolatedColorArray = [];

        color1 = color1.match(/\d+/g).map(Number);
        color2 = color2.match(/\d+/g).map(Number);

        for(var i = 0; i < steps; i++) {
            interpolatedColorArray.push(interpolateColor(color1, color2, stepFactor * i));
        }

        return interpolatedColorArray;
    }

    // convert JSON to feature vector
    function mutationsTolist(mutcodes, obj, chain) {
        let l = [];
        Object.keys(obj.ptm_predictions[chain]).forEach((mut) => {
            Object.keys(obj.ptm_predictions[chain][mut]).forEach((i) => {
                let single_mut = [];
                let single_mut_capsule = [];

                single_mut[0] = obj.ptm_predictions[chain][mut][i][0]
                single_mut[1] = mutcodes[mut];
                single_mut_capsule.push(single_mut);
                single_mut_capsule.push(obj.ptm_predictions[chain][mut][i][1]);

                l.push(single_mut_capsule);
            })
        })
        return l;
    }

    function mutationsToFeatures(seq, rawlist, mutations, prob) {
        let l_mut = [];
        let c_mut = [];
        let l_den = [];
        let c_den = new Array(seq.length).fill(1);
        let palette = interpolateColors('(255, 255, 0)','(255, 0, 0)',5);
        for (let i = 0; i < rawlist.length; i++) {

            if (rawlist[i][1] !== 100) {

                if (l_den.includes(rawlist[i][0][0]) === false) {
                    l_den.push(rawlist[i][0][0]);
                }
                c_den[l_den.indexOf(rawlist[i][0][0])] = c_den[l_den.indexOf(rawlist[i][0][0])]*(1-rawlist[i][1]);

                if (mutations.includes(rawlist[i][0][1]) && rawlist[i][1] > prob) {
                    l_mut.push(rawlist[i][0]);
                    c_mut.push(palette[Math.round(rawlist[i][1]*4)]);
                }

            }
        }

        c_den = c_den.filter((x) => {return x !== 1});
        for (let i = 0; i < c_den.length; i++) {
            c_den[i] = palette[Math.round((1-c_den[i])*4)];
        }

        return [c_mut,l_mut,c_den,l_den];
    }

    // applying gradients to sequence viewer
    function applyGradient(gradient, chain, mutations){

        let j = -1;
        mutations.forEach((mut) => {
            let selectorStr = 'g.feature.' + mut + ' rect.feature';
            let el = document.querySelectorAll(selectorStr);
            for (let i = 1; i <= el.length; i++) {
                el[i-1].style.fill = gradient[j+i];
            }
            j = j + el.length;
        })
    }

    // peritope processing
    function paratopeToList(chain_choice, paratopes) {
        let l = [];
        let c = [];
        let palette = interpolateColors('(255, 255, 255)','(255, 0, 255)',100);
        if (paratopes === true) {
            Object.keys(scheme.parapred_predictions[chain_choice]).forEach((index) => {
                l.push(index);
                c.push(palette[Math.round(scheme.parapred_predictions[chain_choice][index]*100)]);
            })
        }
        return [c,l];
    }

    // cdr track to sequence view
    function cdrToList(chain_choice, crd_choice){

        let cdr_features = [];
        if(crd_choice !== 'default') {
            if (chain_choice === 'H') {
                Object.keys(scheme.cdr_ranges).forEach((str) => {
                    if (str.includes(crd_choice + '_CDRH')) {
                        for (let i = 0; i < Object.keys(scheme.cdr_ranges[str]).length; i++) {
                            cdr_features.push([scheme.cdr_ranges[str][i][0], scheme.cdr_ranges[str][i][1]]);
                        }
                    }
                });


            } else if (chain_choice === 'L') {
                Object.keys(scheme.cdr_ranges).forEach((str) => {
                    if (str.includes(crd_choice + '_CDRL')) {
                        for (let i = 0; i < Object.keys(scheme.cdr_ranges[str]).length; i++) {
                            cdr_features.push([scheme.cdr_ranges[str][i][0], scheme.cdr_ranges[str][i][1]]);
                        }
                    }
                });

            }
        }
        return cdr_features;
    }

    // pulling CDR3 regions
    function CDR3(crd_choice, paratopes){
        let schemeId;
        let baseH = 'darkblue';
        let baseL = 'darkred';


        if (paratopes.value === true) {
            let palette = interpolateColors('(255, 255, 255)','(255, 0, 255)',100);
            let selectionScheme = [];
            Object.keys(scheme.parapred_predictions).forEach((chain) => {
                Object.keys(scheme.parapred_predictions[chain]).forEach((index) => {
                    selectionScheme.push([
                        palette[Math.round(scheme.parapred_predictions[chain][index]*100)],
                        `${index} and :${chain}`
                    ]);
                })
            })
            schemeId = NGL.ColormakerRegistry.addSelectionScheme(selectionScheme);
        } else {

            if (crd_choice.value === 'default') {
                schemeId = NGL.ColormakerRegistry.addSelectionScheme([
                    [baseH, "* and :H"],
                    [baseL, "* and :L"]
                ]);
            } else {
                let scheme_buffer = [];
                Object.keys(scheme.cdr_ranges).forEach((str) => {
                    if (str.includes(crd_choice.value + '_CDRH')) {
                        let str_buffer = ''
                        for(let i = 0; i < Object.keys(scheme.cdr_ranges[str]).length; i++) {
                            str_buffer = str_buffer + ` or ${scheme.cdr_ranges[str][i][0]}-${scheme.cdr_ranges[str][i][1]} and :H`;
                        }
                        str_buffer = str_buffer.slice(4);
                        scheme_buffer.push(["limegreen", str_buffer]);
                        scheme_buffer.push([baseH, "* and :H"]);

                    } else if (str.includes(crd_choice.value + '_CDRL')) {
                        let str_buffer = ''
                        for( let i = 0; i < Object.keys(scheme.cdr_ranges[str]).length; i++) {
                            str_buffer = str_buffer + ` or ${scheme.cdr_ranges[str][i][0]}-${scheme.cdr_ranges[str][i][1]} and :L`;
                        }
                        str_buffer = str_buffer.slice(4);
                        scheme_buffer.push(["limegreen", str_buffer]);
                        scheme_buffer.push([baseL, "* and :L"]);
                    }
                });
                schemeId = NGL.ColormakerRegistry.addSelectionScheme(scheme_buffer);
            }
        }
        return {color: schemeId};
    }

    // sidechain selection
    function sidechain_select(ptm_features, chain_choice) {
        let sidechains = '';
        ptm_features = ptm_features.flat();
        ptm_features = ptm_features.filter(function(el) {return typeof el !== 'string'});
        ptm_features = [...new Set(ptm_features)];

        for (let i=0; i < ptm_features.length; i++) {
            sidechains = sidechains + `${ptm_features[i] + 1} and :${chain_choice.value} and (not backbone or .CA or (PRO and .N))`
            if (i !== ptm_features.length - 1) {
                sidechains = sidechains + ' or ';
            }
        }

        return sidechains;
    }

    // ngl loading
    async function loadPdb(bytes, repChoice, schemeObj, sidechains = '') {
        stage.loadFile(bytes).then(function (o) {
            o.addRepresentation(repChoice.value, schemeObj);
            if (sidechains.length > 0) {
                o.addRepresentation( "ball+stick", { sele: sidechains} );
            }
            o.autoView();
        });
    }

    // sequence loading
    function loadSequence(chain_choice, ptm_choice, ptm_prob, paratopes, crd_choice){

        let seq;
        if (chain_choice.value === 'H') {
            seq = scheme.heavy_seq
        } else {
            seq = scheme.light_seq
        }

        let mutations = []
        let ptm_choices_lst = ptm_choice.value;
        ptm_choices_lst.forEach((ptm) => {
            mutations.push(mutcodes[ptm]);
        })

        let rawlist = mutationsTolist(mutcodes, scheme, chain_choice.value);
        let ml = mutationsToFeatures(seq, rawlist, mutations, ptm_prob.value);
        let pl = paratopeToList(chain_choice.value, paratopes.value);
        let cdr_features = cdrToList(chain_choice.value, crd_choice.value);
        let ptm_gradient = ml[0];
        let ptm_features = ml[1];
        let den_gradient = ml[2];
        let den_features = ml[3];
        let par_gradient = pl[0];
        let par_features = pl[1];


        let seqEntry = new pviz.SeqEntry({
            sequence : seq
        });
        new pviz.SeqEntryAnnotInteractiveView({
            model : seqEntry,
            el : pViz_host
        }).render();

        mutations.forEach((mut) => {
            pviz.FeatureDisplayer.trackHeightPerCategoryType[mut] = 1.5;
            pviz.FeatureDisplayer.setStrikeoutCategory(mut);
        });
        seqEntry.addFeatures(ptm_features.map(function(ft) {
            return {
                groupSet: 'PTMs',
                category : ft[1],
                type : ft[1],
                start : ft[0],
                end : ft[0],
                text : ft[1],
                improbable : true
            }
        }));
        seqEntry.addFeatures(par_features.map(function(pft) {
            return {
                category: 'Paratope predictions',
                type : 'P',
                start : pft,
                end : pft,
                text : '',
                improbable : true
            }
        }));
        seqEntry.addFeatures(den_features.map(function(dft) {
            return {
                category: 'PTM density',
                type : 'D',
                start : dft,
                end : dft,
                text : '',
                improbable : true
            }
        }));
        seqEntry.addFeatures(cdr_features.map(function(cft) {
            return {
                category: 'CDR region',
                type : 'CDR',
                start : cft[0],
                end : cft[1],
                text : '',
                improbable : true
            }
        }));

        applyGradient(ptm_gradient, chain_choice.value, mutations);
        applyGradient(den_gradient, chain_choice.value, ['D']);
        applyGradient(par_gradient, chain_choice.value, ['P']);

        return ptm_features;
    }

    // selection saving
    async function save_load(table, root) {

        async function saveSelectedRows(table, uniqueId, connection, fileToSave) {

            let indexes = table
                .groupBy([`${uniqueId.value}`])
                .whereRowMask(table.selection)
                .aggregate();

            let data = `${table.toString()}\n${uniqueId.stringValue}\n${indexes.col(0).toList()}`;
            await grok.dapi.files.writeAsText(`${connection}${fileToSave.value}.txt`, data);
        }

        async function loadSelectedRows(table, connection, savedFilesList) {

            let res = await grok.dapi.files.readAsText(`${connection}${savedFilesList.value}`);
            res = res.split("\n");
            let uniqueColumnName = res[1];
            let values = JSON.parse("[" + res[2] + "]");
            values = values.map((e) => parseInt(e));
            table.rows.select((row) => values.includes(row[`${uniqueColumnName}`]));

        }


        let fileToSave = ui.stringInput('FileName', 'filename');
        let connection = 'Demo:TestJobs:Files:DemoFiles/';

        let files = await grok.dapi.files.list(connection, false, '');
        files = files.map((e) => e.path);
        let savedFilesList = await ui.choiceInput('Saved Rows', ' ', files)
        // please define here you primary key column
        // let uniqueId = ui.columnInput('Unique id column', table, table.col('tenx_barcode'));
        let uniqueId = ui.stringInput('Unique id column', 'tenx_barcode');


        let saveDialog = () => {
            ui.dialog('Save rows to file')
                .add(uniqueId).add(fileToSave)
                .onOK(() => saveSelectedRows(table, uniqueId, connection, fileToSave)).show();
        };

        let loadDialog = async () => {
            let files = await grok.dapi.files.list(connection, false, '');
            files = files.map((e) => e.path);
            let savedFilesList = await ui.choiceInput('Saved Rows', ' ', files)

            ui.dialog('Load rows from file')
                .add(savedFilesList)
                .onOK(() => loadSelectedRows(table, connection, savedFilesList)).show();
        };


        let saveRowsButton = ui.button('SAVE');
        saveRowsButton.addEventListener("click", saveDialog);

        let loadRowsButton = ui.button('LOAD')
        loadRowsButton.addEventListener("click", loadDialog);

        let acc_save = ui.accordion();
        acc_save.addPane('save row selection', () => ui.divH([saveRowsButton, loadRowsButton]));
        root.append(acc_save.root);
    }

    ///// MAIN BODY ////

    let reps = ['cartoon','backbone','ball+stick','licorice','hyperball', 'surface'];
    let repChoice = ui.choiceInput('Representation', 'cartoon', reps);

    let schemes_lst = extract_schemes();
    let CDR3_choice = ui.choiceInput('CDR3 Scheme', 'default', schemes_lst);

    let chain_choice = ui.choiceInput('Chain', 'H',Object.keys(scheme.ptm_predictions));
    let ptm_predictions = [...new Set([...Object.keys(scheme.ptm_predictions.H), ...Object.keys(scheme.ptm_predictions.L)])];
    let ptm_choice = ui.multiChoiceInput('', [], ptm_predictions);

    let ptm_prob = ui.floatInput('PTM probability', 0.2);

    let paratopes = ui.boolInput('Paratopes', false);

    repChoice.onChanged(async () => {
        let sidechains = sidechain_select(ptm_features_obj.ft, chain_choice);

        $(ngl_host).empty();
        stage = new NGL.Stage(ngl_host);
        let schemeObj = CDR3(CDR3_choice, paratopes);
        await loadPdb(path, repChoice, schemeObj, sidechains);
    });

    CDR3_choice.onChanged(async () => {
        loadSequence(chain_choice, ptm_choice, ptm_prob, paratopes, CDR3_choice);
        let sidechains = sidechain_select(ptm_features_obj.ft, chain_choice);

        $(ngl_host).empty();
        stage = new NGL.Stage(ngl_host);
        let schemeObj = CDR3(CDR3_choice, paratopes);
        await loadPdb(path, repChoice, schemeObj, sidechains);
    });

    paratopes.onChanged(async () => {
        let sidechains;
        if (paratopes.value === true) {
            sidechains = '';
        } else {
            sidechains = sidechain_select(ptm_features_obj.ft, chain_choice);
        }

        $(ngl_host).empty();
        stage = new NGL.Stage(ngl_host);
        let schemeObj = CDR3(CDR3_choice, paratopes);
        await loadPdb(path, repChoice, schemeObj, sidechains);
        loadSequence(chain_choice, ptm_choice, ptm_prob, paratopes, CDR3_choice);
    });

    chain_choice.onChanged(async () => {
        ptm_features_obj.ft = loadSequence(chain_choice, ptm_choice, ptm_prob, paratopes, CDR3_choice);
        let sidechains = sidechain_select(ptm_features_obj.ft, chain_choice);

        $(ngl_host).empty();
        stage = new NGL.Stage(ngl_host);
        let schemeObj = CDR3(CDR3_choice, paratopes);
        await loadPdb(path, repChoice, schemeObj, sidechains);
    });

    ptm_choice.onChanged(async () => {
        ptm_features_obj.ft = loadSequence(chain_choice, ptm_choice, ptm_prob, paratopes, CDR3_choice);
        let sidechains = sidechain_select(ptm_features_obj.ft, chain_choice);

        $(ngl_host).empty();
        stage = new NGL.Stage(ngl_host);
        let schemeObj = CDR3(CDR3_choice, paratopes);
        await loadPdb(path, repChoice, schemeObj, sidechains);
    });

    ptm_prob.onChanged(async () => {
        ptm_features_obj.ft = loadSequence(chain_choice, ptm_choice, ptm_prob, paratopes, CDR3_choice);
        let sidechains = sidechain_select(ptm_features_obj.ft, chain_choice);

        $(ngl_host).empty();
        stage = new NGL.Stage(ngl_host);
        let schemeObj = CDR3(CDR3_choice, paratopes);
        await loadPdb(path, repChoice, schemeObj, sidechains);
    });


    let table = view.table;

    // region Logging ---------------
    let logger = new DG.Logger((m) => m.params['log_param'] = 'ig-repert');
    table.onCurrentRowChanged.subscribe(function () {
         if (table.currentRow.idx >= 0) {
             var row_str = table.currentRow.idx.toString();
             logger.log('row-change', {lparam: table.name, seq_id: row_str}, 'rlog');
         }
    });
    table.onCurrentCellChanged.subscribe(function () {
        if (table.currentRow.idx >= 0) {
            // grok.shell.info(`Cell: ${t.currentCell.rowIndex}, ${t.currentCell.column.name}`));
            logger.log('cell-change', {lparam: table.name, col: table.currentCell.column.name}, 'rlog');
        }
    });
    // endregion

    let root = ui.div();

    root.appendChild(ui.h1('Repertoire viewer options'));

    root.appendChild(ui.h3('NGL settings'));
    root.appendChild(ui.inputs([repChoice, CDR3_choice]));

    root.appendChild(ui.h3('Sequence settings'));
    root.appendChild(ui.inputs([chain_choice, paratopes, ptm_prob]));
    let acc_ptm = ui.accordion();
    acc_ptm.addPane('ptm list', () => ui.inputs([ptm_choice]));
    root.append(acc_ptm.root);

    root.appendChild(ui.h3('Save/Load'));
    await save_load(table, root);
    grok.shell.o = root;

    var ngl_host = ui.div([],'d4-ngl-viewer');
    ngl_host.style.backgroundColor ='black';
    view.box = true;
    view.dockManager.dock(ngl_host, 'right');
    var stage = new NGL.Stage(ngl_host);
    let path = _package.webRoot + 'pdbfiles/' + 'TPP000153303.pdb';
    let schemeObj = CDR3(CDR3_choice, paratopes);
    await loadPdb(path, repChoice, schemeObj);

    let pViz_host = ui.box();
    view.dockManager.dock(pViz_host, 'down');
    var pviz = window.pviz;
    let ptm_features_obj = {};
    ptm_features_obj.ft = loadSequence(chain_choice, ptm_choice, ptm_prob, paratopes ,CDR3_choice);

}