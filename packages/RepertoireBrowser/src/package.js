/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";
import scheme from "./TPP000153303.json";
import mutcodes from "./mutcodes.json";


export let _package = new DG.Package();

//name: launchBrowser
export async function launchBrowser(s) {

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

    // pulling CDR3 regions
    function CDR3(crd_choice){
        let schemeId;
        let baseH = 'darkblue';
        let baseL = 'darkred';


        if (crd_choice.value === 'chothia') {
            schemeId = NGL.ColormakerRegistry.addSelectionScheme([
                ["yellow", "25-31 and :H or 51-56 and :H or 98-106 and :H"],
                [baseH, "* and :H"],
                ["green", "23-38 and :L or 54-60 and :L or 93-101 and :L"],
                [baseL, "* and :L"]
            ]);
        } else if (crd_choice.value === 'contact') {
            schemeId = NGL.ColormakerRegistry.addSelectionScheme([
                ["yellow", "29-34 and :H or 46-58 and :H or 96-105 and :H"],
                [baseH, "* and :H"],
                ["green", "34-40 and :L or 50-59 and :L or 93-100 and :L"],
                [baseL, "* and :L"]
            ]);
        } else if (crd_choice.value === 'imgt') {
            schemeId = NGL.ColormakerRegistry.addSelectionScheme([
                ["yellow", "25-32 and :H or 50-57 and :H or 96-106 and :H"],
                [baseH, "* and :H"],
                ["green", "26-36 and :L or 54-56 and :L or 93-101 and :L"],
                [baseL, "* and :L"]
            ]);
        } else if (crd_choice.value === 'kabat') {
            schemeId = NGL.ColormakerRegistry.addSelectionScheme([
                ["yellow", "30-34 and :H or 49-65 and :H or 98-106 and :H"],
                [baseH, "* and :H"],
                ["green", "23-38 and :L or 54-60 and :L or 93-101 and :L"],
                [baseL, "* and :L"]
            ]);
        } else if (crd_choice.value === 'north') {
            schemeId = NGL.ColormakerRegistry.addSelectionScheme([
                ["yellow", "22-34 and :H or 49-58 and :H or 96-106 and :H"],
                [baseH, "* and :H"],
                ["green", "23-38 and :L or 53-60 and :L or 93-101 and :L"],
                [baseL, "* and :L"]
            ]);
        } else {
            schemeId = NGL.ColormakerRegistry.addSelectionScheme([
                [baseH, "* and :H"],
                [baseL, "* and :L"]
            ]);
        }
        return {color: schemeId};
    }

    // ngl loading
    async function loadPdb(bytes, repChoice, schemeObj) {
        stage.loadFile(bytes).then(function (o) {
            o.addRepresentation(repChoice.value, schemeObj);
            o.autoView();
        });
    }

    // sequence loading
    function loadSequecne(chain_choice, ptm_choice, ptm_prob){

        let seq;
        if (chain_choice.value === 'H') {
            seq = scheme.heavy_seq
        } else {
            seq = scheme.light_seq
        }

        let mutations = [mutcodes[ptm_choice.value]];
        let rawlist = mutationsTolist(mutcodes, scheme, chain_choice.value);
        let cl = mutationsToFeatures(seq, rawlist, mutations, ptm_prob.value);
        let gradient = cl[0];
        let features = cl[1];
        let den_gradient = cl[2];
        let den_feature = cl[3];

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
        seqEntry.addFeatures(features.map(function(ft) {
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
        seqEntry.addFeatures(den_feature.map(function(ft) {
            return {
                category: 'PTM density',
                type : 'D',
                start : ft,
                end : ft,
                text : '',
                improbable : true
            }
        }));

        applyGradient(gradient, chain_choice.value, mutations);
        applyGradient(den_gradient, chain_choice.value, ['D']);

    }

    ///// MAIN BODY ////


    let reps = ['cartoon','backbone','ball+stick','licorice','hyperball', 'surface'];
    let repChoice = ui.choiceInput('Representation', 'cartoon', reps);

    let CDR3_choice = ui.choiceInput('CDR3 Scheme', 'default', ['default','chothia','contact','imgt','kabat','north']);

    let chain_choice = ui.choiceInput('Chain', 'H',Object.keys(scheme.ptm_predictions));
    let ptm_predictions = [...new Set([...Object.keys(scheme.ptm_predictions.H), ...Object.keys(scheme.ptm_predictions.L)])];
    let ptm_choice = ui.choiceInput('PTMs', 'Phosphoserine_Phosphothreonine', ptm_predictions);

    let ptm_prob = ui.floatInput('Pr threshold', 0.2);

    repChoice.onChanged(async () => {
        $(ngl_host).empty();
        stage = new NGL.Stage(ngl_host);
        let schemeObj = CDR3(CDR3_choice);
        await loadPdb(path, repChoice, schemeObj);
    });

    CDR3_choice.onChanged(async () => {
        $(ngl_host).empty();
        stage = new NGL.Stage(ngl_host);
        let schemeObj = CDR3(CDR3_choice);
        await loadPdb(path, repChoice, schemeObj);
    });

    chain_choice.onChanged(() => {
        loadSequecne(chain_choice, ptm_choice, ptm_prob);
    });

    ptm_choice.onChanged(() => {
        loadSequecne(chain_choice, ptm_choice, ptm_prob);
    });

    ptm_prob.onChanged(() => {
        loadSequecne(chain_choice, ptm_choice, ptm_prob);
    });


    let tname = grok.shell.tableNames;
    let view = grok.shell.getTableView(tname[0]);
    let table = view.table;

    // let logger = new DG.Logger((m) => m.params['log_param'] = 'ig-repert');
    // var last_row = 0;
    // table.onCurrentRowChanged.subscribe(function () {
    //     if (table.currentRow.idx >= 0) {
    //         if (last_row !== table.currentRow.idx) { // row really changed
    //             var table_name = table.name;
    //             var row_str = table.currentRow.idx.toString();
    //             // send row selection event to the logger
    //             logger.log('row-selection', {lparam: table_name, seq_id: row_str}, 'rlog');
    //             last_row = table.currentRow.idx; // update last row
    //         }
    //     }
    // });

    let root = ui.div();
    root.appendChild(ui.h2('NGL options'));
    root.appendChild(ui.inputs([repChoice, CDR3_choice, chain_choice, ptm_choice, ptm_prob]));
    grok.shell.o = root;

    var ngl_host = ui.div([],'d4-ngl-viewer');
    ngl_host.style.backgroundColor ='black';
    view.box = true;
    view.dockManager.dock(ngl_host, 'right');
    var stage = new NGL.Stage(ngl_host);
    let path = _package.webRoot + 'pdbfiles/' + 'TPP000153303.pdb';
    let schemeObj = CDR3(CDR3_choice);
    await loadPdb(path, repChoice, schemeObj);

    let pViz_host = ui.box();
    view.dockManager.dock(pViz_host, 'down');
    var pviz = window.pviz;
    loadSequecne(chain_choice, ptm_choice, ptm_prob);

}

// let mutations = [
//     "ADA", "AI", "FR", "HYL", "HYP", "LG",
//     "MA",'MLY',"MeO","N6AL","NG", "NtG","OlG",
//     "PP","PTB","PCA","SPC",'SUMO','TO','UB'
// ];
// let mutations = Object.keys(mutcodes).filter((key) => {return (ptm_choice.value).includes(key)});