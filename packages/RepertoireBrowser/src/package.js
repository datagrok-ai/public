/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";
import json from "./TPP000153303.json";
import mutcodes from "./mutcodes.json";

export let _package = new DG.Package();

//name: Repertoire Browser
//tags: app
export async function RepertoireBrowserApp() {
    // let loaded;
    // if (loaded == undefined) {
    //     let tname = grok.shell.tableNames;
    //     let view = null;
    //     if (tname === null || tname.length === 0) {
    //         let df = (await grok.functions.eval('OpenServerFile("Dskatov:RepertoireBrowser/RepertoireBrowserSample.csv")'))[0];
    //         view = grok.shell.addTableView(df);
    //     } else {
    //         view = grok.shell.getTableView(tname[0]);
    //     }
    //     grok.shell.v = view;
    //     await launchBrowser(view);
    //     loaded = true;
    // }
    let tname = grok.shell.tableNames;
    let view = grok.shell.getTableView(tname[0]);
    grok.shell.v = view;
    launchBrowser(view);
}

//name: launchBrowser
export async function launchBrowser(view) {

    // processes JSON to derive scheme names
    function extract_schemes() {
        let raw_scheme_names = Object.keys(json.cdr_ranges);
        let schemes_lst = ['default'];
        raw_scheme_names.forEach((str) => {
            str = str.split('_')
            if (schemes_lst.includes(str[0]) === false) {
                schemes_lst.push(str[0]);
            }
        })
        return schemes_lst;
    }


    // ---- pViz ----
    // color interpolation
    function interpolateColors(color1, color2, steps) {

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

        var stepFactor = 1 / (steps - 1),
            interpolatedColorArray = [];

        color1 = color1.match(/\d+/g).map(Number);
        color2 = color2.match(/\d+/g).map(Number);

        for(var i = 0; i < steps; i++) {
            interpolatedColorArray.push(interpolateColor(color1, color2, stepFactor * i));
        }

        return interpolatedColorArray;
    }

    // mapping objects for sequence rendering
    function ptmMapping(ptm_choices, prob) {

        let ptmMap = {}
        let chains = Object.keys(pVizParams.seq);
        chains.forEach((chain) => {
            let ptm_feature_map = [];
            let ptm_color_obj = {};
            let ptm_el_obj = {};
            let ptm_prob_obj = {};
            let palette = interpolateColors('(255, 255, 0)','(255, 0, 0)',5);

            ptm_choices.forEach(ptm => {
                let ptm_array = json.ptm_predictions[chain][ptm];
                if (ptm_array !== undefined) {

                    let ptm_color_arr = [];
                    let ptm_el_arr = [];
                    let ptm_prob_arr = [];
                    ptm_array.forEach(point => {
                        if(point[1] > prob) {

                            ptm_feature_map.push({
                                groupSet: 'PTMs',
                                category : mutcodes[ptm],
                                type : mutcodes[ptm],
                                start : point[0],
                                end : point[0],
                                text : mutcodes[ptm],
                                improbable : true
                            })
                            ptm_color_arr.push(palette[Math.round(point[1]*4)])
                            ptm_el_arr.push(point[0]);
                            ptm_prob_arr.push(point[1]);
                        }
                    })
                    if (ptm_color_arr.length > 0) {
                        ptm_color_obj[mutcodes[ptm]] = ptm_color_arr;
                        ptm_el_obj[mutcodes[ptm]] = ptm_el_arr;
                        ptm_prob_obj[mutcodes[ptm]] = ptm_prob_arr;
                    }
                }
            })
            ptmMap[chain] = {ptm_feature_map: ptm_feature_map, ptm_color_obj: ptm_color_obj,
                ptm_el_obj: ptm_el_obj, ptm_prob_obj: ptm_prob_obj};
        })

        return(ptmMap);
    }

    function ptmDenMapping() {

        let denMap = {}
        let chains = Object.keys(pVizParams.seq);
        chains.forEach((chain) => {
            let den_feature_map = [];
            let den_color_arr = new Array(pVizParams.seq[chain].length).fill(1);
            let palette = interpolateColors('(255, 255, 0)', '(255, 0, 0)', 5);

            Object.values(json.ptm_predictions[chain]).forEach((ptm_array) => {
                ptm_array.forEach(point => {
                    if (!(den_feature_map.includes(point[0]))) {
                        den_feature_map.push(point[0]);
                    }
                    den_color_arr[point[0]] = den_color_arr[point[0]] * (1 - point[1])
                })
            })

            den_feature_map.sort((a, b) => a - b);
            den_feature_map = den_feature_map.map(function (ft) {
                return {
                    groupSet: 'PTM density',
                    category: '',
                    type: 'D',
                    start: ft,
                    end: ft,
                    text: '',
                    improbable: true
                }
            });

            den_color_arr = den_color_arr.filter((x) => {
                return x !== 1
            });
            for (let i = 0; i < den_color_arr.length; i++) {
                den_color_arr[i] = palette[Math.round((1 - den_color_arr[i]) * 4)];
            }

            denMap[chain] = {den_feature_map: den_feature_map, den_color_obj: {'D': den_color_arr}}
        })


        return (denMap)
    }

    function paratopeMapping() {

        let parMap = {}
        let chains = Object.keys(pVizParams.seq);
        chains.forEach((chain) => {
            let par_feature_map = [];
            let par_color_arr = [];
            let palette = interpolateColors('(255, 255, 255)', '(255, 0, 255)', 100);

            Object.keys(json.parapred_predictions[chain]).forEach((index) => {
                par_feature_map.push({
                    category: 'Paratope predictions',
                    type: 'P',
                    start: index,
                    end: index,
                    text: '',
                    improbable: true
                })
                par_color_arr.push(palette[Math.round(json.parapred_predictions[chain][index] * 100)]);
            })

            parMap[chain] = {par_feature_map: par_feature_map, par_color_obj: {'P': par_color_arr}}
        })

        return(parMap)
    }

    function cdrMapping(cdr_scheme) {

        let cdrMap = {}
        let chains = Object.keys(pVizParams.seq);
        chains.forEach((chain) => {
            let cdr_feature_map = [];
            let cdr_ranges = json.cdr_ranges[cdr_scheme + '_CDR' + chain + '_ranges']
            if (cdr_scheme !== 'default') {
                Object.values(cdr_ranges).forEach((range) => {
                    cdr_feature_map.push({
                        category: 'CDR region',
                        type: 'CDR',
                        start: range[0],
                        end: range[1],
                        text: '',
                        improbable: true
                    })
                })
            }
            cdrMap[chain] = {cdr_feature_map: cdr_feature_map}
        })

        return (cdrMap)
    }

    function applyGradient(gradient_obj){
        Object.keys(gradient_obj).forEach((ptm_track) => {
            let selectorStr = 'g.feature.' + ptm_track + ' rect.feature';
            let el = document.querySelectorAll(selectorStr);
            for (let i = 0; i < el.length; i++) {
                el[i].style.fill = gradient_obj[ptm_track][i];
            }
        })
    }

    // main sequence rendering func
    async function loadSequence(host, chain, paratopes) {

        if( $(host).width() !== 0) {
            let seq = pVizParams.seq[chain]

            let seqEntry = new pviz.SeqEntry({
                sequence: seq
            });
            new pviz.SeqEntryAnnotInteractiveView({
                model: seqEntry,
                collapsible: true,
                el: host
            }).render();

            let mod_codes = Object.keys(pVizParams.ptmMap[chain].ptm_color_obj);

            mod_codes.forEach((mod) => {
                pviz.FeatureDisplayer.trackHeightPerCategoryType[mod] = 1.5;
                pviz.FeatureDisplayer.setStrikeoutCategory(mod);
            });

            
            let switchObj = {'H':{}, 'L':{}}
            pviz.FeatureDisplayer.addClickCallback (mod_codes, async function(ft) {
                let selectorStr = 'g.feature.' + ft.category + ' rect.feature';
                let el = document.querySelectorAll(selectorStr);
                let el_lst = pVizParams.ptmMap[chain].ptm_el_obj[ft.category];

                let sidechains = `${ft.start + 1} and :${chain} and (not backbone or .CA or (PRO and .N))`
                let r;
                if (switchObj[chain][ft.start] === undefined) {
                    r = stage.compList[0].addRepresentation("ball+stick", {sele: sidechains});
                    switchObj[chain][ft.start] = {};
                    switchObj[chain][ft.start]['state'] = false;
                    switchObj[chain][ft.start]['rep'] = r
                    el[el_lst.indexOf(ft.start)].style.fill = 'black';
                } else {
                    r = switchObj[chain][ft.start]['rep'];
                    r.setVisibility(switchObj[chain][ft.start]['state']);
                    if (switchObj[chain][ft.start]['state'] === false) {
                        el[el_lst.indexOf(ft.start)].style.fill = pVizParams.ptmMap[chain].ptm_color_obj[ft.category][el_lst.indexOf(ft.start)];
                    } else {
                        el[el_lst.indexOf(ft.start)].style.fill = 'black';
                    }
                    switchObj[chain][ft.start]['state'] = !switchObj[chain][ft.start]['state']
                }

                console.log(switchObj);
            })

            pviz.FeatureDisplayer.addMouseoverCallback(mod_codes, async function(ft) {
                let selectorStr = 'g.feature.' + ft.category + ' rect.feature';
                let el = document.querySelectorAll(selectorStr);
                let el_lst = pVizParams.ptmMap[chain].ptm_el_obj[ft.category];
                let prob_lst = pVizParams.ptmMap[chain].ptm_prob_obj[ft.category]
                el = el[el_lst.indexOf(ft.start)];
                let prob =  prob_lst[el_lst.indexOf(ft.start)];

                console.log(el.getBoundingClientRect().top);
                ui.tooltip.show(
                    ui.span([`${ft.category}: Pr ~${prob.toFixed(2)}`]),
                    el.getBoundingClientRect().left + 10,
                    el.getBoundingClientRect().top + 10
                );
            }).addMouseoutCallback(mod_codes, function(ft) {
                ui.tooltip.hide();
            }) ;

            if (paratopes === true) {
                seqEntry.addFeatures(pVizParams.parMap[chain].par_feature_map);
            }
            seqEntry.addFeatures(pVizParams.ptmMap[chain].ptm_feature_map);
            seqEntry.addFeatures(pVizParams.denMap[chain].den_feature_map);
            seqEntry.addFeatures(pVizParams.cdrMap[chain].cdr_feature_map);
            applyGradient(pVizParams.ptmMap[chain].ptm_color_obj);
            applyGradient(pVizParams.denMap[chain].den_color_obj);
            applyGradient(pVizParams.parMap[chain].par_color_obj);

        }
    }


    // ---- MSA ----
    function msaRender(m, msa_fasta, gff_annots = null,
                       gff_aa_annots1 = null,
                       gff_aa_annots2 = null) {
        const seqs = msa.io.fasta.parse(msa_fasta);
        m.seqs.reset(seqs);
        m.seqs.removeAllFeatures();
        if (gff_aa_annots1) {
            const features = gffParser.parseSeqs(gff_aa_annots1);
            m.seqs.addFeatures(features);
        }
        if (gff_aa_annots2) {
            const features = gffParser.parseSeqs(gff_aa_annots2);
            m.seqs.addFeatures(features);
        }
        if (gff_annots) {
            const features = gffParser.parseSeqs(gff_annots);
            m.seqs.addFeatures(features);
        }
        m.render();
    }

    function gffAnnotator(seqid, source, type, start, end, score, strand, phase, attributes) {
        return `${seqid}\t${source}\t${type}\t${start}\t${end}\t${score}\t${strand}\t${phase}\t${attributes}\n`;
    }

    // Translate AA sequence into a GFF formatted annotations
    function getAA_gff(molId, aaStr) {
        let gff = "##gff-version 3\n";
        let pos = 2;
        for (i = 0; i < aaStr.length; i++) {
            let aaBase = aaStr[i];
            let line = `${molId} . p	${pos} 	${pos}	.	.  +  Name=${aaBase};Color=white\n`;
            gff = gff + line;
            pos = pos + 3;
        }
        return gff;
    }

    // Format AA column values into GFF annotation
    function makeColGFF(col, rowIdx, seqId) {
        let gff = null;
        if (col) {
            const seq = col.get(rowIdx);
            gff = getAA_gff(seqId, seq);
        }
        return gff;
    }

    function drawAlignments() {
        const idx = table.currentRow.idx;
        if (seqHeavyCol && germHeavyCol) {
            const seqsH = `>seq_align_heavy\n${seqHeavyCol.get(idx)}\n>germline_align_heavy\n${germHeavyCol.get(idx)}\n`;

            const gffAnnotsH = '##gff-version 3\n' + (
                (vStartHeavy && vEndHeavy) ? gffAnnotator('germline_align_heavy', '.', 'gene',
                    vStartHeavy.get(idx), vEndHeavy.get(idx), '.', '+', '.', 'Name=V region;Color=violet') : '') + (
                (dStartHeavy && dEndHeavy) ? gffAnnotator('germline_align_heavy', '.', 'gene',
                    dStartHeavy.get(idx), dEndHeavy.get(idx), '.', '+', '.', 'Name=D region;Color=green') : '') + (
                (jStartHeavy && jEndHeavy) ? gffAnnotator('germline_align_heavy', '.', 'gene',
                    jStartHeavy.get(idx), jEndHeavy.get(idx), '.', '+', '.', 'Name=J region;Color=blue') : '');

            let gffAASeq = makeColGFF(seqAlignHeavyAACol, idx, "seq_align_heavy");
            let gffAAGerm = makeColGFF(germAlignHeavyAACol, idx, "germline_align_heavy");

            msaRender(msaH, seqsH,
                gffAnnotsH.length > 16 ? gffAnnotsH : null,
                gffAASeq.length > 16   ? gffAASeq : null,
                gffAAGerm.length > 16  ? gffAAGerm : null
            );
        }
        if (seqLightCol && germLightCol) {
            const seqsL = `>seq_align_light\n${seqLightCol.get(idx)}\n>germline_align_light\n${germLightCol.get(idx)}\n`;

            const gffAnnotsL = '##gff-version 3\n' + (
                (vStartLight && vEndLight) ? gffAnnotator('germline_align_light', '.', 'gene',
                    vStartLight.get(idx), vEndLight.get(idx), '.', '+', '.', 'Name=V region;Color=violet') : '') + (
                (dStartLight && dEndLight) ? gffAnnotator('germline_align_light', '.', 'gene',
                    dStartLight.get(idx), dEndLight.get(idx), '.', '+', '.', 'Name=D region;Color=green') : '') + (
                (jStartLight && jEndLight) ? gffAnnotator('germline_align_light', '.', 'gene',
                    jStartLight.get(idx), jEndLight.get(idx), '.', '+', '.', 'Name=J region;Color=blue') : '');

            let gffAASeq = makeColGFF(seqAlignLightAACol, idx, "seq_align_light");
            let gffAAGerm = makeColGFF(germAlignLightAACol, idx, "germline_align_light");

            msaRender(msaL, seqsL,
                gffAnnotsL.length > 16 ? gffAnnotsL : null,
                gffAASeq.length > 16   ? gffAASeq : null,
                gffAAGerm.length > 16  ? gffAAGerm : null
            );
        }
    }

    function drawGroupAlignments(groups, tab = 'heavy', max = 5) {
        let groupsShown = 0;
        const tabHost = (tab === 'heavy') ? msa_group_host_H : msa_group_host_L;
        $(tabHost).empty();
        for (let groupData of Object.values(groups)) {
            if (groupsShown >= max) return;
            let seqs = '';
            for (let i = 0, n = groupData.seqIds.length; i < n; i++) {
                seqs += `>${groupData.seqIds[i]}_seq\n${groupData.alignedSequences[i]}\n`;
                seqs += `>${groupData.seqIds[i]}_germ\n${groupData.alignedGermlines[i]}\n`;
            }
            const groupDiv = ui.div();
            tabHost.appendChild(groupDiv);
            msaOpts.el = groupDiv;
            const m = new msa.msa(msaOpts);
            msaRender(m, seqs);
            groupsShown += 1;
        }
      }
      
    function selectGroups() {
        const indexes = table.selection.getSelectedIndexes();
        const heavyGroups = {};
        const lightGroups = {};
        for (let i of indexes) {
            const alH = germAlignHeavyAACol.get(i);
            const seqIdH = heavyIdCol.get(i);
            const seqH = seqHeavyCol.get(i);
            const germH = germHeavyCol.get(i);
            if (alH in heavyGroups) {
                heavyGroups[alH].aaAlignments.push(alH);
                heavyGroups[alH].alignedSequences.push(seqH);
                heavyGroups[alH].alignedGermlines.push(germH);
                heavyGroups[alH].seqIds.push(seqIdH);
            } else {
                heavyGroups[alH] = {
                    aaAlignments: [alH],
                    alignedSequences: [seqH],
                    alignedGermlines: [germH],
                    seqIds: [seqIdH],
                };
            }
            const alL = germAlignLightAACol.get(i);
            const seqIdL = lightIdCol.get(i);
            const seqL = seqLightCol.get(i);
            const germL = germLightCol.get(i);
            if (alH in lightGroups) {
                lightGroups[alL].aaAlignments.push(alL);
                lightGroups[alL].alignedSequences.push(seqL);
                lightGroups[alL].alignedGermlines.push(germL);
                lightGroups[alL].seqIds.push(seqIdL);
            } else {
                lightGroups[alL] = {
                    aaAlignments: [alL],
                    alignedSequences: [seqL],
                    alignedGermlines: [germL],
                    seqIds: [seqIdL],
                };
            }
        }
        drawGroupAlignments(heavyGroups, 'heavy', 2);
        drawGroupAlignments(lightGroups, 'light', 2);
    }


    // ---- NGL ----
    // create a color scheme for CDR3 regions
    function CDR3(cdr_scheme, paratopes) {
        let schemeId;
        let baseH = 'darkblue';
        let baseL = 'darkred';

        if (paratopes.value === true) {
            let palette = interpolateColors('(255, 255, 255)', '(255, 0, 255)', 100);
            let selectionScheme = [];
            Object.keys(json.parapred_predictions).forEach((chain) => {
                Object.keys(json.parapred_predictions[chain]).forEach((index) => {
                    selectionScheme.push([
                        palette[Math.round(json.parapred_predictions[chain][index] * 100)],
                        `${index} and :${chain}`
                    ]);
                })
            })
            schemeId = NGL.ColormakerRegistry.addSelectionScheme(selectionScheme);
        } else {
            if (cdr_scheme.value === 'default') {
                schemeId = NGL.ColormakerRegistry.addSelectionScheme([
                    [baseH, "* and :H"],
                    [baseL, "* and :L"]
                ]);
            } else {
                let scheme_buffer = [];
                Object.keys(json.cdr_ranges).forEach((str) => {
                    if (str.includes(cdr_scheme.value + '_CDRH')) {
                        let str_buffer = ''
                        for (let i = 0; i < Object.keys(json.cdr_ranges[str]).length; i++) {
                            str_buffer = str_buffer + ` or ${json.cdr_ranges[str][i][0]}-${json.cdr_ranges[str][i][1]} and :H`;
                        }
                        str_buffer = str_buffer.slice(4);
                        scheme_buffer.push(["limegreen", str_buffer]);
                        scheme_buffer.push([baseH, "* and :H"]);

                    } else if (str.includes(cdr_scheme.value + '_CDRL')) {
                        let str_buffer = ''
                        for (let i = 0; i < Object.keys(json.cdr_ranges[str]).length; i++) {
                            str_buffer = str_buffer + ` or ${json.cdr_ranges[str][i][0]}-${json.cdr_ranges[str][i][1]} and :L`;
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

    // load the 3D model
    async function loadPdb(bytes, repChoice, schemeObj, sidechains = '') {
        stage.loadFile(bytes).then(function (o) {
            o.addRepresentation(repChoice.value, schemeObj);
            // if (sidechains.length > 0) {
            //     o.addRepresentation("ball+stick", {sele: sidechains});
            //     // let r = stage.compList[104].addRepresentation("ball+stick", {sele: sidechains});
            //     // r.setVisibility(false);
            // }
            o.autoView();
        });
    }


    // ---- Resizing ----
    function nglResize(host, stage) {
        let canvas = host.querySelector('canvas');

        function resize() {
            canvas.width = Math.floor(host.clientWidth * window.devicePixelRatio);
            canvas.height = Math.floor(host.clientHeight * window.devicePixelRatio);
            stage.handleResize();
        }

        ui.onSizeChanged(host).subscribe((_) => resize());
        resize();
    }

    async function pvizResize(host, chain) {
        ui.onSizeChanged(host).subscribe(async (_) => {
            await loadSequence(host, chain, paratopes.value)
        });
    }

    function setDockSize(node, nodeContent) {
        let nodeContentHeight = 0;
        let rootNodeHeight = view.dockManager.rootNode.container.containerElement.clientHeight;
        let newHeight = 0;
        newHeight = $("#feature-viewer").outerHeight(true) + 55;
        newHeight = 1/(rootNodeHeight/newHeight);
        // newHeight = Math.ceil(newHeight*100)/100;

        return view.dockManager.dock(nodeContent, 'down', node, 'Sequence', newHeight.toFixed(2));
    }


    // ---- Save/Load ----
    // selection saving
    async function save_load(table, acc) {

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
            let values = res[2].split(',');
            console.log(values);
            table.rows.select((row) => values.includes(row[`${uniqueColumnName}`]));

        }


        let fileToSave = ui.stringInput('FileName', 'filename');
        let connection = 'Demo:TestJobs:Files:DemoFiles/';

        let files = await grok.dapi.files.list(connection, false, '');
        files = files.map((e) => e.path);
        let savedFilesList = await ui.choiceInput('Saved Rows', ' ', files)
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

        acc.addPane('Save/Load', () => ui.divH([saveRowsButton, loadRowsButton]));
    }


    ///// MAIN BODY ////

    let table = view.table;
    // tweak the App page properties
    {
        let windows = grok.shell.windows;
        windows.showProperties = false;
        windows.showHelp = false;
        windows.showConsole = false;
    }

    let reps = ['cartoon', 'backbone', 'ball+stick', 'licorice', 'hyperball', 'surface'];
    let repChoice = ui.choiceInput('Representation', 'cartoon', reps);

    let schemes_lst = extract_schemes();
    let cdr_scheme = ui.choiceInput('CDR3 Scheme', 'default', schemes_lst);

    let ptm_predictions = [...new Set([...Object.keys(json.ptm_predictions.H), ...Object.keys(json.ptm_predictions.L)])];
    let ptm_choices = ui.multiChoiceInput('', [], ptm_predictions);

    let ptm_prob = ui.floatInput('PTM probability', 0.2);

    let paratopes = ui.boolInput('Paratopes', false);

    repChoice.onChanged(async () => {

        stage.removeAllComponents();
        let schemeObj = CDR3(cdr_scheme, paratopes);
        await loadPdb(path, repChoice, schemeObj);
    });

    cdr_scheme.onChanged(async () => {

        stage.removeAllComponents();
        let schemeObj = CDR3(cdr_scheme, paratopes);
        await loadPdb(path, repChoice, schemeObj);

        pVizParams.cdrMap = cdrMapping(cdr_scheme.value)
        await loadSequence(pViz_host_H, 'H', paratopes.value)
        await loadSequence(pViz_host_L, 'L', paratopes.value)

        setDockSize(ngl_node, sequence_tabs);
    });

    paratopes.onChanged(async () => {

        stage.removeAllComponents();
        let schemeObj = CDR3(cdr_scheme, paratopes);
        await loadPdb(path, repChoice, schemeObj);
        await loadSequence(pViz_host_H, 'H', paratopes.value)
        await loadSequence(pViz_host_L, 'L', paratopes.value)

        setDockSize(ngl_node, sequence_tabs);
    });

    ptm_choices.onChanged(async () => {

        pVizParams.ptmMap = ptmMapping(ptm_choices.value, ptm_prob.value)
        console.log(pVizParams);
        await loadSequence(pViz_host_H, 'H', paratopes.value)
        await loadSequence(pViz_host_L, 'L', paratopes.value)

        setDockSize(ngl_node, sequence_tabs);
    });

    ptm_prob.onChanged(async () => {

        pVizParams.ptmMap = ptmMapping(ptm_choices.value, ptm_prob.value)
        await loadSequence(pViz_host_H, 'H', paratopes.value)
        await loadSequence(pViz_host_L, 'L', paratopes.value)

        setDockSize(ngl_node, sequence_tabs);
    });


    let seqHeavyCol = table.col('sequence_alignment_heavy');
    let seqLightCol = table.col('sequence_alignment_light');
    let germHeavyCol = table.col('germline_alignment_heavy');
    let germLightCol = table.col('germline_alignment_light');

    let seqAlignHeavyAACol = table.col('sequence_alignment_aa_heavy');
    let seqAlignLightAACol = table.col('sequence_alignment_aa_light');
    let germAlignHeavyAACol = table.col('germline_alignment_aa_heavy');
    let germAlignLightAACol = table.col('germline_alignment_aa_light');

    let heavyIdCol = table.col('tenx_contig_id_heavy');
    let lightIdCol = table.col('tenx_contig_id_light');

    let vStartHeavy = table.col('v_alignment_start_heavy');
    let dStartHeavy = table.col('d_alignment_start_heavy');
    let jStartHeavy = table.col('j_alignment_start_heavy');
    let vEndHeavy = table.col('v_alignment_end_heavy');
    let dEndHeavy = table.col('d_alignment_end_heavy');
    let jEndHeavy = table.col('j_alignment_end_heavy');

    let vStartLight = table.col('v_alignment_start_light');
    let dStartLight = table.col('d_alignment_start_light');
    let jStartLight = table.col('j_alignment_start_light');
    let vEndLight = table.col('v_alignment_end_light');
    let dEndLight = table.col('d_alignment_end_light');
    let jEndLight = table.col('j_alignment_end_light');


    let root = ui.div();
    let acc_options = ui.accordion();
    acc_options.addPane('3D model', () => ui.inputs([repChoice, cdr_scheme]));
    acc_options.addPane('Sequence', () => ui.inputs([paratopes, ptm_prob]));
    acc_options.addPane('PTMs', () => ui.inputs([ptm_choices]));
    await save_load(table, acc_options)
    root.append(acc_options.root);


    var ngl_host = ui.div([],'d4-ngl-viewer');
    ngl_host.style.backgroundColor ='black';
    view.box = true;
    var stage = new NGL.Stage(ngl_host);
    let path = _package.webRoot + 'pdbfiles/' + 'TPP000153303.pdb';
    let schemeObj = CDR3(cdr_scheme, paratopes);
    await loadPdb(path, repChoice, schemeObj);


    let pViz_host_L = ui.box();
    let pViz_host_H = ui.box();
    let msa_host_L = ui.box();
    let msa_host_H = ui.box();
    let msa_group_host_L = ui.divV();
    let msa_group_host_H = ui.divV();
    let sequence_tabs =
        ui.tabControl({
            'HEAVY': pViz_host_H,
            'LIGHT': pViz_host_L,
            'MSA HEAVY': msa_host_H,
            'MSA LIGHT': msa_host_L,
            'GRP MSA HEAVY': msa_group_host_H,
            'GRP MSA LIGHT': msa_group_host_L,
        }).root;


    const msaOpts = {
        el: msa_host_L,
        vis: {
            conserv: false,
            overviewbox: false,
            consensus: true,
            seqlogo: true,
            scaleslider: false,
        },
        conf: {
            dropImport: true
        },
        bootstrapMenu: false,
        zoomer: {
        boxRectHeight: 1,
        boxRectWidth: 1,
        labelNameLength: 110,
        labelFontsize: 10,
        labelIdLength: 20
        }
    };
    const msaL = new msa.msa(msaOpts);
    msaOpts.el = msa_host_H;
    const msaH = new msa.msa(msaOpts);
    const gffParser = msa.io.gff;
    DG.debounce(table.onSelectionChanged, 200).subscribe(selectGroups);
    DG.debounce(table.onCurrentRowChanged, 200).subscribe(drawAlignments);
    drawAlignments();


    var pviz = window.pviz;
    let pVizParams = {};
    pVizParams.seq = {'H':json.heavy_seq, 'L': json.light_seq}
    pVizParams.ptmMap = ptmMapping(ptm_choices.value, ptm_prob.value)
    pVizParams.denMap = ptmDenMapping()
    pVizParams.parMap = paratopeMapping()
    pVizParams.cdrMap = cdrMapping(cdr_scheme.value)

    let panel_node = view.dockManager.dock(root, 'right', null, 'NGL');
    let ngl_node = view.dockManager.dock(ngl_host, 'left', panel_node, 'NGL');
    let sequence_node = view.dockManager.dock(sequence_tabs, 'down', ngl_node, 'Sequence', 0.225);

    await loadSequence(pViz_host_H, 'H', paratopes.value);


    // Start
    // function svgResize() {
    //     let svgParent = document.querySelector('#feature-viewer');
    //     let stopHandle  = ui.tools.handleResize(svgParent, (width, height) => {
    //
    //         console.log(width + height);
    //         svgParent = document.querySelector('#feature-viewer');
    //         if (width + height === 0) {
    //             stopHandle();
    //             svgResize();
    //         }
    //
    //         setDockSize(ngl_node, sequence_tabs);
    //     })
    // }
    // svgResize()






    nglResize(ngl_host, stage);
    pvizResize(pViz_host_H, 'H');
    pvizResize(pViz_host_L, 'L');

    setDockSize(ngl_node, sequence_tabs);

    // msaResize(msa_host_L);
    // msaResize(msa_host_H);
    // $(sequence_tabs).children().css("height","100%");
}