/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export let _package = new DG.Package();

function paramSelector(method) {
    let methodparams;
    if (method === 'CalcBasicStats') {
        let p1 = ui.choiceInput('Weighting method', '',['by clonotype frequency','unweighted']);
        methodparams = [p1];
    }

    if (method === 'CalcSegmentUsage') {
        let p1 = ui.choiceInput('Weighting method', '',['by clonotype frequency','unweighted']);
        methodparams = [p1];
    }

    if (method === 'CalcSpectratype') {
        let p1 = ui.choiceInput('Weighting method', '',['by clonotype frequency','unweighted']);
        let p2 = ui.choiceInput('CDR3 sequence', '',['amino acid', 'nucleotide']);
        methodparams = [p1, p2];

    }
    return methodparams;
}

function shCallCompose(method, methodparams) {
    let defaultCall = 'java -Xmx10G -jar pathToJar routine';
    defaultCall = defaultCall.replace('routine', method);
    let resultsFile = 'out/0.'

    if (method === 'CalcBasicStats') {
        resultsFile = resultsFile + 'basicstats.txt'
        if (methodparams[0].value === 'unweighted') {
            defaultCall = defaultCall + ' -u';
        }
    }

    if (method === 'CalcSegmentUsage') {
        resultsFile = resultsFile + 'segments.wt.V.txt'
        if (methodparams[0].value === 'unweighted') {
            defaultCall = defaultCall + ' -u';
            resultsFile = resultsFile.replace('wt','unwt')
        }
    }

    if (method === 'CalcSpectratype') {
        resultsFile = resultsFile + 'spectratype.nt.wt.txt'
        if (methodparams[0].value === 'unweighted') {
            defaultCall = defaultCall + ' -u';
            resultsFile = resultsFile.replace('wt','unwt')
        }
        if (methodparams[1].value === 'amino acid') {
            defaultCall = defaultCall + ' -a';
            resultsFile = resultsFile.replace('nt','aa')
        }
    }
    defaultCall = defaultCall + ' samplePath out/0'
    return [defaultCall, resultsFile];
}

async function analyze(data, shCall) {
    let metadata = await grok.functions.call('VDJtools:vdjwrapper',
        {
            'clonotypeDf': data,
            'shcall':shCall[0],
            'resultsfile':shCall[1]
        })
    return metadata;
}

//name: vdjtools
export async function vdjtools() {

    //Create dialogue window
    let v = ui.dialog('VDJ Tools');

    let tableName = ui.choiceInput('Table', null, grok.shell.tableNames);
    let dataTable = grok.shell.tableByName(tableName.value);
    let analysisRoutine = ui.choiceInput('Routine','',['CalcBasicStats', 'CalcSegmentUsage', 'CalcSpectratype']);
    let methodparams;

    let container0 = ui.div();
    let inputs0 = ui.inputs([tableName, analysisRoutine]);
    container0.appendChild(inputs0);

    let container_acc = ui.div();
    let acc = ui.accordion();

    analysisRoutine.onChanged(function () {
        $(container_acc).empty()
        methodparams = paramSelector(analysisRoutine.value);
        container_acc.appendChild(ui.inputs(methodparams));
    });
    acc.addPane('parameters',() => container_acc);


    v.add(container0);
    v.add(acc).onOK(async ()=>{

        let shCall = shCallCompose(analysisRoutine.value, methodparams)
        grok.shell.info('call '+ shCall[0]);
        grok.shell.info('output ' + shCall[1]);
        let metadata = await analyze(dataTable, shCall);
        grok.shell.addTableView(metadata);

    }).show();

}
