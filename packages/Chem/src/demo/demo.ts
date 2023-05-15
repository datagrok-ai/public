import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { closeAllAccordionPanes, demoScaffold, getAccordionPane, openMoleculeDataset, openSketcher, scrollTable } from '../utils/demo-utils';
import { DemoScript } from '@datagrok-libraries/tutorials/src/demo-script';
import { awaitCheck, delay } from '@datagrok-libraries/utils/src/test';
import { _importSdf } from '../open-chem/sdf-importer';
import { _package, chemSpaceTopMenu } from '../package';
import { rGroupAnalysis } from '../analysis/r-group-analysis';
import { chemSpace, getEmbeddingColsNames } from '../analysis/chem-space';
import { CLIFFS_DF_NAME, activityCliffsIdx, getActivityCliffs } from '@datagrok-libraries/ml/src/viewers/activity-cliffs';
import { getSimilaritiesMarix } from '../utils/similarity-utils';
import { createPropPanelElement, createTooltipElement } from '../analysis/activity-cliffs';
import { ScaffoldTreeViewer } from '../widgets/scaffold-tree';

export async function _demoChemOverview(): Promise<void> {

    const sketcherType = DG.chem.currentSketcherType;
    DG.chem.currentSketcherType = 'OpenChemLib';

    const demoScript = new DemoScript('Demo', 'Overview of Cheminformatics functionality');
    let table: DG.DataFrame;
    let tv: DG.TableView;
    let propPanel: Element;
    let canvas: HTMLCanvasElement;
    demoScript
        .step('Loading table', async () => {
            tv = await openMoleculeDataset('sar-small.csv');
            grok.shell.windows.showHelp = false;
            table = tv.dataFrame;
        }, { description: 'Load dataset with molecule columns', delay: 3000 })
        .step('Molecule properties', async () => {
            await delay(1000);
            grok.shell.windows.showHelp = false; //for some reason help panel appears again, need to hide it
            propPanel = document.getElementsByClassName('grok-entity-prop-panel')[0];
            closeAllAccordionPanes(propPanel!);
            const structurePaneContent = getAccordionPane('Structure', propPanel!);
            getAccordionPane('3D Structure', structurePaneContent!);
            const biologyPaneContent = getAccordionPane('Biology', propPanel!);
            getAccordionPane('Toxicity', biologyPaneContent!);
            await delay(3000);
            grok.shell.windows.showHelp = false;
            table.currentRowIdx = 5;
            grok.shell.windows.showHelp = false;
            await delay(3000);
            table.currentRowIdx = 3;
            grok.shell.windows.showHelp = true;
        }, { description: 'Molecules properties are re-calculating when changing current molecule', delay: 3000 })
        .step('Fast rendering', async () => {
            await delay(1000);
            canvas = tv.grid.root.getElementsByTagName('canvas')[2];
            await scrollTable(canvas, 300, 15, 100);
        }, { description: 'Molecules are rendered immediately when scrolling dataset', delay: 2000 })
        .step('Filtering', async () => {
            await delay(1000);
            const filters = tv.getFiltersGroup();
            await delay(1000);
            const sketcherDlg = await openSketcher(filters.root, 'sketch-link');
            const sketcherInput = sketcherDlg!.getElementsByClassName('grok-sketcher-input')[0]?.children[0] as HTMLInputElement;
            sketcherInput.value = 'C1CCCCC1';
            await delay(1000);
            sketcherInput.dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter' }));
            await delay(1000);
            await scrollTable(canvas, 200, 10, 200);
            await delay(1000);
            Array.from(sketcherDlg!.getElementsByTagName('span')).find(el => el.textContent === 'CANCEL')?.click();
            delay(500);
            filters.close();
        }, { description: 'Filtering dataset by substructure', delay: 2000 })
        .step('Aligning to scaffold', async () => {
            await delay(1000);
            grok.shell.o = tv.dataFrame.col('smiles');
            await delay(2000);
            closeAllAccordionPanes(propPanel!);
            const chemistryPaneContent = getAccordionPane('Chemistry', propPanel!);
            const renderingPaneContent = getAccordionPane('Rendering', chemistryPaneContent!) as HTMLElement;
            await delay(1000);
            const scaffoldSketcher = await openSketcher(renderingPaneContent, 'sketch-link');
            const scaffoldSketcherInput = scaffoldSketcher!.getElementsByClassName('grok-sketcher-input')[0]?.children[0] as HTMLInputElement;

            let dT = null;
            try { dT = new DataTransfer(); } catch (e) { }
            const evt = new ClipboardEvent('paste', { clipboardData: dT });
            evt.clipboardData!.setData('text/plain', demoScaffold);
            scaffoldSketcherInput.value = demoScaffold;
            await delay(100);
            scaffoldSketcherInput.dispatchEvent(evt);
            await delay(1000);
            await scrollTable(canvas, 200, 10, 200);
            await delay(1000);
            Array.from(scaffoldSketcher!.getElementsByTagName('span')).find(el => el.textContent === 'CANCEL')?.click();
        }, { description: 'Aligning structures by scaffold', delay: 1000 })
        .step('Final', async () => {
            console.log('Finished');
            DG.chem.currentSketcherType = sketcherType;
        })
        .start();
}


export async function _demoSimilaritySearch(): Promise<void> {

    const demoScript = new DemoScript('Demo', 'Searching for molecules most similar to target molecule');
    let table: DG.DataFrame;
    let tv: DG.TableView;
    demoScript
        .step('Loading table', async () => {
            tv = await openMoleculeDataset('sar-small.csv');
            table = tv.dataFrame;
        }, { description: 'Load dataset with molecule columns', delay: 2000 })
        .step('Adding viewer', async () => {
            await delay(1000);
            const similarityViewer = tv.addViewer('Chem Similarity Search');
            grok.shell.o = similarityViewer;
        }, { description: 'Open similarity search viewer. Selected molecule becomes target.', delay: 2000 })
        .step('Changing target molecule', async () => {
            table.currentRowIdx = 2;
            await delay(3000);
            table.currentRowIdx = 10;
            await delay(3000);
            table.currentRowIdx = 3;
        }, { description: 'Fast similarity search re-calculating when changing current molecule', delay: 3000 })
        .step('Final', async () => console.log('Finished'))
        .start();
}


export async function _demoRgroupAnalysis(): Promise<void> {

    const demoScript = new DemoScript('Demo', 'Performing R Group Analysis');
    let table: DG.DataFrame;
    let tv: DG.TableView;
    let sketcherInput: HTMLInputElement;
    let sketcher: Element;
    demoScript
        .step('Loading table', async () => {
            tv = await openMoleculeDataset('sar-small.csv');
            table = tv.dataFrame;
        }, { description: 'Load dataset with molecule columns', delay: 2000 })
        .step('Opening R Group Analysis viewer', async () => {
            await delay(1000);
            rGroupAnalysis(table.col('smiles')!);
            await delay(2000);
            sketcher = document.getElementsByClassName('d4-dialog')[0];
            sketcherInput = sketcher!.getElementsByClassName('grok-sketcher-input')[0]?.children[0] as HTMLInputElement;
            sketcherInput.value = 'O=C1CN=C(c2ccccc2N1)C3CCCCC3';
            sketcherInput.dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter' }));
        }, { description: 'Open R Group Analysis viewer and enter scaffold structure', delay: 2000 })
        .step('Running R Group Analysis', async () => {
            await delay(1000);
            Array.from(sketcher!.getElementsByTagName('span')).find(el => el.textContent === 'OK')?.click();
            await awaitCheck(() => {
                for (let viewer of tv.viewers) {
                    if (viewer.type === DG.VIEWER.TRELLIS_PLOT)
                        return true;
                }
                return false;
            },
                'r group analysis has not been loaded', 10000);
        }, { description: 'Trellis plot is created from R Group Analysis results', delay: 2000 })
        .step('Final', async () => console.log('Finished'))
        .start();
}


export async function _demoActivityCliffs(): Promise<void> {

    const demoScript = new DemoScript('Demo', 'Searching similar structures with significant activity difference');
    let table: DG.DataFrame;
    let tv: DG.TableView;
    let scatterPlot: DG.Viewer;
    demoScript
        .step('Loading table', async () => {
            tv = await openMoleculeDataset('activity_cliffs.csv');
            table = tv.dataFrame;
        }, { description: 'Load dataset with molecule and activity columns', delay: 2000 })
        .step('Running activity cliffs analysis', async () => {
            const molecules = table.col('smiles')!
            const progressBar = DG.TaskBarProgressIndicator.create(`Activity cliffs running...`);
            const axesNames = getEmbeddingColsNames(table);
            scatterPlot = await getActivityCliffs(table, molecules, null as any, axesNames, 'Activity cliffs', table.col('Activity')!, 80, 'Tanimoto',
                't-SNE', DG.SEMTYPE.MOLECULE, { 'units': molecules.tags['units'] }, chemSpace, getSimilaritiesMarix,
                createTooltipElement, createPropPanelElement, undefined, undefined, 0.5);
            progressBar.close();
            await delay(1000);
        }, { description: 'Results are shown on a scatter plot', delay: 2000 })
        .step('Opening table with cliffs', async () => {
            await delay(1000);
            (Array.from(scatterPlot!.root.children)
                .filter((it) => it.className === 'ui-btn ui-btn-ok scatter_plot_link cliffs_grid')[0] as HTMLElement).click();
            await delay(1000);
        }, { description: 'Detected cliffs are available in a separate table', delay: 2000 })
        .step('Selecting cliffs', async () => {
            await delay(1000);
            let cliffsGrid: DG.Viewer | null = null;
            for (const i of tv.viewers) {
                if (i.dataFrame.name === `${CLIFFS_DF_NAME}${activityCliffsIdx}`)
                    cliffsGrid = i;
            }
            cliffsGrid!.dataFrame.currentRowIdx = 0;
            await delay(3000);
            cliffsGrid!.dataFrame.currentRowIdx = 1;
            await delay(3000);
            cliffsGrid!.dataFrame.currentRowIdx = 2;
        }, { description: 'When you select a cliff scatter plot is zoomed to that exact cliff', delay: 3000 })
        .step('Final', async () => console.log('Finished'))
        .start();
}


export async function _demoDatabases(): Promise<void> {
    const ids = ['O=C(C)Oc1ccccc1C(=O)O', 'NC1=NC=NC2=C1N=CN2', 'CC(=O)O'];
    const sketcherType = DG.chem.currentSketcherType;
    DG.chem.currentSketcherType = 'OpenChemLib';

    const properties = [
        {
            "name": "pattern",
            "type": "string",
            "semType": "Molecule"
        },
        {
            "name": "threshold",
            "type": DG.TYPE.FLOAT
        },
    ];

    let props = properties.map((p) => DG.Property.fromOptions(p))

    let object = {
        pattern: '',
        threshold: 0.5,
    };

    const form = ui.input.form(object, props);
    form.classList.add('ui-form-condensed');
    form.style.minWidth = '0px';
    form.style.width = '130px';

    const queryName = ui.divText('Similarity search with threshold', { style: { width: '140px', fontWeight: 'bold' } });
    const runButton = ui.bigButton('RUN', () => { });
    runButton.style.width = '130px';
    const queryDiv = ui.divV([
        queryName,
        form,
        runButton
    ]);

    const gridDiv = ui.box(null, { style: { width: '100%', height: '100%', marginLeft: '30px' } });

    const totalDiv = ui.divH([
        queryDiv,
        gridDiv
    ], { style: { height: '100%', width: '100%' } });

    const loading = (isLoading: boolean) => {
        ui.setUpdateIndicator(gridDiv, isLoading);
    }

    const loadNewQuery = async (id: string) => {
        const parent = document.getElementsByClassName('ui-form-condensed')[0] as HTMLElement;
        const sketcherDlg = await openSketcher(parent, 'd4-input-molecule-canvas-host');
        const sketcherInput = sketcherDlg!.getElementsByClassName('grok-sketcher-input')[0]?.children[0] as HTMLInputElement;
        sketcherInput.value = id;
        await delay(1000);
        sketcherInput.dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter' }));
        await delay(1000);
        Array.from(sketcherDlg!.getElementsByTagName('span')).find(el => el.textContent === 'OK')?.click();
        runButton.classList.add('chem-demo-button-pushed');
        await delay(1000);
        runButton.classList.remove('chem-demo-button-pushed');
        ui.empty(gridDiv);
        loading(true);
    }

    const loadQueryResults = async (t: DG.DataFrame) => {
        await grok.data.detectSemanticTypes(t);
        const grid = t.plot.grid().root;
        loading(false);
        gridDiv.append(grid);
        await delay(1500);
    }

    const view = grok.shell.newView('Databases', [totalDiv]);
    view.helpUrl = `${_package.webRoot}/README.md`

    setTimeout(async () => {
        for (const id of ids) {
            await loadNewQuery(id);
            const t = await grok.data.query("Chembl:patternSimilaritySearchWithThreshold", { 'pattern': id, 'threshold': '0.5' })
            await loadQueryResults(t);
        }
        DG.chem.currentSketcherType = sketcherType;
    }, 500);

}

export async function _demoDatabases2(): Promise<void> { //Databases integration in property panel
    const table = _importSdf(await _package.files.readAsBytes('mol1K.sdf'))[0];
    grok.shell.windows.showProperties = true;
    grok.shell.windows.showHelp = false;

    //1. Open table
    const tv = grok.shell.addTableView(table);
    await delay(2000);

    //2. Open tabs on property panel
    const propPanel = document.getElementsByClassName('grok-entity-prop-panel')[0];
    closeAllAccordionPanes(propPanel!);
    const databasesPaneContent = getAccordionPane('Databases', propPanel!);
    getAccordionPane('ChEMBL (Internal) Substructure Search', databasesPaneContent!) as HTMLElement;
    getAccordionPane('ChEMBL (Internal) Similarity Search', databasesPaneContent!) as HTMLElement;
    await delay(3000);
    table.currentRowIdx = 2;
    await delay(3000);
    table.currentRowIdx = 5;
}



export async function _demoDatabases3(): Promise<void> {
    const ids = ['CHEMBL1827', 'CHEMBL1829', 'CHEMBL1830'];
    const query = `SELECT m.chembl_id AS compound_chembl_id, s.canonical_smiles, act.standard_type, act.standard_value
    FROM compound_structures s, molecule_dictionary m, compound_records r, docs d, activities act, assays a, target_dictionary t
    WHERE s.molregno     = m.molregno
    AND m.molregno       = r.molregno
    AND r.record_id      = act.record_id
    AND r.doc_id         = d.doc_id
    AND act.assay_id     = a.assay_id
    AND a.tid            = t.tid
    AND act.standard_type = 'IC50'
    AND t.chembl_id      = '~id~';`

    const connection = await grok.functions.eval('Chembl:Chembl');
    const queryPanel = ui.box();
    const gridDiv = ui.div();
    const scatterPlot = ui.div();
    const barchart = ui.div();

    const totalDiv = ui.splitV([
        queryPanel,
        gridDiv,
        ui.splitH([
            scatterPlot,
            barchart
        ])
    ], { style: { height: '100%', width: '100%' } });

    const loading = (isLoading: boolean) => {
        ui.setUpdateIndicator(gridDiv, isLoading);
        ui.setUpdateIndicator(scatterPlot, isLoading);
        ui.setUpdateIndicator(barchart, isLoading);
    }

    const loadNewQuery = (id: string) => {
        ui.empty(queryPanel);
        const queryDiv = ui.textInput(`Compound activity details for target = ${id}`, query.replace(`~id~`, id));
        queryDiv.input.style.height = '100%';
        queryPanel.append(queryDiv.root);
        const dBQuery = connection!.query('', query.replace(`~id~`, id));
        dBQuery.adHoc = true;
        ui.empty(gridDiv);
        ui.empty(scatterPlot);
        ui.empty(barchart);
        loading(true);
        return dBQuery;
    }

    const loadQueryResults = async (t: DG.DataFrame) => {
        await grok.data.detectSemanticTypes(t);
        const grid = t.plot.grid().root;
        grid.style.width = '100%';
        loading(false);
        gridDiv.append(grid);
        barchart.append(t.plot.bar().root);
        scatterPlot.append(t.plot.scatter().root);
        await delay(1500);
    }


    grok.shell.newView('Databases', [totalDiv]);

    setTimeout(async () => {
        for (const id of ids) {
            const t = await loadNewQuery(id).executeTable();
            await loadQueryResults(t);
        }
    }, 500);

}



export async function _demoDatabases4(): Promise<void> {
    const ids = ['CHEMBL1827', 'CHEMBL1829', 'CHEMBL1830'];
    const query = `--name: compound activity details for target
    --connection: Chembl
    --input: string target = ~target~
    --meta.cache: true
    --meta.localCache: true
    --meta.invalidate: 0 0 * ? * * *
    SELECT m.chembl_id AS compound_chembl_id, s.canonical_smiles, act.standard_type, act.standard_value
    FROM compound_structures s, molecule_dictionary m, compound_records r, docs d, activities act, assays a, target_dictionary t
    WHERE s.molregno     = m.molregno
    AND m.molregno       = r.molregno
    AND r.record_id      = act.record_id
    AND r.doc_id         = d.doc_id
    AND act.assay_id     = a.assay_id
    AND a.tid            = t.tid
    AND act.standard_type = 'IC50'
    AND t.chembl_id      = @target
    limit 500;
    --end`

    const connection: DG.DataConnection = await grok.functions.eval('Chembl:Chembl');
    let scatterPlot: DG.Viewer;

    const test = async (array: string[]) => {
        const id = array.shift();
        const dBQuery = connection!.query('', query.replace('~target~', `"${id}"`));
        const editor = await dBQuery.prepare().getEditor();
        ui.dialog('compound activity details for target').add(editor).onOK(async () => {
            ui.setUpdateIndicator(tv.root, true);
            let data = await dBQuery.apply({ target: id });
            ui.setUpdateIndicator(tv.root, false);
            data.name = tv.dataFrame.name;
            await grok.data.detectSemanticTypes(data);
            tv.dataFrame = data;
            if (scatterPlot)
                scatterPlot.close();
            await delay(1000);
            scatterPlot = (await chemSpaceTopMenu(data, data.col('canonical_smiles')!, 't-SNE', 'Tanimoto', true))!
            await delay(1000);
            if (array.length)
                test(array)
        }).show().history(() => { }, () => { });

        await delay(1000);
        const dlgFooter = document!.getElementsByClassName('d4-dialog-footer')[0] as HTMLInputElement;
        Array.from(dlgFooter!.getElementsByTagName('span')).find(el => el.textContent === 'OK')?.click();
    }

    test(ids);

    const tv = grok.shell.addTableView(DG.DataFrame.create());
}



export async function _demoDatabases5(): Promise<void> {
    const ids = ['CHEMBL1827', 'CHEMBL1829', 'CHEMBL1830'];
    const query = `SELECT m.chembl_id AS compound_chembl_id,
    s.canonical_smiles,
    r.compound_key,
    coalesce(d.pubmed_id::text, d.doi) AS pubmed_id_or_doi,
    a.description                   AS assay_description,   act.standard_type,
    act.standard_relation,
    act.standard_value,
    act.standard_units,
    act.activity_comment
    FROM compound_structures s, molecule_dictionary m, compound_records r, docs d, activities act, assays a, target_dictionary t
    WHERE s.molregno     = m.molregno
    AND m.molregno       = r.molregno
    AND r.record_id      = act.record_id
    AND r.doc_id         = d.doc_id
    AND act.assay_id     = a.assay_id
    AND a.tid            = t.tid
    AND act.standard_type = 'IC50'
    AND act.standard_relation = '='
    AND act.standard_units = 'nM'
    AND t.chembl_id      = '~id~'
    limit 500;`

    const properties = [
        {
            "name": "target",
            "type": "string",
        },
        {
            "name": "standardType",
            "type": "string"
        },
        {
            "name": "standardRelation",
            "type": "string"
        },
        {
            "name": "standardUnits",
            "type": "string"
        },
        {
            "name": "limit",
            "type": DG.TYPE.FLOAT
        }
    ];

    let props = properties.map((p) => DG.Property.fromOptions(p))

    let object = {
        target: 'CHEMBL1827',
        standardType: 'IC50',
        standardRelation: '=',
        standardUnits: 'nM',
        limit: 500
    };

    const form = ui.input.form(object, props);
    form.classList.add('ui-form-condensed');
    form.style.minWidth = '0px';
    form.style.width = '130px';

    const queryName = ui.divText('Compound details for target', { style: { width: '140px', fontWeight: 'bold' } });
    const runButton = ui.bigButton('RUN', () => { });
    runButton.style.width = '130px';
    const queryDiv = ui.divV([
        queryName,
        form,
        runButton
    ]);

    const gridDiv = ui.box(null, { style: { width: '100%', height: '100%', marginLeft: '30px' } });
    const scatterPlot = ui.box();
    const barchart = ui.box(null, { style: { marginLeft: '50px' } });

    const totalDiv = ui.splitV([
        ui.divH([
            queryDiv,
            gridDiv
        ]),
        ui.splitH([
            scatterPlot,
            barchart
        ])
    ], { style: { height: '100%', width: '100%' } });

    const loading = (isLoading: boolean) => {
        ui.setUpdateIndicator(gridDiv, isLoading);
        ui.setUpdateIndicator(scatterPlot, isLoading);
        ui.setUpdateIndicator(barchart, isLoading);
    }

    const connection: DG.DataConnection = await grok.functions.eval('Chembl:Chembl');

    const loadNewQuery = async (id: string) => {
        const idInput = form.children[0].getElementsByClassName('ui-input-editor')[0];
        //@ts-ignore
        idInput.value = id;
        idInput.classList.add('chem-demo-chembl-id-selected');
        await delay(1000);
        runButton.classList.add('chem-demo-button-pushed');
        await delay(1000);
        idInput.classList.remove('chem-demo-chembl-id-selected');
        runButton.classList.remove('chem-demo-button-pushed');
        const dBQuery = connection!.query('', query.replace(`~id~`, id));
        dBQuery.adHoc = true;
        ui.empty(gridDiv);
        ui.empty(scatterPlot);
        ui.empty(barchart);
        loading(true);
        return dBQuery;
    }

    const loadQueryResults = async (t: DG.DataFrame) => {
        await grok.data.detectSemanticTypes(t);
        const grid = t.plot.grid().root;
        loading(false);
        gridDiv.append(grid);
        barchart.append(t.plot.bar().root);
        scatterPlot.append(t.plot.scatter().root);
        await delay(1500);
    }

    grok.shell.newView('Databases', [totalDiv]);

    setTimeout(async () => {
        for (const id of ids) {
            const t = await (await loadNewQuery(id)).executeTable();
            await loadQueryResults(t);
        }
    }, 500);

}

export async function _demoScaffoldTree(): Promise<void> {
    const tv: DG.TableView = await openMoleculeDataset('mol1K.csv');
    grok.shell.windows.showHelp = true;
    grok.shell.windows.help.showHelp('/help/domains/chem/scaffold-tree');
    const table: DG.DataFrame = tv.dataFrame;
    const tree = await _package.files.readAsText('scaffold-tree.json');
    const viewer = new ScaffoldTreeViewer();
    viewer.autoGenerate = false;
    viewer.dataFrame = table;
    viewer.size = 'small';
    await viewer.loadTreeStr(tree);
    tv.dockManager.dock(viewer, DG.DOCK_TYPE.LEFT);
}
