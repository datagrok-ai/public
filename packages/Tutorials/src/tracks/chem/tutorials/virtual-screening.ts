import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import $ from 'cash-dom';
import { filter, map } from 'rxjs/operators';
import { Tutorial } from "../../../tutorial";
import { _package } from '../../../package';
import { interval } from 'rxjs';
import wu from "wu";


export class VirtualScreeningTutorial extends Tutorial {
  get name() {
    return 'Virtual Screening';
  }

  get description() {
    return 'In this tutorial, you are playing the role of an in silico researcher. ' +
    'The department of chemical synthesis has provided a number of new compounds, and biological ' +
    'laboratories measured their in vitro activities against one of the molecular targets, which ' +
    'is engaged in viral replication. Your mission is to determine if there are compounds of these ' +
    'chemical classes which could be synthesized and investigated at the next iteration, and thus reduce ' +
    'the costs of following synthesis and research by excluding less potent molecules from the list.';
  }

  get steps() {
    return 46;
  }

  demoTable: string = 'chem/tutorials/training-data.csv';
  helpUrl: string = 'https://datagrok.ai/help/domains/chem/chem-curate';

  protected async _run(): Promise<void> {
    // Save the view with experimental data
    const v = grok.shell.v;
    // Generated data
    const screeningData = await grok.data.getDemoTable('chem/tutorials/screening-data.csv');
    const screeningSetView = grok.shell.addTableView(screeningData);
    grok.shell.v = v;

    this.header.textContent = this.name;
    this.describe('<p>In this tutorial, you are playing the role of an <i>in silico</i> researcher. ' +
      'The department of chemical synthesis has provided a number of new compounds, and biological ' +
      'laboratories measured their <i>in vitro</i> activities against one of the molecular targets, which ' +
      'is engaged in viral replication.</p><p>Your mission is to determine if there are compounds of these ' +
      'chemical classes which could be synthesized and investigated at the next iteration, and thus reduce ' +
      'the costs of following synthesis and research by excluding less potent molecules from the list.</p>' +
      '<p>You have two datasets:' + 
      `<ol>
        <li>A dataset with experimental results: structural information on synthesized 
          compounds and related biological activity in the <b>Ki</b> column.</li>
        <li>A generated dataset with all structures of the same chemical classes which were not yet synthesized.</li>
      </ol></p>`);

    // this.describe(ui.link('More about chemical structures curation', this.helpUrl).outerHTML);

    const curatedMolColName = 'curated_molecule';
    const activityColName = 'pKi';
    const predictedColName = 'outcome';

    this.title('Process biological activity data');

    const addNCDlg = await this.openAddNCDialog();
    await this.dlgInputAction(addNCDlg, `Name the new column "${activityColName}"`, '', activityColName);

    const pKiInfo = '<b>Ki</b> is a binding constant for each structure represented in nanomolar concentration. ' +
      'To curate the experimental activity, we change the units to molar, perform a log transformation as ' +
      'log-transformed concentrations are commonly normally distributed, and invert the values, since the ' +
      'compounds with a lower <b>Ki</b> are more active.';
    const formulaRegex = /^(9\s*-\s*Log10\(\$\{Ki\}\)|Sub\(9,\s*Log10\(\$\{Ki\}\)\))$/;
    await this.action('Transform the column "Ki" according to the formula "9 - Log10(${Ki})" and click "OK"',
      this.t!.onColumnsAdded.pipe(filter((data) => data.args.columns.some((col: DG.Column) => {
        return col.name === activityColName &&
          col.tags.has(DG.TAGS.FORMULA) &&
          formulaRegex.test(col.tags[DG.TAGS.FORMULA]);
      }))), null, pKiInfo);

    const descriptors = ['MolWt', 'HeavyAtomMolWt', 'NumValenceElectrons', 'NumRadicalElectrons',
      'MaxPartialCharge', 'MinPartialCharge', 'MaxAbsPartialCharge', 'MinAbsPartialCharge',
      'FpDensityMorgan1', 'FpDensityMorgan2', 'FpDensityMorgan3', 'BalabanJ', 'BertzCT', 'Chi0',
      'Chi0n', 'Chi0v', 'Chi1', 'Chi1n', 'Chi1v', 'Chi2n', 'Chi2v', 'Chi3n', 'Chi3v', 'Chi4n',
      'Chi4v', 'HallKierAlpha', 'Ipc', 'Kappa1', 'Kappa2', 'Kappa3', 'LabuteASA', 'PEOE_VSA1',
      'PEOE_VSA10', 'PEOE_VSA11', 'PEOE_VSA12', 'PEOE_VSA13', 'PEOE_VSA14', 'PEOE_VSA2', 'PEOE_VSA3',
      'PEOE_VSA4', 'PEOE_VSA5', 'PEOE_VSA6', 'PEOE_VSA7', 'PEOE_VSA8', 'PEOE_VSA9', 'SMR_VSA1',
      'SMR_VSA10', 'SMR_VSA2', 'SMR_VSA3', 'SMR_VSA4', 'SMR_VSA5', 'SMR_VSA6', 'SMR_VSA7', 'SMR_VSA8',
      'SMR_VSA9', 'SlogP_VSA1', 'SlogP_VSA10', 'SlogP_VSA11', 'SlogP_VSA12', 'SlogP_VSA2', 'SlogP_VSA3',
      'SlogP_VSA4', 'SlogP_VSA5', 'SlogP_VSA6', 'SlogP_VSA7', 'SlogP_VSA8', 'SlogP_VSA9', 'TPSA',
      'EState_VSA1', 'EState_VSA10', 'EState_VSA11', 'EState_VSA2', 'EState_VSA3', 'EState_VSA4',
      'EState_VSA5', 'EState_VSA6', 'EState_VSA7', 'EState_VSA8', 'EState_VSA9', 'VSA_EState1',
      'VSA_EState10', 'VSA_EState2', 'VSA_EState3', 'VSA_EState4', 'VSA_EState5', 'VSA_EState6',
      'VSA_EState7', 'VSA_EState8', 'VSA_EState9', 'PMI1', 'PMI2', 'PMI3', 'NPR1', 'NPR2', 'RadiusOfGyration',
      'InertialShapeFactor', 'Eccentricity', 'Asphericity', 'SpherocityIndex'];

    const computeDescriptors = async (descriptors: string[], hasHistory = false) => {
      if (!hasHistory) this.title('Compute molecular descriptors');

      const chemMenu = this.getMenuItem('Chem');
      const curationInfo = 'At first glance at the provided chemical data, all compounds were extracted as salts ' +
        'with different counterions, and there might be different inconsistencies in the raw data. Let\'s curate ' +
        'the chemical structures given in the dataset. We will assume that all these compounds are present in the ' +
        'body in a neutralized form.';
      const curationDlg = await this.openDialog('Open "Chem | Curate"', 'CurateChemStructures', chemMenu, curationInfo);

      const neutralizationInfo = 'Perform "Neutralization" to remove charges.';
      await this.dlgInputAction(curationDlg, 'Check "Neutralization"', 'Neutralization', 'true', neutralizationInfo);

      const tautomerInfo = 'Perform "Tautomerization" to transform all tautomers to a unified form.';
      await this.dlgInputAction(curationDlg, 'Check "Tautomerization"', 'Tautomerization', 'true', tautomerInfo);

      const mainFragmentInfo = 'Choose "Main Fragment" to remove counterions.';
      await this.dlgInputAction(curationDlg, 'Check "Main Fragment"', 'Main Fragment', 'true', mainFragmentInfo);

      const outputComment = '';
      await this.action('Click "OK" and wait for the procedures to complete',
        grok.functions.onAfterRunAction.pipe(filter((call) => {
          return call.func.name === 'CurateChemStructures' &&
          call.inputs.get('neutralization') &&
          call.inputs.get('tautomerization') &&
          call.inputs.get('mainFragment');
        })), null, outputComment);

      let ppColumn: DG.Accordion;
      const ppDescription = 'In the property panel, you should now see all the actions ' +
        'applicable to the column under the "Actions" section. Let\'s calculate some ' +
        'molecular descriptors for the given dataset.';

      await this.action('Click on the header of the column with curated molecules',
        grok.events.onAccordionConstructed.pipe(filter((acc) => {
          if (acc.context instanceof DG.Column && acc.context?.name === curatedMolColName) {
            ppColumn = acc;
            return true;
          }
          return false;
        })), null, ppDescription);

      ppColumn!.getPane('Actions').expanded = true;
      const descriptorsLabel = $(ppColumn!.root)
        .find('label.d4-link-action')
        .filter((idx, el) => el.textContent == 'Chem | Descriptors...')[0];

      const descDlg = 'To characterize the molecules, we should calculate molecular descriptors – ' +
        'useful values that reflect the molecule\'s structure and thus have structural interpretation. ' +
        'Each chemical structure will have a facilitated representation as vector of descriptors. There ' +
        'are different types of descriptors, you can select them either by category or individually. We ' +
        'will combine the ones that describe molecular graph itself, as well as the surface of the whole ' +
        'molecule and its physico-chemical properties.';
      const descriptorDlg = await this.openDialog('Click on "Chem | Descriptors..." action in the property panel',
        'Descriptors', descriptorsLabel, descDlg);

      const groupHints = $(descriptorDlg.root)
        .find('.d4-tree-view-node')
        .filter((idx, el) => ['Descriptors', 'GraphDescriptors', 'MolSurf', 'EState VSA', 'Descriptors 3D']
          .some((group) => group === el.textContent))
        .get();

      if (hasHistory) groupHints.push($(descriptorDlg.root).find('i.fa-history.d4-command-bar-icon')[0]​!);

      const groupDescription = 'To exclude <b>ExactMolWt</b>, expand the <b>Descriptors</b> group. ' +
        'The fastest way is to check the box for the entire group and unselect this particular descriptor. ' +
        'The results of computation will be added to the main dataframe and can be used to train a model.';

      const descriptorDlgHistory = 'Find the previously entered parameters in the dialog\'s history.';

      await this.action('Select groups "Descriptors" (excluding "ExactMolWt"), "GraphDescriptors", ' +
        '"MolSurf", "EState VSA" and "Descriptors 3D" and press the "OK" button',
        grok.functions.onAfterRunAction.pipe(filter((call) => {
          const inputs = call.inputs.get('descriptors');
          return call.func.name === 'ChemDescriptors' &&
            descriptors.length == inputs.length &&
            descriptors.every((d) => inputs.includes(d));
        })), groupHints, hasHistory ? descriptorDlgHistory : groupDescription);
    };

    await computeDescriptors(descriptors);

    this.title('Train a model to predict activity based on molecule structure')
    const pmv = await this.openViewByType('Click on "ML | Train Model..."', 'PredictiveModel', this.getMenuItem('ML'));

    // UI generation delay
    await new Promise((resolve) => setTimeout(resolve, 2000));
    await this.choiceInputAction(pmv.root, `Set "Table" to "${this.t!.name}"`, 'Table', this.t!.name);
    await this.columnInpAction(pmv.root, `Set "Predict" to "${activityColName}"`, 'Predict', activityColName);

    const ignoredColumnNames = ['SMILES', 'ID', 'Ki', activityColName, curatedMolColName];
    const featureSelectionTip = 'Select <b>All</b> columns in the popup dialog and uncheck ' +
      `the first five columns: ${ignoredColumnNames}.`;
    await this.columnsInpAction(pmv.root, 'Set "Features" to all calculated descriptors',
      'Features', `(${descriptors.length}) ${this.t!.columns.names()
        .filter((name: string) => !ignoredColumnNames.includes(name))
        .join(', ')}`, featureSelectionTip);

    await this.textInpAction(pmv.root, 'Increase maximum runtime to 500 seconds', 'Max runtime secs', '500');

    const modelInfo = 'Now we are ready to develop a model to predict activity. All the specified ' +
      'descriptors are now used as features and the first dataset is now a training set. Build the ' +
      'model and create an "observed versus predicted" plot to make sure that the model adequately ' +
      'predicts activity in the training set.';
    await this.buttonClickAction(pmv.root, 'Click the "Train" button', 'TRAIN', modelInfo);
    await this.contextMenuAction('Right-click on the trained model and select "Apply to | ' +
      `${this.t!.toString()}"`, this.t!.toString(), null, 'The menu opens both from the status bar with the model ' +
      'name and from <b>Functions | Models</b>. The result will be available in the selected ' +
      `table as a column named <b>${predictedColName}</b>.`);

    grok.shell.v = v;
    const sp = await this.openPlot('scatter plot', (x) => x.type === DG.VIEWER.SCATTER_PLOT);
    const info = <{ [key: string]: any }>sp.getInfo();
    const colSelectionTip = 'Simply start typing the column name when the column list opens. ' +
      'When the column is selected, look at how close the predictions are to the actual values.';
    await this.action(`Set X to "${activityColName}"`, info.xColSelector.onChanged.pipe(filter((name: string) =>
      name === activityColName)), info.xColSelector.root);
    await this.action(`Set Y to "${predictedColName}"`, info.yColSelector.onChanged.pipe(filter((name: string) =>
      name === predictedColName)), info.yColSelector.root, colSelectionTip);

    const tablesPane = grok.shell.sidebar.getPane('Tables');
    const tabPaneHints = [tablesPane.header, $(tablesPane.content)
      .find('div.d4-toggle-button')
      .filter((idx, el) => Array.from(el.children).some((c) => c.textContent === screeningData.name))[0]!,
    ];
    const generatedDataInfo = 'Now use the model to screen the generated structures. To predict activities for ' +
      'these structures, repeat the curation procedure for them and generate descriptors for the curated dataset.';

    this.title('Virtual screening');
    this.describe('This procedure is called virtual screening — we are now trying to find ' +
      'the possible hits among structures that were not yet synthesized and tested.');

    await this.action(`Find "${screeningData.name}" in the tables tab`, grok.events.onCurrentViewChanged.pipe(
      filter((_) => grok.shell.v.name === screeningData.name && grok.shell.v.type === DG.VIEW_TYPE.TABLE_VIEW)),
      tabPaneHints, generatedDataInfo);

    await computeDescriptors(descriptors, true);

    const pmBrowserDescription = 'This is Predictive Model Browser. Here, you can browse ' +
      'models that you trained or that were shared with you. It\'s time to apply the model ' +
      'to another dataset, which has been added to your open tables.';

    await this.openViewByType('Click on "Functions | Models" to open the Model Browser',
      DG.View.MODELS, this.getSidebarHints('Functions', DG.View.MODELS), pmBrowserDescription);
    await this.contextMenuAction('Right-click on the trained model and select "Apply to | ' +
      `${screeningData.toString()}"`, screeningData.toString());

    await this.action(`Open "${screeningData.name}"`, grok.events.onCurrentViewChanged.pipe(
      filter((_) => grok.shell.v.name === screeningData.name && grok.shell.v.type === DG.VIEW_TYPE.TABLE_VIEW)),
      tabPaneHints);

    const sortInfo = `Sort the values in the <b>${predictedColName}</b> column in descending order by double-clicking ` +
      'the column header. Drag the column to the left if you want it to be closer to the columns with molecules.';
    await this.action('Sort the data to find the compounds with the highest activity',
      screeningSetView.grid.onRowsSorted, null, sortInfo);

    this.describe('Let us check that the given structures have similar activities due to same chemical groups.');

    const outcomeCol = screeningData.getCol(predictedColName);
    await this.action('Click on the most active compound', screeningData.onCurrentCellChanged.pipe(
      filter(() => screeningData.currentCol.name === curatedMolColName &&
        outcomeCol.get(screeningData.currentCell.rowIndex) === outcomeCol.max)),
      null, 'It should be in the first row after sorting. Make sure to pick the curated form.');
    const firstIndex = screeningData.currentCell.rowIndex;

    const filterDescription = 'Find the funnel icon <i class="grok-icon far fa-filter grok-icon-filter"></i> ' +
      'in the toolbox or press <b>Ctrl+F</b>. The panel should show a sketcher for better filtering on molecular data.';
    const filters = await this.openPlot('filters', (x) => x.type === DG.VIEWER.FILTERS, filterDescription);

    const filterColumns = filters.props.columnNames;
    if (!filterColumns.includes(curatedMolColName)) {
      filters.setOptions({columnNames: [curatedMolColName]});
    }

    // UI generation delay
    await new Promise((resolve) => setTimeout(resolve, 1000));
    const sketcherInput = $(filters.root).find('div.grok-sketcher-input > input')[0];
    if (!sketcherInput) return;

    await this.action('Filter the compounds by the substituents "CCl.OCc1ccccc1C"', interval(1000).pipe(
      map(() => wu(screeningData.rows.filters).toArray()),
      filter((filters: string[]) => filters.some((f) => f === `${curatedMolColName}: contains Cc1ccccc1CO.CCl`))), // TODO: check if this string replacement is ok
      sketcherInput, 'As you might have noticed, some of the most potent molecules in the data happen ' +
      'to have them. The next step is to identify the scaffolds at the core of these molecules and look ' +
      'them up in the ChEMBL database that contains bioactive molecules with drug-like properties.');

    this.describe('Predicted active compounds have common groups but different scaffolds. We should now ' +
      'estimate if these scaffolds have already been investigated. The filter panel may be closed.');

    const scaffoldColName = 'scaffolds';
    await this.dlgInputAction(await this.openAddNCDialog('Add an empty column for scaffolds (we will fill it in the next step)'),
      `Name it "${scaffoldColName}"`, '', scaffoldColName);
    await this.action('Click "OK" and drag the new column to the beginning of the grid', screeningData.onColumnsAdded.pipe(filter((data) =>
      data.args.columns.some((col: DG.Column) => col.name === scaffoldColName))));

    const scaffoldCol = screeningData.getCol(scaffoldColName);
    const mol1 = screeningData.getCol(curatedMolColName).get(firstIndex);
    let scaffold1 = 'O=C1CC2=C(N1)C=CC=C2';
    let scaffold2 = 'C1OC(c2cccs2)c2cccnc12';
    if (mol1.includes('s')) {
      // swap the scaffolds in case the predictions differ
      [scaffold1, scaffold2] = [scaffold2, scaffold1];
    }

    await this.action(`Paste the value "${scaffold1}" to the first row in "${scaffoldColName}"`,
      screeningData.onValuesChanged.pipe(filter(() => scaffoldCol.get(firstIndex) === scaffold1)));
    await this.action(`Paste the value "${scaffold2}" to the second row in "${scaffoldColName}"`,
      screeningData.onValuesChanged);
    // TODO: check according to outcomeCol.getSortedOrder()[1]
    scaffoldCol.semType = DG.SEMTYPE.MOLECULE;
    // TODO: invalidate the grid properly
    const testSetView = <DG.TableView>grok.shell.v;
    testSetView.loadLayout(testSetView.saveLayout());
    testSetView.grid.sort([predictedColName], [false]);

    let chemblPane: DG.AccordionPane;
    await this.action('Click on the first scaffold', grok.events.onAccordionConstructed.pipe(filter((acc) => {
      if (acc.context.value === scaffold1) {
        chemblPane = acc.getPane('ChEMBL search');
        return true;
      }
      return false;
    })));

    const chemblSearch = 'Expand the "ChEMBL search" pane and wait for the query to complete. ' +
      'If the search outputs a list of compounds containing this scaffold, it means that ' +
      'similar molecules have been successfully synthesized and investigated. It makes ' +
      'the further research less risky, but also less profitable due to license agreements. ' +
      'Otherwise, we can say that we have found a novel compound, which has more potential ' +
      'profits, but may involve more risks during synthesis and further research stages.';
    chemblPane!.expanded = false;
    await this.action('Search the scaffold in the ChEMBL database',
      interval(3000).pipe(filter(() => chemblPane!.expanded)),
      chemblPane!.root,
      chemblSearch);
    
    await this.action('Click on the second scaffold and run "ChEMBL search"',
      grok.events.onAccordionConstructed.pipe(filter((acc) => acc.context.value === scaffold2)));
    

    this.title('Conclusions');
    this.describe('We have performed a virtual screening which yielded two compounds with ' +
      'different scaffolds. The one is a brand-new scaffold and could possibly be used to ' +
      'create new drugs as new chemical entities with high benefits, though such approach ' +
      'have its limitations and complexities. The other is already used, as was shown in ' +
      'similarity search in ChEMBL, and it facilitates the development of its derivatives, ' +
      'but additional investigation is required to be sure that the potent hit is not ' +
      'protected by any umbrella patent.');
  }
}
