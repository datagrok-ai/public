import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import $ from 'cash-dom';
import { filter } from 'rxjs/operators';
import { Tutorial } from "../../../tutorial";
import { _package } from '../../../package';


export class VirtualScreeningTutorial extends Tutorial {
  get name() {
    return 'Virtual Screening';
  }

  get description() {
    return 'Tutorial description in a card';
  }

  get steps() {
    return 15;
  }

  //demoTable: string = 'chem/smiles_only.csv';
  helpUrl: string = 'https://datagrok.ai/help/domains/chem/chem-curate';

  protected async _run(): Promise<void> {
    // TODO: add the dataset to demo files
    // Experimental data
    const t = await grok.data.loadTable(`${_package.webRoot}src/tracks/chem/tables/chem-tutorial-1-1.csv`);
    const v = grok.shell.addTableView(t);
    // Generated data
    const smiles = await grok.data.loadTable(`${_package.webRoot}src/tracks/chem/tables/chem-tutorial-1-2.csv`);
    grok.shell.addTableView(smiles);
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

    const addNCIcon = null;
    const addNCDlg = await this.openDialog('Open the "Add New Column" dialog', 'Add New Column', addNCIcon);
    await this.dlgInputAction(addNCDlg, 'Name the new column "pKi"', '', 'pKi');

    const pKiInfo = '<b>Ki</b> is a binding constant for each structure represented in nanomolar concentration. ' +
      'To curate the experimental activity, we change the units to molar, perform a log transformation as ' +
      'log-transformed concentrations are commonly normally distributed, and invert the values, since the ' +
      'compounds with a lower <b>Ki</b> are more active.';
    const formulaRegex = /^(9\s*-\s*Log10\(\$\{Ki\}\)|Sub\(9,\s*Log10\(\$\{Ki\}\)\))$/;
    await this.action('Transform the activity column "Ki" according to the formula "9 - Log10(${Ki})"',
      t.onColumnsAdded.pipe(filter((data) => data.args.columns.some((col: DG.Column) => {
        return col.name === 'pKi' &&
          col.tags.has(DG.TAGS.FORMULA) &&
          formulaRegex.test(col.tags[DG.TAGS.FORMULA]);
      }))), addNCIcon, pKiInfo);

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

    const computeDescriptors = async (descriptors: string[]) => {
      const chemMenu = null;
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

      await this.action('Click on the header of a column with curated structures',
        grok.events.onAccordionConstructed.pipe(filter((acc) => {
          if (acc.context instanceof DG.Column && acc.context?.name === 'curated_structures') {
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

      const groupDescription = 'To exclude <b>ExactMolWt</b>, expand the <b>Descriptors</b> group. ' +
        'The fastest way is to check the box for the entire group and unselect this particular descriptor. ' +
        'The results of computation will be added to the main dataframe and can be used to train a model.';

      await this.action('Select groups "Descriptors" (excluding "ExactMolWt"), "GraphDescriptors", ' +
        '"MolSurf", "EState VSA" and "Descriptors 3D" and press the "OK" button',
        grok.functions.onAfterRunAction.pipe(filter((call) => {
          const inputs = call.inputs.get('descriptors');
          return call.func.name === 'ChemDescriptors' &&
            descriptors.length == inputs.length &&
            descriptors.every((d) => inputs.includes(d));
        })), groupHints, groupDescription);
    };

    await computeDescriptors(descriptors);

    const pmv = await this.openViewByType('Click on "ML | Train Model..."', 'PredictiveModel');

    // UI generation delay
    await new Promise((resolve) => setTimeout(resolve, 2000));
    await this.choiceInputAction(pmv.root, `Set "Table" to "${t.name}"`, 'Table', t.name);
    await this.columnInpAction(pmv.root, 'Set "Predict" to "pKi"', 'Predict', 'pKi');

    const ignoredColumnNames = ['SMILES', 'ID', 'Ki', 'pKi', 'curated_structures'];
    const featureSelectionTip = 'Select <b>All</b> columns in the popup dialog and uncheck ' +
      `the first five columns: ${ignoredColumnNames}.`;
    await this.columnsInpAction(pmv.root, 'Set "Features" to all calculated descriptors',
      'Features', `(${descriptors.length}) ${t.columns.names()
        .filter((name: string) => !ignoredColumnNames.includes(name))
        .join(', ')}`, featureSelectionTip);

    await this.textInpAction(pmv.root, 'Increase maximum runtime to 500 seconds', 'Max runtime secs', '500');

    const modelInfo = 'Now we are ready to develop a model to predict activity. All the specified ' +
      'descriptors are now used as features and the first dataset is now a training set. Build the ' +
      'model and create an "observed versus predicted" plot to make sure that the model adequately ' +
      'predicts activity in the training set.';
    await this.buttonClickAction(pmv.root, 'Click the "Train" button', 'TRAIN', modelInfo);
    await this.contextMenuAction('Right-click on the trained model and select "Apply to | ' +
      `${t.toString()}"`, t.toString(), null, 'The result will be available in the selected ' +
      'table as a column named "Outcome".');

    grok.shell.v = v;
    const sp = await this.openPlot('scatter plot', (x) => x.type === DG.VIEWER.SCATTER_PLOT);
    const info = <{ [key: string]: any }>sp.getInfo();
    const colSelectionTip = 'Simply start typing the column name when the column list opens. ' +
      'When the column is selected, look at how close the predictions are to the actual values.';
    await this.action('Set X to "pKi"', info.xColSelector.onChanged.pipe(filter((name: string) =>
      name === 'pKi')), info.xColSelector.root);
    await this.action('Set Y to "outcome"', info.yColSelector.onChanged.pipe(filter((name: string) =>
      name === 'outcome')), info.yColSelector.root, colSelectionTip);

    const tablesPane = grok.shell.sidebar.getPane('Tables');
    const tabPaneHints = [tablesPane.header, $(tablesPane.content)
      .find('div.d4-toggle-button')
      .filter((idx, el) => Array.from(el.children).some((c) => c.textContent === smiles.name))[0]!,
    ];
    const generatedDataInfo = 'Now use the model to screen the generated structures. To predict activities for ' +
      'these structures, repeat the curation procedure for them and generate descriptors for the curated dataset.';

    await this.action(`Find "${smiles.name}" in the tables tab`, grok.events.onCurrentViewChanged.pipe(
      filter((_) => grok.shell.v.name === smiles.name && grok.shell.v.type === DG.VIEW_TYPE.TABLE_VIEW)),
      tabPaneHints, generatedDataInfo);

    await computeDescriptors(descriptors);

    const funcPane = grok.shell.sidebar.getPane('Functions');
    const funcPaneHints = [funcPane.header, $(funcPane.content)
      .find(`div.d4-toggle-button[data-view=${DG.View.MODELS}]`)[0]!];

    const pmBrowserDescription = 'This is Predictive Model Browser. Here, you can browse ' +
      'models that you trained or that were shared with you. It\'s time to apply the model ' +
      'to another dataset, which has been added to your open tables.';

    await this.openViewByType('Click on "Functions | Models" to open the Model Browser',
      DG.View.MODELS, funcPaneHints, pmBrowserDescription);
    await this.contextMenuAction('Right-click on the trained model and select "Apply to | ' +
      `${smiles.toString()}"`, smiles.toString());
  }
}
