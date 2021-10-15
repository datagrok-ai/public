import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import $ from 'cash-dom';
import { filter } from 'rxjs/operators';
import { Tutorial } from "../../../tutorial";
import { _package } from '../../../package';


export class ActivityPredictionTutorial extends Tutorial {
  get name() {
    return 'Activity Prediction';
  }

  get description() {
    return 'Tutorial description in a card';
  }

  get steps() {
    return 15;
  }

  //demoTable: string = 'chem/smiles_only.csv';
  helpUrl: string = '';

  protected async _run(): Promise<void> {
    // TODO: add the dataset to demo files
    const t = await grok.data.loadTable(`${_package.webRoot}src/tracks/chem/tables/chem-tutorial-1-1.csv`);
    grok.shell.addTableView(t);

    this.header.textContent = this.name;
    this.describe('Introduction to this tutorial');

    this.describe(ui.link('More about ' + this.name, this.helpUrl).outerHTML);

    const chemMenu = null;
    const curationInfo = '';
    const curationDlg = await this.openDialog('Open "Chem | Curate"', 'CurateChemStructures', chemMenu, curationInfo);

    const neutralizationInfo = '';
    await this.dlgInputAction(curationDlg, 'Check "Neutralization"', 'Neutralization', 'true', neutralizationInfo);

    const tautomerInfo = '';
    await this.dlgInputAction(curationDlg, 'Check "Tautomerization"', 'Tautomerization', 'true', tautomerInfo);

    const mainFragmentInfo = '';
    await this.dlgInputAction(curationDlg, 'Check "Main Fragment"', 'Main Fragment', 'true', mainFragmentInfo);

    const outputComment = '';
    await this.action('Click "OK" and wait for the procedures to complete',
      grok.functions.onAfterRunAction.pipe(filter((call) => {
        return call.func.name === 'CurateChemStructures' &&
        call.inputs.get('neutralization') &&
        call.inputs.get('tautomerization') &&
        call.inputs.get('mainFragment');
      })), null, outputComment);

    const addNCIcon = null;
    const addNCDlg = await this.openDialog('Open the "Add New Column" dialog', 'Add New Column', addNCIcon);
    await this.dlgInputAction(addNCDlg, 'Name the new column "pKi"', '', 'pKi');

    const pKiInfo = '';
    const formulaRegex = /^(9\s*-\s*Log10\(\$\{Ki\}\)|Sub\(9,\s*Log10\(\$\{Ki\}\)\))$/;
    await this.action('Transform the activity column "Ki" according to the formula "9 - Log10(${Ki})"',
      t.onColumnsAdded.pipe(filter((data) => data.args.columns.some((col: DG.Column) => {
        return col.name === 'pKi' &&
          col.tags.has(DG.TAGS.FORMULA) &&
          formulaRegex.test(col.tags[DG.TAGS.FORMULA]);
      }))), addNCIcon, pKiInfo);

    let ppColumn: DG.Accordion;
    const ppDescription = 'In the property panel, you should now see all the actions ' +
      'applicable to the column under the "Actions" section. Let\'s calculate some ' +
      'molecular descriptors for the given dataset.';

    await this.action('Click on the header of a column with curated structures',
      grok.events.onAccordionConstructed.pipe(filter((acc) => {
        if (acc.context instanceof DG.Column && acc.context?.name == 'curated_structures') {
          ppColumn = acc;
          return true;
        }
        return false;
      })), null, ppDescription);

    ppColumn!.getPane('Actions').expanded = true;
    const descriptorsLabel = $(ppColumn!.root)
      .find('label.d4-link-action')
      .filter((idx, el) => el.textContent == 'Chem | Descriptors...')[0];

    const descDlg = 'The platform supports various groups of ' +
      'descriptors. You can select them either by category or individually.';
    const descriptorDlg = await this.openDialog('Click on "Chem | Descriptors..." action in the property panel',
      'Descriptors', descriptorsLabel, descDlg);

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
  }
}
