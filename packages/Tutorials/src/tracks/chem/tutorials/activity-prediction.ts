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
    return 6;
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
    const dlg = await this.openDialog('Open "Chem | Curate"', 'CurateChemStructures', chemMenu, curationInfo);

    const neutralizationInfo = '';
    await this.dlgInputAction(dlg, 'Check "Neutralization"', 'Neutralization', 'true', neutralizationInfo);

    const tautomerInfo = '';
    await this.dlgInputAction(dlg, 'Check "Tautomerization"', 'Tautomerization', 'true', tautomerInfo);

    const outputComment = '';
    await this.action('Click "OK" and wait for the procedures to complete',
      grok.functions.onAfterRunAction.pipe(filter((call) => {
        return call.func.name === 'CurateChemStructures' &&
        call.inputs.get('neutralization') &&
        call.inputs.get('tautomerization');
      })), null, outputComment);

    const addNCIcon = null;
    const pKiInfo = '';
    const formulaRegex = /^(9\s*-\s*Log10\(\$\{Ki\}\)|Sub\(9,\s*Log10\(\$\{Ki\}\)\))$/;
    await this.action('Transform the activity column into column "pKi"', t.onColumnsAdded.pipe(
      filter((data) => data.args.columns.some((col: DG.Column) => {
        return col.name === 'pKi' &&
          col.tags.has(DG.TAGS.FORMULA) &&
          formulaRegex.test(col.tags[DG.TAGS.FORMULA]);
      }))), addNCIcon, pKiInfo);
  }
}
