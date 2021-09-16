import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import { fromEvent } from 'rxjs';
import { filter, map } from 'rxjs/operators';
import { Tutorial } from "../../../tutorial";


export class PredictiveModelingTutorial extends Tutorial {
  get name() { return 'Predictive Modeling'; }
  get description() {
    return 'Predictive modeling uses statistics to predict outcomes';
  }

  protected async _run() {

    /** Train model actions */
    const trainModel = async (method: string): Promise<void> => {
      let pmv: DG.View;

      await this.action('Click on "ML | Train Model..." to open a dialog for training models',
        grok.events.onViewAdded.pipe(filter((v: DG.View) => {
          if (v.type === 'PredictiveModel') {
            pmv = v;
            return true;
          }
          return false;
        })));

      // UI generation delay
      await new Promise((resolve) => setTimeout(resolve, 1500));

      this.describe('In this view, you can set parameters to train a model (press F1 to get ' +
        'help on the Predictive Modeling plugin). In the next steps, we will train, apply, ' +
        'analyze performance, and compare models to predict sex by age, weight and height.');

      const viewInputAction = async (instructions: string, caption: string, value: string) => {
        let inputRoot: HTMLDivElement;
        let select: HTMLSelectElement;
        $(pmv.root).find('.ui-input-root .ui-input-label span').each((idx, el) => {
          if (el.innerText == caption) {
            inputRoot = el.parentElement?.parentElement as HTMLDivElement;
            select = $(inputRoot).find('select')[0] as HTMLSelectElement;
          }
        });
        const source = fromEvent(select!, 'change');
        await this.action(instructions,
          source.pipe(map((event: any) => select.value), filter((v: string) => v === value)),
          inputRoot!);
      }

      // TODO: transform to actions
      //await viewInputAction('Set "Table" to "demog"', 'Table', 'demog');
      this.describe('Set "Predict" to "SEX"');
      this.describe('Set "Features" to ["AGE", "WEIGHT", "HEIGHT"]');

      this.describe('You will now see that the column names are highlighted in red. ' +
        'Click on the question mark to the right of the input to see a warning message. ' +
        'It says that some of the columns contain missing values, and suggests two ways to address this.');

      let dlg: DG.Dialog;
      await this.action("Let's use missing values imputation.",
        grok.events.onDialogShown.pipe(filter((dialog: any) => {
          if (dialog.title === 'Missing Values Imputation') {
            dlg = dialog;
            return true;
          }
          return false;
      })));

      await this.dlgInputAction(dlg!, 'Set "Table" to "demog"', 'Table', 'demog');
      await this.dlgInputAction(dlg!, 'Choose [AGE, WEIGHT, HEIGHT] for "Impute"', 'Impute', 'AGE,HEIGHT,WEIGHT');
      await this.dlgInputAction(dlg!, 'Set "Data" to DIS_POP', 'Data', 'DIS_POP');
      await this.dlgInputAction(dlg!, 'Set the number of nearest neighbors to "5"', 'Nearest Neighbours', '5');

      await this.action('Click "OK" and wait for the values to be calculated.',
        grok.functions.onAfterRunAction.pipe(
          filter((call: DG.FuncCall) => call.func.name === 'MissingValuesImputation')
        ));

      await viewInputAction(`Set "Method" to "${method}"`, 'Method', method);
    };

    this.title('Train a model');
    await trainModel('Distributed Random Forest');

    const browseModels = async () => {
      await this.action('Click on "ML | Browse Models" to open the Models Browser',
        grok.events.onViewAdded.pipe(filter((v: DG.View) => v.type === DG.View.MODELS)));
    };

    this.title('Model performance');
    await browseModels();

  }
}
