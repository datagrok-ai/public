import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import { fromEvent, interval } from 'rxjs';
import { filter, map } from 'rxjs/operators';
import { Tutorial } from '../../../tutorial';


export class PredictiveModelingTutorial extends Tutorial {
  get name() {
    return 'Predictive Modeling';
  }
  get description() {
    return 'Predictive modeling uses statistics to predict outcomes';
  }
  get steps(){ return 12} //set thew number of steps
  private async buttonClickAction(root: HTMLElement, instructions: string, description: string, caption: string) {
    const btn = $(root).find('button.ui-btn').filter((idx, btn) => btn.textContent === caption)[0];
    if (btn == null) return;
    const source = fromEvent(btn, 'click');
    await this.action(instructions,description, source, btn);
  };

  private async columnInpAction(root: HTMLElement, instructions: string, description: string, caption: string, value: string) {
    const columnInput = $(root)
      .find('div.ui-input-root.ui-input-column')
      .filter((idx, inp) => $(inp).find('label.ui-label.ui-input-label')[0]?.textContent === caption)[0];
    if (columnInput == null) return;
    const source = interval(1000).pipe(
      map((_) => $(columnInput).find('div.d4-column-selector-column')[0]?.textContent),
      filter((val) => val === value));
    await this.action(instructions, description, source, columnInput);
  };

  private async columnsInpAction(root: HTMLElement, instructions: string, description: string, caption: string, value: string) {
    const columnsInput = $(root)
      .find('div.ui-input-root.ui-input-columns')
      .filter((idx, inp) => $(inp).find('label.ui-label.ui-input-label')[0]?.textContent === caption)[0];
    if (columnsInput == null) return;
    const source = interval(1000).pipe(
      map((_) => $(columnsInput).find('div.ui-input-editor > div.ui-input-column-names')[0]?.textContent),
      filter((val) => val === value));
    await this.action(instructions, description, source, columnsInput);
  };

  protected async _run() {
    /** Train model actions */

    this.describe('Predictive modeling uses statistics to predict outcomes.');

    const trainModel = async (method: string): Promise<void> => {
      const pmv = await this.openViewByType(
        'Click on "ML | Train Model..." to open a dialog for training models',
        'PredictiveModel', ''
      );

      // UI generation delay
      await new Promise((resolve) => setTimeout(resolve, 1500));

      const viewInputAction = async (instructions: string, description:string, caption: string, value: string) => {
        let inputRoot: HTMLDivElement;
        let select: HTMLSelectElement;
        $(pmv.root).find('.ui-input-root .ui-input-label span').each((idx, el) => {
          if (el.innerText == caption) {
            inputRoot = el.parentElement?.parentElement as HTMLDivElement;
            select = $(inputRoot).find('select')[0] as HTMLSelectElement;
          }
        });
        const source = fromEvent(select!, 'change');
        await this.action(instructions, description,
          source.pipe(map((event: any) => select.value), filter((v: string) => v === value)),
          inputRoot!);
      };

      //await viewInputAction('Set "Table" to "demog"', 'Table', 'demog');
      await this.columnInpAction(pmv.root, 'Set "Predict" to "SEX"','', 'Predict', 'SEX');
      await this.columnsInpAction(pmv.root, 'Set "Features" to ["AGE", "HEIGHT", "WEIGHT"]',
      'You will see that the column names are highlighted in red. Click on the question mark to the right of the input to see a warning message. It says that some of the columns contain missing values, and suggests two ways to address this.',
        'Features', '(3) AGE, HEIGHT, WEIGHT');

      const dlg = await this.openDialog('Let\'s use missing values imputation.','', 'Missing Values Imputation');

      await this.dlgInputAction(dlg, 'Set "Table" to "demog"', '', 'Table', 'demog');
      await this.dlgInputAction(dlg, 'Choose [AGE, HEIGHT, WEIGHT] for "Impute"', '', 'Impute', 'AGE,HEIGHT,WEIGHT');
      await this.dlgInputAction(dlg, 'Set "Data" to DIS_POP', '', 'Data', 'DIS_POP');
      await this.dlgInputAction(dlg, 'Set the number of nearest neighbors to "5"', '', 'Nearest Neighbours', '5');

      await this.action('Click "OK" and wait for the values to be calculated.', '',
        grok.functions.onAfterRunAction.pipe(
          filter((call: DG.FuncCall) => call.func.name === 'MissingValuesImputation'),
        ));

      await viewInputAction(`Set "Method" to "${method}"`, '', 'Method', method);
      await this.buttonClickAction(pmv.root, 'Click the "Train" button', '', 'TRAIN');
    };

    //this.title('Train a model');
    await trainModel('Distributed Random Forest');

    const browseModels = async () => {
      await this.openViewByType('Click on "Functions | Models" to open the Models Browser', '',  DG.View.MODELS);
    };

    //this.title('Model performance');
    await browseModels();

  }
}
