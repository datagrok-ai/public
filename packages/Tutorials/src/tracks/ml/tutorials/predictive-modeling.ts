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
        grok.events.onViewAdded.pipe(filter((v) => {
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
          source.pipe(map((event) => select.value), filter((v) => v === value)),
          inputRoot!);
      }

      //await viewInputAction('Set "Table" to "demog"', 'Table', 'demog');
      await viewInputAction(`Set "Method" to "${method}"`, 'Method', method);
    };

    this.title('Train a model');
    await trainModel('Distributed Random Forest');

    const browseModels = async () => {
      await this.action('Click on "ML | Browse Models" to open the Models Browser',
        grok.events.onViewAdded.pipe(filter((v) => v.type === DG.View.MODELS)));
    };

    this.title('Model performance');
    await browseModels();

  }
}
