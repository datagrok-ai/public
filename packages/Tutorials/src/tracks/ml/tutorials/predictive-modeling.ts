import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import $ from 'cash-dom';
import { filter, map } from 'rxjs/operators';
import { Tutorial, TutorialPrerequisites } from '@datagrok-libraries/tutorials/src/tutorial';
import {delay} from '@datagrok-libraries/test/src/test';
import { interval, Observable } from 'rxjs';


export class PredictiveModelingTutorial extends Tutorial {
  get name() {
    return 'Predictive Modeling';
  }
  get description() {
    return 'Predictive modeling is a statistical technique used to predict outcomes based on historical data.';
  }
  get steps() { return 26; }

  get icon() {
    return '📈🔮';
  }

  helpUrl: string = 'https://datagrok.ai/help/learn';
  prerequisites: TutorialPrerequisites = {jupyter: true, h2o: true};

  protected async _run() {
    this.header.textContent = this.name;
    this.describe('Predictive modeling is a statistical technique used to predict outcomes ' +
      'based on historical data. In the next steps, we will train a few models, look at their ' +
      'performance, learn how to apply a model to a dataset and share it with others, and ' +
      'lastly, compare the models we have trained.');

    this.describe(ui.link('More about ' + this.name, this.helpUrl).outerHTML); 

    /** Train model actions */
    const trainModel = async (method: string, skipPMVOpening = false): Promise<void> => {
      const pmv = await this.openViewByType(
        skipPMVOpening ? 'Return to the "Predictive model" tab' :
        'Click on "ML | Models | Train Model..." to open a dialog for training models',
        'PredictiveModel',
        skipPMVOpening ? [] : this.getMenuItem('ML'),
      );

      // UI generation delay
      if (!skipPMVOpening) {
        await new Promise((resolve) => setTimeout(resolve, 2000));
      }

      const hasNulls = this.t!.columns
        .byNames(['AGE', 'HEIGHT', 'WEIGHT'])
        .some((col) => col?.stats.missingValueCount > 0);

      await this.choiceInputAction(pmv.root, 'Set "Table" to "demog"', 'Table', 'demog');
      await this.columnInpAction(pmv.root, 'Set "Predict" to "SEX"', 'Predict', 'SEX');
      const missingValuesWarning = 'You will see that the column names are highlighted in red. ' +
        'Click on the question mark to the right of the input to see a warning message. It says ' +
        'that some of the columns contain missing values, and suggests two ways to address this. ' +
        'Let\'s use missing values imputation.';
      await this.columnsInpAction(pmv.root, 'Set "Features" to ["AGE", "HEIGHT", "WEIGHT"]',
        'Features', '(3) AGE, HEIGHT, WEIGHT');

      if (hasNulls) {
        const dlg = await this.openDialog('Click on "Missing values imputation"', 'Missing Values Imputation',
          $(pmv.root).find('label.d4-link-label.d4-validator-helper-link')[0], missingValuesWarning);

        await this.dlgInputAction(dlg, 'Set "Table" to "demog"', 'Table', 'demog');
        await this.dlgInputAction(dlg, 'Choose ["AGE", "HEIGHT", "WEIGHT"] for "Impute"', 'impute', 'AGE,HEIGHT,WEIGHT');
        await this.dlgInputAction(dlg, 'Set "Data" to "DIS_POP"', 'data', 'DIS_POP');
        await this.dlgInputAction(dlg, 'Set the number of nearest neighbors to "5"', 'Nearest Neighbours', '5');
  
        await this.action('Click "OK" and wait for the values to be calculated.',
          grok.functions.onAfterRunAction.pipe(
            filter((call: DG.FuncCall) => call.func.name === 'MissingValuesImputation'),
          ));
      }

      await this.choiceInputAction(pmv.root, `Set "Method" to "${method}"`, 'Method', method);
      await this.buttonClickAction(pmv.root, 'Click the "Train" button', 'TRAIN');

      await delay(3000); // model training delay
    };

    this.title('Train a model');
    const method = 'Distributed Random Forest';
    await trainModel(method);

    this.title('Model performance, application, and sharing');

    const funcPaneHints = grok.shell.sidebar.getPane('Browse').header;
    let v = grok.shell.view('Browse');
    if (v) v.close();
    await this.action('Open Browse view', new Observable((subscriber: any) => {
      $(funcPaneHints).one('click', () => subscriber.next(true));
    }), funcPaneHints);
    const pmBrowserDescription = 'This is Predictive Models Browser. Here, you can browse ' +
      'models that you trained or that were shared with you. In the next steps, we will look ' +
      'at model performance, apply a model to a dataset, and share the model.';
    let currentUrl = location.href;
    await this.action('In Browse, click on Platform | Models to open the Models Browser',
    new Observable((subscriber: any) => {
      const int = setInterval(() => {
        if (location.href !== currentUrl) {
          currentUrl = location.href;
          if (grok.shell.v.path.startsWith('/models')) {
            subscriber.next(true);
            clearInterval(int);
          }
        }
      }, 100);
    }), [], pmBrowserDescription);

    let modelPP: DG.Accordion;
    const modelName = `Predict SEX by AGE, HEIGHT, WEIGHT using ${method}`;
    const modelsViewInfo = 'Search for the model you trained (applicable models are always ' +
      'at the top of the list, also, you can search for it by the previously defined name). ' +
      'Now, select the model by clicking on it.';
    await this.action('Find your model', grok.events.onAccordionConstructed.pipe(
      map((acc) => {
        modelPP = acc;
        return acc.context;
      }), filter((model) => model instanceof DG.Model && model.friendlyName === modelName)),
      null, modelsViewInfo);

    await this.action('Expand the "Performance" pane in the context panel',
      interval(1000).pipe(filter(() => modelPP?.getPane('Performance')?.expanded)));

    await this.contextMenuAction('Right-click on the trained model and select "Apply to | ' +
      `${this.t!.toString()}"`, this.t!.toString(), null, 'The result will be available in ' +
      'the selected table as a column named "Outcome". Check it out in the opened table view.');

    await this.contextMenuAction('Right-click on the model and select "Share...".', 'Share...',
      null, 'The model you trained is only visible to you, but it is easy to share it with other ' +
      'users from this dialog: just add a few users (e.g., "All users") and click "OK".');

    this.title('Compare models');

    await trainModel('Gradient Boosting Machine', true);
    v = grok.shell.view('Browse');
    if (v) v.close();
    await this.action('Open Browse view', new Observable((subscriber: any) => {
      $(funcPaneHints).one('click', () => subscriber.next(true));
    }), funcPaneHints);
    await this.action('In Browse, click on Platform | Models to open the Models Browser',
    new Observable((subscriber: any) => {
      const int = setInterval(() => {
        if (location.href !== currentUrl) {
          currentUrl = location.href;
          if (grok.shell.v.path.startsWith('/models')) {
            subscriber.next(true);
            clearInterval(int);
          }
        }
      }, 100);
    }));

    await this.action('Find the trained models and compare them',
      grok.events.onViewAdded.pipe(filter((view) => view.name == 'Compare models')),
      $('div.d4-accordion-pane-header').filter((idx, el) => el.textContent == 'Commands')[0],
      'Search for the trained models and select them holding <b>Shift</b>. ' +
      'Then run the "Compare" command given at the context panel.');
  }
}
