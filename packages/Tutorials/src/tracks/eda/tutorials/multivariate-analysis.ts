import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { filter } from 'rxjs/operators';
import { Tutorial } from '../../../tutorial';


export class MultivariateAnalysisTutorial extends Tutorial {
  get name() { return 'Multivariate Analysis'; }
  get description() {
    return 'Multivariate analysis (MVA) is based on the statistical ' +
    'principle of multivariate statistics, which involves observation ' +
    'and analysis of more than one statistical outcome variable at a time';
  }

  demoTable: string = 'cars.csv';

  protected async _run() {
    this.header.textContent = 'Multivariate Analysis';
    let dlg: DG.Dialog;

    await this.action('Click on "ML | Multivariate Analysis (PLS)..."',
      grok.events.onDialogShown.pipe(filter((dialog: any) => {
        if (dialog.title === 'Multivariate Analysis (PLS)') {
          dlg = dialog;
          return true;
        }
        return false;
      })));

    this.describe('This is multivariate analysis dialog. It uses a statistical method that resembles ' +
      'principal components regression. Instead of finding hyperplanes of maximum variance ' +
      'between the response and independent variables, it finds a linear regression model by projecting ' +
      'the predicted variables and the observable variables to a new space.');

    this.describe('In the following example, we will predict a car price by its attributes, and see ' +
      'how different measures are related to each other and to the outcome.');

    await this.dlgInputAction(dlg!, 'Set "Table" to "cars"', 'Table', 'cars');
    await this.dlgInputAction(dlg!, 'Set "Predict" to "price"', 'Predict', 'price');

    this.describe('Set "Features" to "all". Click "All" in the column selection dialog.');
    // await this.dlgInputAction(dlg!, 'Set "Features" to "all". Click "All" in the column selection dialog',
    //   'Features', this.t.columns.names().join(','));
    await this.dlgInputAction(dlg!, 'Set the number of components to "3"', 'Components', '3');
    this.describe(`Once you run the analysis, the following visualizations will appear:
    <ul style="list-style-type:disc">
      <li><code>Scores</code> is a scatter plot that shows correlation between observations</li>
      <li><code>Explained variance</code> is a bar chart with variable variance explained</li>
      <li><code>Correlation loadings</code> is a scatter plot that shows correlation between variables</li>
      <li><code>Predicted vs. reference</code> is a scatter plot with a regression line comparing predicted vs reference outcomes</li>
      <li><code>Regression coefficients</code> is a bar chart with regression coefficients</li>
    </ul>`);

    await this.action('Click "OK" and wait for the analysis to complete.',
      grok.functions.onAfterRunAction.pipe(filter((call: DG.FuncCall) => call.func.name === 'MultivariateAnalysis')),
    );

    this.describe('Look into these findings to figure out which variables correlate with each other ' +
      'and which of them have a greater impact on the outcome.');
    this.describe('See the context help on the right to learn more about multivariate analysis, ' +
      'or navigate to a documentation page by clicking the "Open in new tab" button.');
  }
}
