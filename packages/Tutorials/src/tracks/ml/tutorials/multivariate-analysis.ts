import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { filter } from 'rxjs/operators';
import { Tutorial, TutorialPrerequisites } from '@datagrok-libraries/tutorials/src/tutorial';


export class MultivariateAnalysisTutorial extends Tutorial {
  get name() { return 'Multivariate Analysis'; }
  get description() {
    return `Multivariate analysis reveals complex interactions and patterns within a dataset. It uses 
    statistical techniques to explore the relationships among multiple variables.`;
  }
  get steps() { return 6; }
    
  demoTable: string = 'cars.csv';
  helpUrl: string = 'https://datagrok.ai/help/explore/multivariate-analysis/pls';
  prerequisites: TutorialPrerequisites = {jupyter: true};

  protected async _run() {
    this.header.textContent = this.name;

    this.describe(`Partial Least Squares Regression (PLS) models the relationship 
    between independent variables (predictors) and dependent variables (responses). Use it 
    when there are a large number of predictors, multicollinearity among them, and relatively few observations.`);

    this.describe(ui.link('More about ' + this.name, this.helpUrl).outerHTML);

    const plsDlg = await this.openDialog('Click on "ML | Analyze | Multivariate Analysis..."',
      'Multivariate Analysis (PLS)', this.getMenuItem('ML'));

    plsDlg.root.hidden = true;
    
    // We create fake dialog that runs analysis (since inputs of the "main" dialog are added using ui.form).
    const dlg = ui.dialog({title: 'Multivariate Analysis (PLS)', helpUrl: this.helpUrl});

    dlg.add(ui.columnInput('Predict', this.t!, null, () => {}, {
      filter: (col: DG.Column) => (col.type === DG.COLUMN_TYPE.INT) || (col.type === DG.COLUMN_TYPE.FLOAT)
    }));
  
    dlg.add(ui.columnsInput('Using', this.t!, () => {}, {
      available: this.t!.columns.toList().filter((col) => 
        (col.type === DG.COLUMN_TYPE.INT) || (col.type === DG.COLUMN_TYPE.FLOAT)
      ).map((col) => col.name),
      checked: [],
    }));
  
    dlg.add(ui.intInput('Components', null, () => {}));
  
    dlg.add(ui.columnInput('Names', this.t!, null, () => {}, {
      filter: (col: DG.Column) => (col.type === DG.COLUMN_TYPE.STRING)
    }));

    let viewers = [] as DG.Viewer[];
  
    dlg.addButton('RUN', () => {
      dlg.close();
      plsDlg.getButton('RUN').click();

      setTimeout(() => {
        viewers = [...(grok.shell.v as DG.TableView).viewers];
        console.log(viewers);
      }, 2000);
    }, undefined, 'Perform multivariate analysis');
  
    dlg.show();
      
    await this.dlgInputAction(dlg, 'Set "Predict" to "price"', 'Predict', 'price',
      'Specify column with the response variable.');
      
    await this.dlgInputAction(dlg, 'Select all columns, except "price", as "Using"', 'Using',
      this.t!.columns.names().filter((n: string) => n !== 'model' && n !== 'price').join(','),
      `Set columns with predictors' values. Click "All" in the column selection dialog and uncheck "price".`);

    await this.dlgInputAction(dlg, 'Set the number of components to "3"', 'Components', '3',
      'Define the number of the latent factors.');
    
    await this.dlgInputAction(dlg, 'Set "Names" to "model"', 'Names', 'model',
      'Select column with data samples names.');

    await this.action('Click "RUN" and wait for the analysis to complete.', dlg.onClose);

    setTimeout(() => {
      const wizard = ui.hints.addTextHint({title: 'Explore', pages: [
        {
          caption: ui.h2('Observed vs. Predicted'),
          text: 'Closer to the line means better price prediction.',
          showNextTo: viewers[1].root,
        },
        {
          caption: ui.h2('Scores'),
          text: 'Similarities & dissimilarities: alfaromeo and mercedes are different.',
          showNextTo: viewers[4].root,
        },
        {
          caption: ui.h2('Loadings'),
          text: 'The impact of each feature on the latent factors: higher loading means stronger influence.',
          showNextTo: viewers[3].root,
        },
        {
          caption: ui.h2('Regression Coefficients'),
          text: 'Parameters of the obtained linear model: the "diesel" feature affects the price the most.',
          showNextTo: viewers[2].root,
        },
        {
          caption: ui.h2('Explained Variance'),
          text: 'How well the latent components fit source data: closer to one means better fit.',
          showNextTo: viewers[5].root,
        },
      ]});

      wizard.helpUrl = 'https://datagrok.ai/help/explore/multivariate-analysis/pls';
    }, 2100);
  }
}
