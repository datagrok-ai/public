import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { Observable } from 'rxjs';
import { Tutorial } from '@datagrok-libraries/tutorials/src/tutorial';


export class MultivariateAnalysisTutorial extends Tutorial {
  get name() { return 'Multivariate Analysis'; }
  get description() {
    return `Multivariate analysis models a response variable from many predictors at once, including predictors
    that correlate with each other. Learn to run partial least squares (PLS) regression and interpret its results.`;
  }
  get steps() { return 7; }

  get icon() {
    return '📊🔀';
  }

  demoTable: string = 'cars.csv';
  helpUrl: string = 'https://datagrok.ai/help/explore/multivariate-analysis';

  protected async _run() {
    this.header.textContent = this.name;

    grok.shell.windows.context.visible = false;
    grok.shell.windows.help.visible = false;

    this.describe(`Partial least squares (PLS) regression models a response variable from many predictors at once.
    It builds a linear model out of <b>latent factors</b> - combinations of the predictors that maximize the
    covariance with the response. This is why PLS works where ordinary multiple regression breaks down: numerous
    predictors, correlated predictors, or few observations.`);

    this.describe(ui.link('More about ' + this.name, this.helpUrl).outerHTML);

    const plsDlg = await this.openDialog('Click on "ML | Analyze | Multivariate Analysis..."',
      'Multivariate Analysis (PLS)', this.getMenuItem('ML', true));

    plsDlg.root.hidden = true;

    // We create fake dialog that runs analysis (since inputs of the "main" dialog are added using ui.form).
    const dlg = ui.dialog({title: 'Multivariate Analysis (PLS)', helpUrl: this.helpUrl});

    dlg.add(ui.input.column('Predict', {table: this.t!,
      filter: (col: DG.Column) => (col.type === DG.COLUMN_TYPE.INT) || (col.type === DG.COLUMN_TYPE.FLOAT)
    }));

    dlg.add(ui.input.columns('Using', {table: this.t!,
      value: [], available: this.t!.columns.toList().filter((col) =>
        (col.type === DG.COLUMN_TYPE.INT) || (col.type === DG.COLUMN_TYPE.FLOAT)
      ).map((col) => col.name),
    }));

    dlg.add(ui.input.int('Components'));
    dlg.add(ui.input.bool('Quadratic', {value: false}));

    dlg.add(ui.input.column('Names', {table: this.t!, filter: (col: DG.Column) => (col.type === DG.COLUMN_TYPE.STRING)}));

    let viewers = [] as DG.Viewer[];

    dlg.addButton('RUN', () => {
      dlg.close();
      plsDlg.getButton('RUN').click();

      setTimeout(() => {
        viewers = [...(grok.shell.v as DG.TableView).viewers];
        //console.log(viewers);
      }, 2000);
    }, undefined, 'Perform multivariate analysis');

    dlg.show();

    await this.dlgInputAction(dlg, 'Set "Predict" to "price"', 'Predict', 'price',
      'The response variable - the column the model predicts.');

    await this.dlgInputAction(dlg, 'Select all columns, except "price", as "Using"', 'Using',
      this.t!.columns.names().filter((n: string) => n !== 'model' && n !== 'price').join(','),
      `The predictors - the columns the model learns from. Click "All" in the column selection dialog, then uncheck "price".`);

    await this.dlgInputAction(dlg, 'Set the number of components to "3"', 'Components', '3',
      'The number of latent factors. Too few underfit the data, too many fit the noise.');

    await this.dlgInputAction(dlg, 'Set "Names" to "model"', 'Names', 'model',
      'The column that labels the points on the plots. Here, each row is a car model.');

    await this.action('Click "RUN" and wait for the analysis to complete', dlg.onClose);

    let viewerRoots: HTMLElement[];

    const viewerMd = [
      '# Observed vs. Predicted\n\n' +
      'The actual price against the predicted one. The closer the points lie to the diagonal, ' +
      'the better the prediction.',
      '# Scores\n\n' +
      'Similar cars sit close together, dissimilar ones far apart:\n\n' +
      '* `Volvos` are close to each other\n' +
      '* `Porsche` and `Mercedes` are far apart',
      '# Loadings\n\n' +
      'How strongly each latent factor describes each feature. Features sitting together carry the same information, ' +
      'features on opposite sides are inversely related:\n\n' +
      '* `width`, `length` and `curb.weight` form one group\n' +
      '* `city.mpg` and `highway.mpg` are opposite to them - bigger cars burn more fuel',
      '# Variable Importance\n\n' +
      'How much each feature helps the model explain the price, measured on standardized data - so features in ' +
      'different units are comparable:\n\n' +
      '* above 1 means more than an average predictor: `eng.size` leads with 1.52\n' +
      '* below 0.8 means weak: `peak.rpm` and `two.doors` are candidates for removal',
      '# Explained Variance\n\n' +
      'The share of the price variation the latent factors capture, added up over the components:\n\n' +
      '* closer to one means a better fit\n' +
      '* one component already covers 70%, three cover 92%',
    ];

    let idx = 0;
    let hint: HTMLElement;
    let msg: HTMLDivElement;
    let popup: HTMLDivElement;
    const nextBtn = ui.button('next', () => hint.click(), 'Go to the next viewer');
    const prevBtn = ui.button('prev', () => {
      idx -= 2;
      hint.click();
    }, 'Go to the previous viewer');
    const doneBtn = ui.button('done', () => hint.click(), 'Complete this tutorial');
    const btnsDiv = ui.divH([prevBtn, nextBtn, doneBtn]);
    btnsDiv.style.marginLeft = 'auto';
    btnsDiv.style.marginRight = '0';

    const step = () => {
      if (idx < viewerRoots.length) {
        msg = ui.divV([ui.markdown(viewerMd[idx]), btnsDiv]);
        popup = ui.hints.addHint(viewerRoots[idx], msg, 'left');
        doneBtn.hidden = (idx < viewerRoots.length - 1);
        nextBtn.hidden = (idx === viewerRoots.length - 1);
        prevBtn.hidden = (idx < 1);
        hint = ui.hints.addHintIndicator(popup, undefined, 3000);
        hint.onclick = () => {
          popup.remove();
          ++idx;
          step();
        };
      }
    };

    setTimeout(async () => {
      viewerRoots = [
        viewers[1].root, // Observed vs. Predicted scatterplot
        viewers[4].root, // Scores scatterplot
        viewers[3].root, // Loadings scatterplot
        viewers[6].root, // Variable importance bar chart
        viewers[5].root, // Explained variances bar chart
      ];

      step();
    }, 2100);

    await this.action('Explore each viewer', new Observable((subscriber: any) => {
      //@ts-ignore
      $(doneBtn).one('click', () => subscriber.next(true));
    }), undefined, 'Click "Next" to switch to the next viewer');
  }
}
