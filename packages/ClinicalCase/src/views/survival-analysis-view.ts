import * as DG from "datagrok-api/dg";
import * as grok from 'datagrok-api/grok';
import * as ui from "datagrok-api/ui";
import { dataframeContentToRow } from "../data-preparation/utils";

export class SurvivalAnalysisView extends DG.ViewBase {

  survivalPlotDiv = ui.box();
  covariatesPlotDiv = ui.box();
  survivalColumns: string[];
  confIntervals = [ 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99 ];
  confInterval = 0.7;
  strata = '';
  endpoint = '';
  survivalDataframe: DG.DataFrame;
  selectedCovariates: string[];

  constructor() {
    super();

    this.survivalDataframe = grok.shell.table('survival');
    this.survivalColumns = this.survivalDataframe.columns.names();
    let survivalOpions = [ '' ].concat(this.survivalColumns);
    let covariatesOptions = this.survivalColumns.filter(it => it !== 'time' && it !== 'status');


    let confIntChoices = ui.choiceInput('Confidence Intreval', this.confIntervals[ 0 ], this.confIntervals);
    confIntChoices.onChanged((v) => {
      this.confInterval = confIntChoices.value;
      this.updateSurvivalPlot();
    });

    let strataChoices = ui.choiceInput('Strata', survivalOpions[ 0 ], survivalOpions);
    strataChoices.onChanged((v) => {
      this.strata = strataChoices.value;
      this.updateSurvivalPlot();
    });

    let covariatesChoices = ui.multiChoiceInput('Covariates', null, covariatesOptions);
    covariatesChoices.onChanged((v) => {
      this.selectedCovariates = covariatesChoices.value;
      this.updateCovariatesPlot();
    });

    let filters = DG.Viewer.fromType('Filters', this.survivalDataframe, {
      'columnNames': this.survivalColumns,
      'showContextMenu': false,
    });

    let applyFilters = ui.bigButton('Apply to curves', () => { });
    applyFilters.addEventListener('click', (event) => {
      this.updateSurvivalPlot();
      if (this.selectedCovariates) {
        this.updateCovariatesPlot();
      }
    });

    grok.functions.call(
      "Clinicalcase:survivalPlot", {
      "survivalDf": this.survivalDataframe,
      "inputStrata": this.strata,
      "confInt": this.confInterval.toString()

    }).then((survivalResult) => {

      this.root.className = 'grok-view ui-box';
      this.root.append(
        ui.tabControl({
          'Survival data':
            ui.splitV([
              filters.root,
              ui.box(ui.div([ applyFilters ]), { style: { maxHeight: '40px' } }),
              this.survivalDataframe.plot.grid().root,
            ]),
          'Survival chart':
            ui.splitV([
              ui.box(ui.panel([
                ui.divH([ confIntChoices.root,
                strataChoices.root ])
              ]), { style: { maxHeight: '80px' } }),
              this.survivalPlotDiv ]),
          'Co-variates':
            ui.splitH([
              ui.panel([
                covariatesChoices.root
              ], { style: { maxWidth: '150px' } }),
              this.covariatesPlotDiv ])
        }).root
      );
      this.updatePlotDiv(survivalResult[ 'plot' ], this.survivalPlotDiv);
      console.warn(dataframeContentToRow(survivalResult[ 'diagnostics' ]));
    });

  }

  private updatePlotDiv(img: string, div: HTMLDivElement) {
    div.innerHTML = '';
    //@ts-ignore
    div.append(ui.image(`data:image/png;base64,${img}`));
  }

  private updateSurvivalPlot() {
    grok.functions.call(
      "Clinicalcase:survivalPlot", {
      "survivalDf": this.survivalDataframe,
      "inputStrata": this.strata,
      "confInt": this.confInterval.toString()
    }).then((result) => {
      this.updatePlotDiv(result[ 'plot' ], this.survivalPlotDiv);
      console.warn(dataframeContentToRow(result[ 'diagnostics' ]));
    });
  }

  private updateCovariatesPlot() {
    grok.functions.call(
      "Clinicalcase:covariates", {
      "covariatesDf": this.survivalDataframe,
      "coVariates": this.selectedCovariates.join(' + ')
    }).then((result) => {
      this.updatePlotDiv(result[ 'plot' ], this.covariatesPlotDiv);
      console.warn(dataframeContentToRow(result[ 'diagnostics' ]));
    });
  }
}