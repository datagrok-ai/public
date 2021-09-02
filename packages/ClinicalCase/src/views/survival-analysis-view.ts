import * as DG from "datagrok-api/dg";
import * as grok from 'datagrok-api/grok';
import * as ui from "datagrok-api/ui";
import { dataframeContentToRow } from "../data-preparation/utils";

export class SurvivalAnalysisView extends DG.ViewBase {

  survivalPlotDiv = ui.div();
  covariatesPlotDiv = ui.div();
  survivalColumns: string[];
  confIntervals = [ 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99 ];
  confInterval = 0.7;
  strata = '';
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
      'showContextMenu':false,
    });

    let applyFilters = ui.button('Apply to curves', () => { });
    applyFilters.addEventListener('click', (event) => {
      let test = this.survivalDataframe;
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

      this.root.appendChild(
        ui.splitV([
          ui.tabControl({
            'Survival data':
              ui.splitV([
                ui.divH([
                  this.survivalDataframe.plot.grid().root,
                  applyFilters
                ]),
                filters.root
              ], { style: { width: '100%' } }),
            'Survival chart':
              ui.divH([
                ui.divV([
                  confIntChoices.root,
                  strataChoices.root,
                ]),
                this.survivalPlotDiv ]),
            'Co-variates':
              ui.divH([
                ui.divV([
                  covariatesChoices.root
                ]),
                this.covariatesPlotDiv ])
          }).root
        ], { style: { width: '100%', height: '100%' } })
      );
      this.updatePlotDiv(survivalResult['plot'], this.survivalPlotDiv);
      console.warn(dataframeContentToRow(survivalResult['diagnostics']));
    });

  }

  private updatePlotDiv(img: string, div: HTMLDivElement) {
    div.innerHTML = '';
    div.append(ui.image(`data:image/png;base64,${img}`, 700, 700));
  }

  private updateSurvivalPlot() {
    grok.functions.call(
      "Clinicalcase:survivalPlot", {
      "survivalDf": this.survivalDataframe,
      "inputStrata": this.strata,
      "confInt": this.confInterval.toString()
    }).then((result) => {
      this.updatePlotDiv(result['plot'], this.survivalPlotDiv);
      console.warn(dataframeContentToRow(result['diagnostics']));
    });
  }

  private updateCovariatesPlot() {
    grok.functions.call(
      "Clinicalcase:covariates", {
      "covariatesDf": this.survivalDataframe,
      "coVariates": this.selectedCovariates.join(' + ')
    }).then((result) => {
      this.updatePlotDiv(result['plot'], this.covariatesPlotDiv);
      console.warn(dataframeContentToRow(result['diagnostics']));
    });
  }
}