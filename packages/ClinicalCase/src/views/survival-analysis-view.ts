import * as DG from "datagrok-api/dg";
import * as grok from 'datagrok-api/grok';
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";
import { createSurvivalData } from "../data-preparation/data-preparation";
import { dataframeContentToRow } from "../data-preparation/utils";

export class SurvivalAnalysisView extends DG.ViewBase {

  survivalPlotDiv = ui.box();
  covariatesPlotDiv = ui.box();
  survivalColumns: string[];
  confIntervals = [ 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99 ];
  covariatesOptions = [ 'AGE', 'SEX', 'RACE', 'ACTARM' ];
  endpointOptions = { 'DEATH': 'DTHDTC', 'TIME TO FIRST AE': 'DTHDTC' };
  confInterval = 0.7;
  strata = '';
  endpoint = '';
  covariates = [];
  survivalDataframe: DG.DataFrame;
  plotCovariates: string[];

  constructor() {
    super();

    this.survivalDataframe = grok.shell.table('survival');
    this.survivalColumns = this.survivalDataframe.columns.names();
    this.endpoint = Object.keys(this.endpointOptions)[ 0 ];
    let survivalOpions = [ '' ].concat(this.survivalColumns);
    let plotCovariatesOptions = this.survivalColumns.filter(it => it !== 'time' && it !== 'status');

    let endpointChoices = ui.choiceInput('Endpoint', Object.keys(this.endpointOptions)[ 0 ], Object.keys(this.endpointOptions));
    endpointChoices.onChanged((v) => {
      this.endpoint = this.endpointOptions[endpointChoices.value];
    });

    let covariatesChoices = ui.multiChoiceInput('Covariates', null, this.covariatesOptions);
    covariatesChoices.onChanged((v) => {
      this.covariates = covariatesChoices.value;
    });

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

    let plotCovariatesChoices = ui.multiChoiceInput('Covariates', null, plotCovariatesOptions);
    plotCovariatesChoices.onChanged((v) => {
      this.plotCovariates = plotCovariatesChoices.value;
      this.updateCovariatesPlot();
    });

    let filters = DG.Viewer.fromType('Filters', this.survivalDataframe, {
      'columnNames': this.survivalColumns,
      'showContextMenu': false,
    });

    let applyFilters = ui.bigButton('Apply to curves', () => { });
    applyFilters.addEventListener('click', (event) => {
      this.updateSurvivalPlot();
      if (this.plotCovariates) {
        this.updateCovariatesPlot();
      }
    });

    let createSurvivalDataframe = ui.bigButton('Create dataframe', () => { });
    createSurvivalDataframe.addEventListener('click', (event) => {
     this.survivalDataframe = createSurvivalData(study.domains.dm.clone(), this.endpointOptions[this.endpoint], this.covariates);
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
          'Create dataset':
            ui.splitV([
              ui.box(ui.panel([
                ui.divH([ endpointChoices.root,
                covariatesChoices.root ])
              ]), { style: { maxHeight: '150px' } }),
              ui.box(ui.div([ createSurvivalDataframe ]), { style: { maxHeight: '40px' } }),
              this.survivalDataframe ? this.survivalDataframe.plot.grid().root : null, ]),
          'Survival data':
            ui.splitV([
              filters.root,
              ui.box(ui.div([ applyFilters ]), { style: { maxHeight: '40px' } }),
              this.survivalDataframe ? this.survivalDataframe.plot.grid().root : null,
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
                plotCovariatesChoices.root
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
      "coVariates": this.plotCovariates.join(' + ')
    }).then((result) => {
      this.updatePlotDiv(result[ 'plot' ], this.covariatesPlotDiv);
      console.warn(dataframeContentToRow(result[ 'diagnostics' ]));
    });
  }
}