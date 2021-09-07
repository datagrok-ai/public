import * as DG from "datagrok-api/dg";
import * as grok from 'datagrok-api/grok';
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";
import { createSurvivalData } from "../data-preparation/data-preparation";
import { dataframeContentToRow } from "../data-preparation/utils";

export class SurvivalAnalysisView extends DG.ViewBase {

  survivalPlotDiv = ui.box();
  covariatesPlotDiv = ui.box();
  survivalGridDivCreate = ui.box();
  survivalGridDivFilter = ui.box();
  survivalFilterDiv = ui.box();
  strataChoicesDiv = ui.div();
  plotCovariatesChoicesDiv = ui.div();
  strataChoices: DG.InputBase;
  plotCovariatesChoices: DG.InputBase;
  survivalColumns = [];
  confIntervals = [ 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99 ];
  survivalOptions = [''];
  covariatesOptions = [ 'AGE', 'SEX', 'RACE', 'ACTARM' ];
  endpointOptions = { 'TIME TO FIRST SAE': 'AESTDTC', 'DEATH': 'DTHDTC' };
  confInterval = 0.7;
  strata = '';
  endpoint = 'TIME TO FIRST SAE';
  covariates = [];
  survivalDataframe: DG.DataFrame;
  plotCovariates: string[];

  constructor() {
    super();

    this.endpoint = Object.keys(this.endpointOptions)[ 0 ];

    let endpointChoices = ui.choiceInput('Endpoint', Object.keys(this.endpointOptions)[ 0 ], Object.keys(this.endpointOptions));
    endpointChoices.onChanged((v) => {
      this.endpoint = endpointChoices.value;
    });

    let covariatesChoices = ui.multiChoiceInput('Covariates', null, this.covariatesOptions);
    covariatesChoices.onChanged((v) => {
      this.covariates = covariatesChoices.value;
    });

    let confIntChoices = ui.choiceInput('Confidence Interval', this.confIntervals[ 0 ], this.confIntervals);
    confIntChoices.onChanged((v) => {
      this.confInterval = confIntChoices.value;
      this.updateSurvivalPlot();
    });

    this.updateStrataChoices();
    this.updatePlotCovariatesChoices();

    let applyFilters = ui.bigButton('Apply to curves', () => { });
    applyFilters.addEventListener('click', (event) => {
      this.updateSurvivalPlot();
      if (this.plotCovariates) {
        this.updateCovariatesPlot();
      }
    });

    let createSurvivalDataframe = ui.bigButton('Create dataframe', () => { });
    createSurvivalDataframe.addEventListener('click', (event) => {
     this.survivalDataframe = createSurvivalData(this.endpointOptions[this.endpoint], this.covariates);
     this.survivalColumns = this.survivalDataframe.columns.names();
     this.survivalOptions = [''].concat(this.survivalColumns.filter(it => it !== 'time' && it !== 'status' && it !== 'USUBJID'));
     this.updateStrataChoices();
     this.updatePlotCovariatesChoices();
     this.updateDivInnerHTML(this.strataChoicesDiv, this.strataChoices.root);
     this.updateDivInnerHTML(this.plotCovariatesChoicesDiv, this.plotCovariatesChoices.root);
     this.updateDivInnerHTML(this.survivalGridDivCreate, this.survivalDataframe.plot.grid().root);
     this.updateDivInnerHTML(this.survivalGridDivFilter, this.survivalDataframe.plot.grid().root);
     this.updateDivInnerHTML(this.survivalFilterDiv, this.getFilters());
     this.updateSurvivalPlot();
    });


      this.root.className = 'grok-view ui-box';
      this.root.append(
        ui.tabControl({
          'Create dataset':
            ui.splitH([
              ui.box(ui.panel([
                ui.divV([ endpointChoices.root,
                covariatesChoices.root,
                ui.box(ui.div([ createSurvivalDataframe ]), { style: { maxHeight: '40px' } })])
              ]), { style: { maxWidth: '250px' } }),
              this.survivalGridDivCreate ]),
          'Survival data':
            ui.splitV([
              this.survivalFilterDiv,
              ui.box(ui.div([ applyFilters ]), { style: { maxHeight: '40px' } }),
              this.survivalGridDivFilter
            ]),
          'Survival chart':
            ui.splitV([
              ui.box(ui.panel([
                ui.divH([ confIntChoices.root,
                this.strataChoicesDiv ])
              ]), { style: { maxHeight: '80px' } }),
              this.survivalPlotDiv ]),
          'Co-variates':
            ui.splitH([
              ui.panel([
                this.plotCovariatesChoicesDiv
              ], { style: { maxWidth: '180px' } }),
              this.covariatesPlotDiv ])
        }).root
      );

  }

  private updateDivInnerHTML(div: HTMLDivElement, content: any){
    div.innerHTML = '';
    div.append(content);
  }

  private updateSurvivalPlot() {
    grok.functions.call(
      "Clinicalcase:survivalPlot", {
      "survivalDf": this.survivalDataframe,
      "inputStrata": this.strata,
      "confInt": this.confInterval.toString()
    }).then((result) => {
      //@ts-ignore
      this.updateDivInnerHTML(this.survivalPlotDiv, ui.image(`data:image/png;base64,${result[ 'plot' ]}`));
      console.warn(dataframeContentToRow(result[ 'diagnostics' ]));
    });
  }

  private updateCovariatesPlot() {
    grok.functions.call(
      "Clinicalcase:covariates", {
      "covariatesDf": this.survivalDataframe,
      "coVariates": this.plotCovariates.join(' + ')
    }).then((result) => {
      //@ts-ignore
      this.updateDivInnerHTML(this.covariatesPlotDiv, ui.image(`data:image/png;base64,${result[ 'plot' ]}`));
      console.warn(dataframeContentToRow(result[ 'diagnostics' ]));
    });
  }

  private getFilters(){
    return DG.Viewer.fromType('Filters', this.survivalDataframe, {
      'columnNames': this.survivalColumns.filter(it => it !== 'time' && it !== 'status' && it !== 'USUBJID'),
      'showContextMenu': false,
    }).root
  }

  private updateStrataChoices(){
    this.strataChoices = ui.choiceInput('Strata', this.survivalOptions[ 0 ], this.survivalOptions);
    this.strataChoices.onChanged((v) => {
      this.strata = this.strataChoices.value;
      this.updateSurvivalPlot();
    });
  }

  private updatePlotCovariatesChoices(){
    this.plotCovariatesChoices = ui.multiChoiceInput('Covariates', null, this.survivalColumns.filter(it => it !== 'time' && it !== 'status' && it !== 'USUBJID'));
    this.plotCovariatesChoices.onChanged((v) => {
      this.plotCovariates = this.plotCovariatesChoices.value;
      this.updateCovariatesPlot();
    });
  }
}