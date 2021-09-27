import * as DG from "datagrok-api/dg";
import { InputBase } from "datagrok-api/dg";
import * as grok from 'datagrok-api/grok';
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";
import { SURVIVAL_ANALYSIS_GUIDE, TREATMENT_ARM } from "../constants";
import { createSurvivalData } from "../data-preparation/data-preparation";
import { dataframeContentToRow } from "../data-preparation/utils";
import { updateDivInnerHTML } from "./utils";

export class SurvivalAnalysisView extends DG.ViewBase {

  survivalPlotDiv = ui.box(ui.divText('Create dataset',{style:{color:'var(--grey-3)',marginTop:'30px', alignItems:'center'}}));
  covariatesPlotDiv = ui.box(ui.divText('Create dataset',{style:{color:'var(--grey-3)',marginTop:'30px', alignItems:'center'}}));
  survivalGridDivCreate = ui.box();
  survivalGridDivFilter = ui.box();
  survivalFilterDiv = ui.box(ui.divText('Create dataset',{style:{color:'var(--grey-3)',marginTop:'30px', alignItems:'center'}}));
  strataChoicesDiv = ui.div();
  plotCovariatesChoicesDiv = ui.div();
  strataChoices: DG.InputBase;
  plotCovariatesChoices: DG.InputBase;
  survivalColumns = [];
  confIntervals = [ 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99 ];
  survivalOptions = [''];
  covariatesOptions = [ 'AGE', 'SEX', 'RACE', TREATMENT_ARM ];
  endpointOptions = { 'SAE': 'AESTDTC', 'DEATH': 'DTHDTC', 'HOSPITALIZATION': 'AESTDTC', 'DRUG RELATED AE': 'AESTDTC' };
  confInterval = 0.7;
  strata = '';
  endpoint = 'TIME TO FIRST SAE';
  covariates = [];
  survivalDataframe: DG.DataFrame;
  plotCovariates: string[];
  filterChanged = false;

  constructor(name) {
    super(name);

    this.name = name;
    this.endpoint = Object.keys(this.endpointOptions)[ 0 ];

    let endpointChoices = ui.choiceInput('Endpoint', Object.keys(this.endpointOptions)[ 0 ], Object.keys(this.endpointOptions));
    endpointChoices.onChanged((v) => {
      this.endpoint = endpointChoices.value;
    });

    let covariatesChoices = ui.multiChoiceInput('Covariates', null, this.covariatesOptions);
    covariatesChoices.onChanged((v) => {
      this.covariates = covariatesChoices.value;
    });

    let confIntChoices = ui.choiceInput('Confidence', this.confIntervals[ 0 ], this.confIntervals);
    confIntChoices.onChanged((v) => {
      this.confInterval = confIntChoices.value;
      this.updateSurvivalPlot();
    });

    this.updateStrataChoices();
    this.updatePlotCovariatesChoices();

    let createSurvivalDataframe = ui.bigButton('Create dataset', () => { 
     this.refreshDataframe();
     this.filterChanged = true;
     this.survivalDataframe.onFilterChanged.subscribe((_) => {
       this.filterChanged = true;
      });
    });

    let tabControl = ui.tabControl(null, false);
    tabControl.addPane('Dataset', () => ui.splitV([this.survivalFilterDiv, this.survivalGridDivCreate]));
    tabControl.addPane('Survival Chart', () => {
      this.updateChartsAfterFiltering();
      this.survivalPlotDiv;
    });
    tabControl.addPane('Covariates', () => {
      this.updateChartsAfterFiltering();
      this.covariatesPlotDiv;
    });

    let customTitle = {style:{
      'color':'var(--grey-6)',
      'margin-top':'8px',
      'font-size':'16px',
    }};

      this.root.className = 'grok-view ui-box';
      this.setRibbonPanels([
        [
          ui.icons.info(()=>{
            updateDivInnerHTML(guide, ui.info(SURVIVAL_ANALYSIS_GUIDE,'Survival Analysis Quick Guide', false));
          })
        ]
      ])

      let guide = ui.info(SURVIVAL_ANALYSIS_GUIDE,'Survival Analysis Quick Guide', false);
      this.root.append(ui.splitV([
        guide,
        ui.splitH([
          ui.box(ui.div([
            ui.panel([
              ui.divText('Dataset', customTitle),
              ui.inputs([ 
                endpointChoices,
                covariatesChoices,
                //@ts-ignore
                ui.buttonsInput([createSurvivalDataframe])
              ])
            ]),
            ui.panel([
              ui.divText('Survival Parameters', customTitle),
              ui.inputs([ 
                confIntChoices,
                //@ts-ignore
                this.strataChoicesDiv,
              ])
            ]),
            ui.panel([
              ui.divText('Co-Variates', customTitle),
              ui.inputs([ 
                //@ts-ignore
                this.plotCovariatesChoicesDiv
              ])
            ])
          ]), { style: { maxWidth: '300px' }}),
          tabControl.root
        ])
      ]))
      //@ts-ignore
      guide.parentNode.style.flexGrow = '0';
      //@ts-ignore
      guide.parentNode.classList = 'ui-div';

  }

  private updateSurvivalPlot() {
    grok.functions.call(
      "Clinicalcase:survivalPlot", {
      "survivalDf": this.survivalDataframe,
      "inputStrata": this.strata,
      "confInt": this.confInterval.toString()
    }).then((result) => {
      //@ts-ignore
      updateDivInnerHTML(this.survivalPlotDiv, ui.image(`data:image/png;base64,${result[ 'plot' ]}`));
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
      updateDivInnerHTML(this.covariatesPlotDiv, ui.image(`data:image/png;base64,${result[ 'plot' ]}`));
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

  private refreshDataframe(){
    this.survivalDataframe = createSurvivalData(this.endpoint, this.endpointOptions[this.endpoint], this.covariates);
     this.survivalColumns = this.survivalDataframe.columns.names();
     this.survivalOptions = [''].concat(this.survivalColumns.filter(it => it !== 'time' && it !== 'status' && it !== 'USUBJID'));
     this.strata = '';
     this.updateStrataChoices();
     this.updatePlotCovariatesChoices();
     updateDivInnerHTML(this.strataChoicesDiv, this.strataChoices.root);
     updateDivInnerHTML(this.plotCovariatesChoicesDiv, this.plotCovariatesChoices.root);
     updateDivInnerHTML(this.survivalGridDivCreate, this.survivalDataframe.plot.grid().root);
     updateDivInnerHTML(this.survivalGridDivFilter, this.survivalDataframe.plot.grid().root);
     updateDivInnerHTML(this.survivalFilterDiv, this.getFilters());
     this.updateSurvivalPlot();
  }

  private updateChartsAfterFiltering(){
    if(this.filterChanged){
      this.updateSurvivalPlot();
      if (this.plotCovariates) {
        this.updateCovariatesPlot();
      }
      this.filterChanged = false;
    }
  }
}