import * as DG from "datagrok-api/dg";
import { InputBase } from "datagrok-api/dg";
import * as grok from 'datagrok-api/grok';
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";
import { TREATMENT_ARM } from "../constants";
import { createSurvivalData } from "../data-preparation/data-preparation";
import { dataframeContentToRow } from "../data-preparation/utils";

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

    let applyFilters = ui.button('Apply to curves', () => { 
      this.updateSurvivalPlot();
      if (this.plotCovariates) {
        this.updateCovariatesPlot();
      }
    });

    let createSurvivalDataframe = ui.bigButton('Create dataset', () => { 
     this.refreshDataframe();
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
            guide.innerHTML = '';
            guide.append(ui.info(`1. Select dataset paramenters and click 'Create dataset'
            2. Filter data and apply them to curves. Click on 'Apply to curves',
            3. Set the survival chart parameters. To see the chart clicn on 'Survival Chart' tab on right.
            4. Set the co-variates. On the right click on 'Co-Variates' tab for see them.
            `,'Survival Analysis Quick Guide', false))
          })
        ]
      ])

      let guide = ui.info(`1. Select dataset paramenters and click 'Create dataset'
      2. Filter data and apply them to curves. Click on 'Apply to curves',
      3. Set the survival chart parameters. To see the chart clicn on 'Survival Chart' tab on right.
      4. Set the co-variates. On the right click on 'Co-Variates' tab for see them.
      `,'Survival Analysis Quick Guide', false);
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
          ui.tabControl({
            'Dataset': ui.splitV([this.survivalFilterDiv, this.survivalGridDivCreate]),
            'Survival Chart': ()=>{
              return this.survivalPlotDiv
            },
            'Co-Variates': ()=>{
              return this.covariatesPlotDiv
            }
          }).root
        ])
      ]))
      //@ts-ignore
      guide.parentNode.style.flexGrow = '0';
      //@ts-ignore
      guide.parentNode.classList = 'ui-div';

      /*
      this.root.append( ui.tabControl({
        '1. Dataset': ui.splitH([
          ui.box(ui.panel([
            ui.inputs([ 
              endpointChoices,
              covariatesChoices,
              //@ts-ignore
              ui.buttonsInput([createSurvivalDataframe, applyFilters])
            ])
          ]), { style: { maxWidth: '300px' }}),
          this.survivalGridDivCreate,
          this.survivalFilterDiv,
        ]),
        '2. Results': ui.splitV([
          ui.splitH([
            ui.box(ui.panel([
              ui.h1('Survival data'),
              ui.inputs([ 
                confIntChoices,
                //@ts-ignore
                this.strataChoicesDiv
              ])
            ])),
            ui.box(ui.panel([
              ui.h1('Co-variates'),
              this.plotCovariatesChoicesDiv
            ])),
         ] , { style: { maxHeight: '150px' } }),
          ui.splitH([
            this.survivalPlotDiv,
            this.covariatesPlotDiv 
          ])
        ]),
      }).root)
      */
      /*
      this.root.append(
        ui.tabControl({
          'Create dataset':
            ui.splitH([
              ui.box(ui.panel([
                ui.inputs([ endpointChoices,
                covariatesChoices,
                //@ts-ignore
                ui.buttonsInput([ createSurvivalDataframe ])
              ])
              ]), { style: { maxWidth: '300px' } }),
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
              ui.box(ui.panel([
                this.plotCovariatesChoicesDiv
              ]), { style: { maxWidth: '180px' } }),
              this.covariatesPlotDiv ])
        }).root
      );*/

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
      let img = ui.image(`data:image/png;base64,${result[ 'plot' ]}`);
      //img.style.backgroundPosition = 'top'
      this.updateDivInnerHTML(this.survivalPlotDiv, img);
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
      let img = ui.image(`data:image/png;base64,${result[ 'plot' ]}`);
      //img.style.backgroundPosition = 'top'
      this.updateDivInnerHTML(this.covariatesPlotDiv, img);
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
     this.updateDivInnerHTML(this.strataChoicesDiv, this.strataChoices.root);
     this.updateDivInnerHTML(this.plotCovariatesChoicesDiv, this.plotCovariatesChoices.root);
     this.updateDivInnerHTML(this.survivalGridDivCreate, this.survivalDataframe.plot.grid().root);
     this.updateDivInnerHTML(this.survivalGridDivFilter, this.survivalDataframe.plot.grid().root);
     this.updateDivInnerHTML(this.survivalFilterDiv, this.getFilters());
     this.updateSurvivalPlot();
  }
}