import * as DG from "datagrok-api/dg";
import { InputBase } from "datagrok-api/dg";
import * as grok from 'datagrok-api/grok';
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";
import { AE_CAUSALITY, AE_REQ_HOSP, AE_SEQ, AE_SEVERITY, AGE, DEATH_DATE, RACE, SEX, SUBJECT_ID, SUBJ_REF_ENDT } from "../constants/columns-constants";
import { SURVIVAL_ANALYSIS_GUIDE } from "../constants/constants";
import { createSurvivalData } from "../data-preparation/data-preparation";
import { dataframeContentToRow } from "../data-preparation/utils";
import { ClinicalCaseViewBase } from "../model/ClinicalCaseViewBase";
import { _package } from "../package";
import { SURVIVAL_ANALYSIS_VIEW_NAME } from "../constants/view-names-constants";
import { AE_START_DAY_FIELD, TRT_ARM_FIELD, VIEWS_CONFIG } from "../views-config";
import { updateDivInnerHTML } from "../utils/utils";

export class SurvivalAnalysisView extends ClinicalCaseViewBase {

  survivalPlotDiv = ui.box();
  covariatesPlotDiv = ui.box();
  survivalGridDivCreate = ui.box();
  survivalFilterDiv = ui.box(ui.divText('Create dataset',{style:{color:'var(--grey-3)',marginTop:'30px', alignItems:'center'}}));
  strataChoicesDiv = ui.div();
  plotCovariatesChoicesDiv = ui.div();
  parametersPanel = ui.div();
  strataChoices: DG.InputBase;
  plotCovariatesChoices: HTMLDivElement;
  endpointChoices: DG.InputBase;
  covariatesChoices: DG.InputBase;
  confIntChoices: DG.InputBase;
  createSurvivalDataframe: HTMLButtonElement;
  survivalColumns = [];
  confIntervals = [ 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99 ];
  survivalOptions = [''];
  covariatesOptions: any;
  endpointOptions = { 'RETAIN IN STUDY': SUBJ_REF_ENDT };
  confInterval = 0.7;
  strata = '';
  endpoint = '';
  covariates = [];
  survivalDataframe: DG.DataFrame;
  plotCovariates = [];
  filterChanged = false;
  colsRequiredForEndpoints: any;

  constructor(name) {
    super({});
    this.name = name;
    this.helpUrl = `${_package.webRoot}/views_help/survival_analysis.md`;
  }

  createView(): void {
    this.colsRequiredForEndpoints = this.getColsRequiredForEndpoints();
    this.updateEndpointOptions();
    this.covariatesOptions = [ AGE, SEX, RACE, VIEWS_CONFIG[SURVIVAL_ANALYSIS_VIEW_NAME][TRT_ARM_FIELD] ].filter(it => study.domains.dm.columns.names().includes(it));
    this.endpoint = Object.keys(this.endpointOptions)[0];
    this.endpointChoices = ui.choiceInput('Endpoint', Object.keys(this.endpointOptions)[0], Object.keys(this.endpointOptions));
    this.endpointChoices.onChanged((v) => {
      this.endpoint = this.endpointChoices.value;
    });

    this.covariatesChoices = ui.multiChoiceInput('Covariates', null, this.covariatesOptions);
    this.covariatesChoices.onChanged((v) => {
      this.covariates = this.covariatesChoices.value;
    });

    this.confIntChoices = ui.choiceInput('Confidence', this.confIntervals[0], this.confIntervals);
    this.confIntChoices.onChanged((v) => {
      this.confInterval = this.confIntChoices.value;
      this.updateSurvivalPlot();
    });

    this.updateStrataChoices();
    this.updatePlotCovariatesChoices();

    this.createSurvivalDataframe = ui.bigButton('Create dataset', () => {
      this.refreshDataframe();
      this.filterChanged = true;
      this.updateChartsAfterFiltering();
      this.survivalDataframe.onFilterChanged.subscribe((_) => {
        this.filterChanged = true;
      });
    });

    let guide = ui.info(SURVIVAL_ANALYSIS_GUIDE, 'Survival Analysis Quick Guide', false);

    let tabControl = ui.tabControl(null, false);
    tabControl.addPane('Dataset', () =>
      ui.splitV([
        ui.splitH([
          ui.box(ui.panel([
            ui.inputs([
              this.endpointChoices,
              this.covariatesChoices,
              //@ts-ignore
              ui.buttonsInput([this.createSurvivalDataframe])
            ])
          ]), { style: { maxWidth: '300px' } }),
          ui.splitV([this.survivalFilterDiv, this.survivalGridDivCreate])
        ])
      ]));

    tabControl.addPane('Survival Chart', () => ui.splitV([
      ui.splitH([
        ui.box(ui.panel([
          ui.inputs([
            this.confIntChoices,
            //@ts-ignore
            this.strataChoicesDiv,
          ])
        ]), { style: { maxWidth: '300px' } }),
        this.survivalPlotDiv
      ])
    ]));
    tabControl.getPane('Survival Chart').header.addEventListener('click', () => {
      this.updateChartsAfterFiltering();

    });

    tabControl.addPane('Covariates', () => ui.splitV([
      ui.box(
        //@ts-ignore
        this.plotCovariatesChoicesDiv,
        { style: { maxHeight: '50px' } }),
      this.covariatesPlotDiv
    ]));
    tabControl.getPane('Covariates').header.addEventListener('click', () => {
      this.updateChartsAfterFiltering();
    });
    this.root.className = 'grok-view ui-box';
    this.setRibbonPanels([
      [
        ui.icons.info(() => {
          updateDivInnerHTML(guide, ui.info(SURVIVAL_ANALYSIS_GUIDE, 'Survival Analysis Quick Guide', false));
        })
      ]
    ])

    this.root.append(ui.splitV([
      guide,
      tabControl.root
    ]))
    //@ts-ignore
    guide.parentNode.style.flexGrow = '0';
    //@ts-ignore
    guide.parentNode.classList = 'ui-div';

    // this.updateParameterPanel('Dataset');
  }

  private updateEndpointOptions(){
    Object.keys(this.colsRequiredForEndpoints).forEach(key => {
      let missingCols = false;
      let domains = Object.keys(this.colsRequiredForEndpoints[key]);
      domains.forEach(dom => {
        this.colsRequiredForEndpoints[key][dom].map(it => {
          if (!study.domains[dom] || !study.domains[dom].columns.names().includes(it)) {
            missingCols = true;
          }
        })
        if(!missingCols){
          let domain = Object.keys(this.colsRequiredForEndpoints[key])[0]
          this.endpointOptions[key] = this.colsRequiredForEndpoints[key][domain][0];
        }
      })
    })
  }

  private updateSurvivalPlot() {
    ui.setUpdateIndicator(this.survivalPlotDiv, true);
    grok.functions.call(
      "Clinicalcase:survivalPlot", {
      "survivalDf": this.survivalDataframe,
      "inputStrata": this.strata,
      "confInt": this.confInterval.toString()
    }).then((result) => {
      ui.setUpdateIndicator(this.survivalPlotDiv, false);
      //@ts-ignore
      updateDivInnerHTML(this.survivalPlotDiv, ui.image(`data:image/png;base64,${result[ 'plot' ]}`));
      console.warn(dataframeContentToRow(result[ 'diagnostics' ]));
    });
  }

  private updateCovariatesPlot() {
    ui.setUpdateIndicator(this.covariatesPlotDiv, true);
    grok.functions.call(
      "Clinicalcase:covariates", {
      "covariatesDf": this.survivalDataframe,
      "coVariates": this.plotCovariates.join(' + ')
    }).then((result) => {
      ui.setUpdateIndicator(this.covariatesPlotDiv, false);
      //@ts-ignore
      let img = ui.image(`data:image/png;base64,${result['plot']}`);
      img.style.backgroundSize = '100% 100%';
      updateDivInnerHTML(this.covariatesPlotDiv, img);
      console.warn(dataframeContentToRow(result['diagnostics']));
    });
  }

  private getFilters(){
    return DG.Viewer.fromType('Filters', this.survivalDataframe, {
      'columnNames': this.survivalColumns.filter(it => it !== 'time' && it !== 'status' && it !== SUBJECT_ID),
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
    this.plotCovariatesChoices = ui.divH([], {style: {'margin-left': '20px'}});
    this.covariates.forEach(it => {
      let covariateCheckbox = ui.boolInput(`${it}`, false);
      this.plotCovariatesChoices.append(covariateCheckbox.root);
      covariateCheckbox.onChanged((v) => {
        covariateCheckbox.value ? 
        this.plotCovariates.push(covariateCheckbox.caption) : 
        this.plotCovariates = this.plotCovariates.filter(it => it!== covariateCheckbox.caption);
        this.updateCovariatesPlot();
      });
    })
/*     this.plotCovariatesChoices = ui.multiChoiceInput(' ', null, this.survivalColumns.filter(it => it !== 'time' && it !== 'status' && it !== SUBJECT_ID));
    this.plotCovariatesChoices.onChanged((v) => {
      this.plotCovariates = this.plotCovariatesChoices.value;
      this.updateCovariatesPlot();
    }); */
  }

  private refreshDataframe(){
     this.survivalDataframe = createSurvivalData(study.domains.dm, study.domains.ae ?? null, this.endpoint, this.endpointOptions[this.endpoint], this.covariates);
     this.survivalColumns = this.survivalDataframe.columns.names();
     this.survivalOptions = [''].concat(this.survivalColumns.filter(it => it !== 'time' && it !== 'status' && it !== SUBJECT_ID));
     this.strata = '';
     this.plotCovariates = [];
     this.updateStrataChoices();
     this.updatePlotCovariatesChoices();
     updateDivInnerHTML(this.strataChoicesDiv, this.strataChoices.root);
     updateDivInnerHTML(this.plotCovariatesChoicesDiv, this.plotCovariatesChoices);
     updateDivInnerHTML(this.survivalGridDivCreate, this.survivalDataframe.plot.grid().root);
     updateDivInnerHTML(this.survivalFilterDiv, this.getFilters());
     this.updateSurvivalPlot();
  }

  private updateChartsAfterFiltering(){
    if(this.filterChanged){
      updateDivInnerHTML(this.covariatesPlotDiv, '');
      updateDivInnerHTML(this.survivalPlotDiv, '');
      this.updateSurvivalPlot();
      if (this.plotCovariates.length) {
        this.updateCovariatesPlot();
      }
      this.filterChanged = false;
    }
  }

  private getColsRequiredForEndpoints() {
    return {
      'SAE': { 'ae': [VIEWS_CONFIG[SURVIVAL_ANALYSIS_VIEW_NAME][AE_START_DAY_FIELD], SUBJECT_ID, AE_SEVERITY, AE_SEQ] },
      'DEATH': { 'dm': [DEATH_DATE] },
      'HOSPITALIZATION': { 'ae': [VIEWS_CONFIG[SURVIVAL_ANALYSIS_VIEW_NAME][AE_START_DAY_FIELD], SUBJECT_ID, AE_REQ_HOSP, AE_SEVERITY, AE_SEQ] },
      'DRUG RELATED AE': { 'ae': [VIEWS_CONFIG[SURVIVAL_ANALYSIS_VIEW_NAME][AE_START_DAY_FIELD], SUBJECT_ID, AE_CAUSALITY, AE_SEVERITY, AE_SEQ] },
    }
  }
}