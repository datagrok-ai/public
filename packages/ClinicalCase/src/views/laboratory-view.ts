import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {ClinRow} from '../clinical-study';
import {studies} from '../package';
import {createBaselineEndpointDataframe, createHysLawDataframe,
  createLabValuesByVisitDataframe,
  createVisitDayStrCol} from '../data-preparation/data-preparation';
import {ALT, BILIRUBIN} from '../constants/constants';
import {createBaselineEndpointScatterPlot, createHysLawScatterPlot} from '../custom-scatter-plots/custom-scatter-plots';
import {updateDivInnerHTML} from '../utils/utils';
import {_package} from '../package';
import {getUniqueValues} from '../data-preparation/utils';
import {LAB_HI_LIM_N, LAB_LO_LIM_N, LAB_TEST, VISIT_DAY,
  SUBJECT_ID, LAB_RES_N,
  VISIT_DAY_STR,
  VISIT} from '../constants/columns-constants';
import {ClinicalCaseViewBase} from '../model/ClinicalCaseViewBase';
import {TRT_ARM_FIELD} from '../views-config';
import {checkColumnsAndCreateViewer} from '../utils/views-validation-utils';
import {CDISC_STANDARD} from '../utils/types';

const resultsTabName = 'Results';
const distrTabName = 'Laboratory distribution';
const hysLawTabName = 'Hy\'s law';
const baselineEndpointTabName = 'Baseline endpoint';
export class LaboratoryView extends ClinicalCaseViewBase {
  hysLawDiv = ui.box();
  baselineEndpointDiv = ui.box();
  distributionDiv = ui.box();
  uniqueLabValues: any;
  uniqueVisits: any;
  uniqueTreatmentArms: any;
  selectedALT = '';
  selectedAST = '';
  selectedBLN = '';
  selectedLabBlEp = '';
  selectedBl = '';
  selectedEp = '';
  selectedLabDistr = '';
  selectedArm = '';
  hysLawScatterPlot: DG.ScatterPlotViewer;
  baselineEndpointPlot: DG.ScatterPlotViewer;
  distributionPlot: DG.Viewer;
  dm: DG.DataFrame;
  lb: DG.DataFrame;
  selectedTab: string;
  isSend = false;
  propPanelPanes: DG.AccordionPane[] = [];
  numVisDayColDict: {[key: string]: string} = {'lb': VISIT_DAY};

  constructor(name, studyId) {
    super(name, studyId);
    this.name = name;
    this.helpUrl = `${_package.webRoot}/views_help/laboratory.md`;
  }

  createView(): void {
    this.isSend = studies[this.studyId].config.standard === CDISC_STANDARD.SEND;
    if (this.isSend)
      createVisitDayStrCol(studies[this.studyId].domains.lb, this.numVisDayColDict);
    this.lb = studies[this.studyId].domains.lb.clone();
    if (studies[this.studyId].domains.dm)
      this.dm = studies[this.studyId].domains.dm.clone();


    this.uniqueLabValues = this.lb.col(LAB_TEST) ? Array.from(getUniqueValues(this.lb, LAB_TEST)) : [];
    this.uniqueVisits = this.lb.col(this.isSend ? VISIT_DAY_STR : VISIT) ?
      Array.from(getUniqueValues(this.lb, this.isSend ? VISIT_DAY_STR : VISIT)) : [];
    this.uniqueTreatmentArms = this.dm &&
      this.dm.col(studies[this.studyId].viewsConfig.config[this.name][TRT_ARM_FIELD]) ?
      Array.from(getUniqueValues(this.dm, studies[this.studyId].viewsConfig.config[this.name][TRT_ARM_FIELD])) : [];
    this.uniqueTreatmentArms = [''].concat(this.uniqueTreatmentArms);
    this.selectedLabBlEp = this.uniqueLabValues.length ? this.uniqueLabValues[0] : null;
    this.selectedBl = this.uniqueVisits.length ? this.uniqueVisits[0] : null;
    this.selectedEp = this.uniqueVisits.length ? this.uniqueVisits[1] : null;
    this.selectedLabDistr = this.uniqueLabValues.length ? this.uniqueLabValues[0] : null;
    this.selectedArm = this.uniqueTreatmentArms.length ? this.uniqueTreatmentArms[0] : null;

    const grid = this.lb.plot.grid();
    this.lb.onCurrentRowChanged.subscribe((_) => {
      grok.shell.o = new ClinRow(this.lb.currentRow);
    });

    // if (studies[this.studyId].domains.dm) {
    //   grok.data.linkTables(studies[this.studyId].domains.dm, this.lb,
    //     [SUBJECT_ID], [SUBJECT_ID],
    //     [DG.SYNC_TYPE.FILTER_TO_FILTER]);
    // }

    this.root.className = 'grok-view ui-box';

    const tabControl = ui.tabControl(null, false);

    const distributionCreated = checkColumnsAndCreateViewer(
      studies[this.studyId].domains.lb,
      [SUBJECT_ID, LAB_TEST, LAB_RES_N, this.numVisDayColDict['lb']],
      this.distributionDiv, () => {
        this.updateDistributionPlot();
      },
      distrTabName);

    tabControl.addPane(distrTabName, () => this.distributionDiv);

    const resultsPane = tabControl.addPane(resultsTabName, () => grid.root);

    checkColumnsAndCreateViewer(
      studies[this.studyId].domains.lb,
      [SUBJECT_ID, LAB_TEST, LAB_RES_N, this.isSend ? VISIT_DAY_STR : VISIT,
        LAB_LO_LIM_N, LAB_HI_LIM_N],
      this.baselineEndpointDiv, () => {
        this.updateBaselineEndpointPlot();
      },
      baselineEndpointTabName);

    tabControl.addPane(baselineEndpointTabName, () => this.baselineEndpointDiv);

    //SDTM specific plots (for human clinical trials)
    if (!this.isSend) {
      const hysLawGuide = ui.info('Please select values for ALT/AST and Bilirubin in a context panel', '', false);
      checkColumnsAndCreateViewer(
        studies[this.studyId].domains.lb,
        [SUBJECT_ID, LAB_RES_N, LAB_HI_LIM_N, LAB_TEST],
        this.hysLawDiv, () => {
          updateDivInnerHTML(this.hysLawDiv, hysLawGuide);
        },
        hysLawTabName);

      tabControl.addPane(hysLawTabName, () => this.hysLawDiv);
    }

    tabControl.onTabChanged.subscribe((tab: DG.TabPane) => {
      this.propPanelPanes.forEach((it) => it.root.style.display = it.name === tab.name ? 'flex' : 'none');
    });

    //this is done to show valid tab initially
    if (!distributionCreated)
      tabControl.currentPane = resultsPane;

    this.root.appendChild(
      tabControl.root,
    );
  }

  updateHysLawScatterPlot() {
    if (this.selectedALT && this.selectedAST && this.selectedBLN) {
      this.createHysLawScatterPlot();
      updateDivInnerHTML(this.hysLawDiv, this.hysLawScatterPlot.root);
    }
  }

  private createHysLawScatterPlot() {
    const hysLawDataframe = createHysLawDataframe(this.lb, this.dm, this.selectedALT,
      this.selectedAST, this.selectedBLN, studies[this.studyId].viewsConfig.config[this.name][TRT_ARM_FIELD]);
    // if (studies[this.studyId].domains.dm) {
    //   grok.data.linkTables(studies[this.studyId].domains.dm, hysLawDataframe,
    //     [SUBJECT_ID], [SUBJECT_ID],
    //     [DG.SYNC_TYPE.FILTER_TO_FILTER]);
    // }
    this.hysLawScatterPlot = createHysLawScatterPlot(hysLawDataframe, ALT,
      BILIRUBIN, studies[this.studyId].viewsConfig.config[this.name][TRT_ARM_FIELD]);
  }

  updateBaselineEndpointPlot() {
    const visitCol = this.isSend ? VISIT_DAY_STR : VISIT;
    const blNumCol = `${this.selectedLabBlEp}_BL`;
    const epNumCol = `${this.selectedLabBlEp}_EP`;
    const baselineEndpointDataframe = createBaselineEndpointDataframe(this.lb, this.dm,
      [studies[this.studyId].viewsConfig.config[this.name][TRT_ARM_FIELD]], LAB_TEST, LAB_RES_N,
      [LAB_LO_LIM_N, LAB_HI_LIM_N], this.selectedLabBlEp, this.selectedBl,
      this.selectedEp, visitCol, blNumCol, epNumCol);
    // if (studies[this.studyId].domains.dm) {
    //   grok.data.linkTables(studies[this.studyId].domains.dm, baselineEndpointDataframe,
    //     [SUBJECT_ID], [SUBJECT_ID],
    //     [DG.SYNC_TYPE.FILTER_TO_FILTER]);
    // }
    if (!baselineEndpointDataframe.rowCount) {
      updateDivInnerHTML(this.baselineEndpointDiv,
        ui.info(`No data found for ${this.selectedLabBlEp} and selected visits`));
      return;
    }
    this.baselineEndpointPlot = createBaselineEndpointScatterPlot(baselineEndpointDataframe, blNumCol,
      epNumCol, studies[this.studyId].viewsConfig.config[this.name][TRT_ARM_FIELD],
      baselineEndpointDataframe.get(LAB_LO_LIM_N, 0), baselineEndpointDataframe.get(LAB_HI_LIM_N, 0));
    updateDivInnerHTML(this.baselineEndpointDiv, this.baselineEndpointPlot.root);
  }

  updateDistributionPlot() {
    const labValue = this.selectedLabDistr;
    const labValueNumColumn = `${labValue} values`;
    const disributionDataframe = createLabValuesByVisitDataframe(this.lb, this.dm, labValue,
      studies[this.studyId].viewsConfig.config[this.name][TRT_ARM_FIELD],
      this.selectedArm, labValueNumColumn, this.numVisDayColDict['lb']);
    // if (studies[this.studyId].domains.dm) {
    //   grok.data.linkTables(studies[this.studyId].domains.dm, disributionDataframe,
    //     [SUBJECT_ID], [SUBJECT_ID],
    //     [DG.SYNC_TYPE.FILTER_TO_FILTER]);
    // }
    this.distributionPlot = DG.Viewer.boxPlot(disributionDataframe, {
      categoryColumnNames: [this.numVisDayColDict['lb']],
      value: labValueNumColumn,
    });
    this.distributionPlot.setOptions({category: this.numVisDayColDict['lb']});
    updateDivInnerHTML(this.distributionDiv, this.distributionPlot.root);
  }

  override async propertyPanel() {
    const panelDiv = ui.div();
    const acclab = this.createAccWithTitle(this.name);

    //distributions panel
    const distrChoices = ui.input.choice('Value', {value: this.selectedLabDistr, items: this.uniqueLabValues});
    const treatmentArmsChoices = ui.input.choice('Treatment arm',
      {value: this.selectedArm, items: this.uniqueTreatmentArms});
    distrChoices.onChanged.subscribe((value) => {
      this.selectedLabDistr = value;
      this.updateDistributionPlot();
    });
    treatmentArmsChoices.onChanged.subscribe((value) => {
      this.selectedArm = value;
      this.updateDistributionPlot();
    });
    distrChoices.input.style.width = '150px';
    this.propPanelPanes.push(acclab
      .addPane(distrTabName, () => ui.divV([distrChoices.root, treatmentArmsChoices.root]), true));

    const labValueChoices = ui.input.choice('Value', {value: this.selectedLabBlEp, items: this.uniqueLabValues});
    const blVisitChoices = ui.input.choice('BL', {value: this.selectedBl, items: this.uniqueVisits});
    const epVisitChoices = ui.input.choice('EP', {value: this.selectedEp, items: this.uniqueVisits});
    labValueChoices.onChanged.subscribe((value) => {
      this.selectedLabBlEp = value;
      this.updateBaselineEndpointPlot();
    });
    blVisitChoices.onChanged.subscribe((value) => {
      this.selectedBl = value;
      this.updateBaselineEndpointPlot();
    });
    epVisitChoices.onChanged.subscribe((value) => {
      this.selectedEp = value;
      this.updateBaselineEndpointPlot();
    });
    labValueChoices.input.style.width = '150px';
    const bePane = acclab.addPane(baselineEndpointTabName,
      () => ui.divV([labValueChoices.root, blVisitChoices.root, epVisitChoices.root]), true);
    //set invisible by default
    bePane.root.style.display = 'none';
    this.propPanelPanes.push(bePane);


    //SDTM specific options
    if (!this.isSend) {
      //Hy's law panel
      const altChoices = ui.input.choice('ALT', {value: this.selectedALT, items: this.uniqueLabValues});
      altChoices.onChanged.subscribe((value) => {
        this.selectedALT = value;
        this.updateHysLawScatterPlot();
      });
      altChoices.input.style.width = '150px';

      const astChoices = ui.input.choice('AST', {value: this.selectedAST, items: this.uniqueLabValues});
      astChoices.onChanged.subscribe((value) => {
        this.selectedAST = value;
        this.updateHysLawScatterPlot();
      });
      astChoices.input.style.width = '150px';

      const blnChoices = ui.input.choice('BLN', {value: this.selectedBLN, items: this.uniqueLabValues});
      blnChoices.onChanged.subscribe((value) => {
        this.selectedBLN = value;
        this.updateHysLawScatterPlot();
      });
      blnChoices.input.style.width = '150px';
      const hysLawPane = acclab
        .addPane(hysLawTabName, () => ui.divV([altChoices.root, astChoices.root, blnChoices.root]), true);
      //set invisible by default
      hysLawPane.root.style.display = 'none';
      this.propPanelPanes.push(hysLawPane);
    }

    panelDiv.append(acclab.root);
    return panelDiv;
  }
}
