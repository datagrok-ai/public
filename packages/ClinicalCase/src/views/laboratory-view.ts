import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { ClinRow, study } from "../clinical-study";
import { createBaselineEndpointDataframe, createHysLawDataframe, createLabValuesByVisitDataframe } from '../data-preparation/data-preparation';
import { ALT, BILIRUBIN } from '../constants';
import { createBaselineEndpointScatterPlot, createHysLawScatterPlot } from '../custom-scatter-plots/custom-scatter-plots';
import { updateDivInnerHTML } from './utils';
import { _package } from '../package';
import { getUniqueValues } from '../data-preparation/utils';
import { LAB_HI_LIM_N, LAB_LO_LIM_N, LAB_TEST, VISIT_DAY, VISIT_NAME, SUBJECT_ID, TREATMENT_ARM, LAB_RES_N } from '../columns-constants';
import { ClinicalCaseViewBase } from '../model/ClinicalCaseViewBase';

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

  constructor(name) {
    super({});
    this.name = name;
    this.helpUrl = `${_package.webRoot}/views_help/laboratory.md`;
  }

  createView(): void {
    this.lb = study.domains.lb.clone();
    this.dm = study.domains.dm.clone();

    this.uniqueLabValues = Array.from(getUniqueValues(this.lb, LAB_TEST));
    this.uniqueVisits = Array.from(getUniqueValues(this.lb, VISIT_NAME));
    this.uniqueTreatmentArms = Array.from(getUniqueValues(this.dm, TREATMENT_ARM));
    this.selectedLabBlEp = this.uniqueLabValues[0];
    this.selectedBl = this.uniqueVisits[0];
    this.selectedEp = this.uniqueVisits[1];
    this.selectedLabDistr = this.uniqueLabValues[0];
    this.selectedArm = this.uniqueTreatmentArms[0];

    let grid = this.lb.plot.grid();
    this.lb.onCurrentRowChanged.subscribe((_) => {
      grok.shell.o = new ClinRow(this.lb.currentRow);
    });

    this.root.className = 'grok-view ui-box';

    let tabControl = ui.tabControl(null, false);

    let hysLawGuide = ui.info('Please select values for ALT/AST and Bilirubin in a property panel', '', false);
    updateDivInnerHTML(this.hysLawDiv, hysLawGuide);

    tabControl.addPane('Hy\'s law', () => this.hysLawDiv);
    
    tabControl.addPane('Baseline endpoint', () => this.baselineEndpointDiv);

    tabControl.addPane('Laboratory distribution', () => this.distributionDiv);

    tabControl.addPane('Results', () => grid.root);

    this.root.appendChild(
      tabControl.root
    );
    this.updateBaselineEndpointPlot();
    this.updateDistributionPlot();

  }


  updateHysLawScatterPlot() {
    if (this.selectedALT && this.selectedAST && this.selectedBLN) {
      this.createHysLawScatterPlot();
      updateDivInnerHTML(this.hysLawDiv, this.hysLawScatterPlot.root);
    }
  }

  private createHysLawScatterPlot() {
    let hysLawDataframe = createHysLawDataframe(this.lb, this.dm, this.selectedALT, this.selectedAST, this.selectedBLN);
    grok.data.linkTables(this.lb, hysLawDataframe,
      [SUBJECT_ID], [SUBJECT_ID],
      [DG.SYNC_TYPE.CURRENT_ROW_TO_ROW, DG.SYNC_TYPE.CURRENT_ROW_TO_SELECTION]);
    grok.data.linkTables(hysLawDataframe, this.lb,
      [SUBJECT_ID], [SUBJECT_ID],
      [DG.SYNC_TYPE.SELECTION_TO_SELECTION, DG.SYNC_TYPE.SELECTION_TO_SELECTION]);
    this.hysLawScatterPlot = createHysLawScatterPlot(hysLawDataframe, ALT, BILIRUBIN, TREATMENT_ARM);
  }

  updateBaselineEndpointPlot() {
    let visitCol = VISIT_NAME;
    let blNumCol = `${this.selectedLabBlEp}_BL`;
    let epNumCol = `${this.selectedLabBlEp}_EP`;
    let baselineEndpointDataframe = createBaselineEndpointDataframe(this.lb, this.dm, [TREATMENT_ARM], LAB_TEST, LAB_RES_N,
      [LAB_LO_LIM_N, LAB_HI_LIM_N], this.selectedLabBlEp, this.selectedBl, this.selectedEp, visitCol, blNumCol, epNumCol);
    grok.data.linkTables(this.lb, baselineEndpointDataframe,
      [SUBJECT_ID], [SUBJECT_ID],
      [DG.SYNC_TYPE.CURRENT_ROW_TO_ROW, DG.SYNC_TYPE.CURRENT_ROW_TO_SELECTION]);
    this.baselineEndpointPlot = createBaselineEndpointScatterPlot(baselineEndpointDataframe, blNumCol, epNumCol, TREATMENT_ARM,
      baselineEndpointDataframe.get(LAB_LO_LIM_N, 0), baselineEndpointDataframe.get(LAB_HI_LIM_N, 0));
    updateDivInnerHTML(this.baselineEndpointDiv, this.baselineEndpointPlot.root);
  }

  updateDistributionPlot() {
    let labValue = this.selectedLabDistr;
    let labValueNumColumn = `${labValue} values`;
    let disributionDataframe = createLabValuesByVisitDataframe(this.lb, this.dm, labValue, this.selectedArm, labValueNumColumn, VISIT_DAY);
    this.distributionPlot = DG.Viewer.boxPlot(disributionDataframe, {
      category: VISIT_DAY,
      value: labValueNumColumn,
    });
    updateDivInnerHTML(this.distributionDiv, this.distributionPlot.root);
  }

  override async propertyPanel() {
    let panelDiv = ui.div();
    let acclab = this.createAccWithTitle(this.name);

    let altChoices = ui.choiceInput('ALT', this.selectedALT, this.uniqueLabValues);
    altChoices.onChanged((v) => {
      this.selectedALT = altChoices.value;
      this.updateHysLawScatterPlot();
    });
    //@ts-ignore
    altChoices.input.style.width = '150px';

    let astChoices = ui.choiceInput('AST', this.selectedAST, this.uniqueLabValues);
    astChoices.onChanged((v) => {
      this.selectedAST = astChoices.value;
      this.updateHysLawScatterPlot();
    });
    //@ts-ignore
    astChoices.input.style.width = '150px';

    let blnChoices = ui.choiceInput('BLN', this.selectedBLN, this.uniqueLabValues);
    blnChoices.onChanged((v) => {
      this.selectedBLN = blnChoices.value;
      this.updateHysLawScatterPlot();
    });
    //@ts-ignore
    blnChoices.input.style.width = '150px';
    acclab.addPane('Hy\'s law', () => ui.divV([altChoices.root, astChoices.root, blnChoices.root]), true);

    let labValueChoices = ui.choiceInput('Value', this.selectedLabBlEp, this.uniqueLabValues);
    let blVisitChoices = ui.choiceInput('BL', this.selectedBl, this.uniqueVisits);
    let epVisitChoices = ui.choiceInput('EP', this.selectedEp, this.uniqueVisits);
    labValueChoices.onChanged((v) => {
      this.selectedLabBlEp = labValueChoices.value;
      this.updateBaselineEndpointPlot();
    });
    blVisitChoices.onChanged((v) => {
      this.selectedBl = blVisitChoices.value;
      this.updateBaselineEndpointPlot();
    });
    epVisitChoices.onChanged((v) => {
      this.selectedEp = epVisitChoices.value;
      this.updateBaselineEndpointPlot();
    });
    //@ts-ignore
    labValueChoices.input.style.width = '150px';
    acclab.addPane('Baseline endpoint', () => ui.divV([labValueChoices.root, blVisitChoices.root, epVisitChoices.root]), true);


    let distrChoices = ui.choiceInput('Value', this.selectedLabDistr, this.uniqueLabValues);
    let treatmentArmsChoices = ui.choiceInput('Treatment arm', this.selectedArm, this.uniqueTreatmentArms);
    distrChoices.onChanged((v) => {
      this.selectedLabDistr = distrChoices.value;
      this.updateDistributionPlot();
    });
    treatmentArmsChoices.onChanged((v) => {
      this.selectedArm = treatmentArmsChoices.value;
      this.updateDistributionPlot();
    });
    //@ts-ignore
    distrChoices.input.style.width = '150px';
    acclab.addPane('Laboratory distribution', () => ui.divV([distrChoices.root, treatmentArmsChoices.root]), true);
    panelDiv.append(acclab.root);
    return panelDiv;
  }

}