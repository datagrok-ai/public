import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { ClinRow, study } from "../clinical-study";
import { createBaselineEndpointDataframe, createHysLawDataframe, createLabValuesByVisitDataframe } from '../data-preparation/data-preparation';
import { ALT, AP, BILIRUBIN, requiredColumnsByView, } from '../constants';
import { createBaselineEndpointScatterPlot, createHysLawScatterPlot } from '../custom-scatter-plots/custom-scatter-plots';
import { ILazyLoading } from '../lazy-loading/lazy-loading';
import { checkMissingDomains } from './utils';
import { _package } from '../package';
import { getUniqueValues } from '../data-preparation/utils';
import { LAB_HI_LIM_N, LAB_LO_LIM_N, LAB_TEST, LAB_VISIT_DAY, LAB_VISIT_NAME, SUBJECT_ID, TREATMENT_ARM } from '../columns-constants';

export class LaboratoryView extends DG.ViewBase implements ILazyLoading {

  constructor(name) {
    super({});
    this.name = name;
    this.helpUrl = `${_package.webRoot}/views_help/laboratory.md`;
 }
  loaded: boolean;

  load(): void {
    checkMissingDomains(requiredColumnsByView[this.name], false, this);
 }
  
  createView(): void {
    let lb = study.domains.lb;
    let dm = study.domains.dm;

    let grid = lb.plot.grid();
    lb.onCurrentRowChanged.subscribe((_) => {
      grok.shell.o = new ClinRow(lb.currentRow);
    });

    let hysLawScatterPlot = this.hysLawScatterPlot(dm, lb);

    let uniqueLabValues = Array.from(getUniqueValues(lb, LAB_TEST));
    let uniqueVisits = Array.from(getUniqueValues(lb, LAB_VISIT_NAME));
    let baselineEndpointPlot = this.baselineEndpointPlot(dm, lb, uniqueLabValues[ 0 ], uniqueVisits[ 0 ], uniqueVisits[ 1 ]);

    let uniqueTreatmentArms = Array.from(getUniqueValues(dm, TREATMENT_ARM));
    let disributionBoxPlot = this.labValuesDistributionPlot(dm, lb, uniqueLabValues[ 0 ], uniqueTreatmentArms[ 0 ]);

    this.generateUI(dm, lb, grid, hysLawScatterPlot, baselineEndpointPlot, uniqueLabValues, uniqueVisits, uniqueTreatmentArms, disributionBoxPlot);
 
  }


  private generateUI(dm: DG.DataFrame, lb: DG.DataFrame, grid: DG.Grid,
    hysLawScatterPlot: DG.ScatterPlotViewer,
    baselineEndpointPlot: DG.ScatterPlotViewer,
    labValues: string[] & any,
    visits: string[] & any,
    treatmentArms: string[] & any,
    disributionBoxPlot: DG.Viewer) {


    let baselineEndpointDiv = ui.box(baselineEndpointPlot ? baselineEndpointPlot.root : null);
    let labValue = labValues[ 0 ];
    let bl = visits[ 0 ];
    let ep = visits[ 1 ];
    let labValueChoices = ui.choiceInput('Value', labValue, labValues);
    let blVisitChoices = ui.choiceInput('BL', bl, visits);
    let epVisitChoices = ui.choiceInput('EP', ep, visits);
    labValueChoices.onChanged((v) => {
      labValue = labValueChoices.value;
      this.updateBaselineEndpointPlot(baselineEndpointPlot, baselineEndpointDiv, dm, lb, labValue, bl, ep);
    });
    blVisitChoices.onChanged((v) => {
      bl = blVisitChoices.value;
      this.updateBaselineEndpointPlot(baselineEndpointPlot, baselineEndpointDiv, dm, lb, labValue, bl, ep);
    });
    epVisitChoices.onChanged((v) => {
      ep = epVisitChoices.value;
      this.updateBaselineEndpointPlot(baselineEndpointPlot, baselineEndpointDiv, dm, lb, labValue, bl, ep);
    });

    let distributionDiv = ui.box(disributionBoxPlot ? disributionBoxPlot.root : null);
    let labValBoxPlot = labValues[ 0 ];
    let trArm = treatmentArms[ 0 ];
    let labValueChoicesBoxPlot = ui.choiceInput('Value', labValBoxPlot, labValues);
    let treatmentArmsChoices = ui.choiceInput('Treatment arm', trArm, treatmentArms);
    labValueChoicesBoxPlot.onChanged((v) => {
      labValBoxPlot = labValueChoicesBoxPlot.value;
      this.updateDistributionBoxPlot(disributionBoxPlot, distributionDiv, dm, lb, labValBoxPlot, trArm);
    });
    treatmentArmsChoices.onChanged((v) => {
      trArm = treatmentArmsChoices.value;
      this.updateDistributionBoxPlot(disributionBoxPlot, distributionDiv, dm, lb, labValBoxPlot, trArm);
    });

    this.root.className = 'grok-view ui-box';
   
    this.root.appendChild(
      ui.tabControl({
        "Hy's law":hysLawScatterPlot.root,
        "Baseline endpoint":ui.splitV([
          ui.box(ui.panel([
            ui.divH([
              labValueChoices.root,blVisitChoices.root,epVisitChoices.root
            ])
          ]),{style:{maxHeight:'80px'}}),
          baselineEndpointDiv
        ]),
        "Laboratory distribution":ui.splitV([
          ui.box(ui.panel([
            ui.divH([
              labValueChoicesBoxPlot.root ,treatmentArmsChoices.root
            ])
          ]),{style:{maxHeight:'80px'}}),
          distributionDiv
        ]),
        "Results":grid.root
      }).root
    );
  
  }


  private updateBaselineEndpointPlot(baselineEndpointPlot: DG.ScatterPlotViewer, baselineEndpointDiv: HTMLElement,
    dm: DG.DataFrame, lb: DG.DataFrame, labValue: string, bl: string, ep: string) {
    baselineEndpointDiv.innerHTML = '';
    baselineEndpointPlot = this.baselineEndpointPlot(dm, lb, labValue, bl, ep);
    baselineEndpointDiv.append(baselineEndpointPlot.root);
  }


  private updateDistributionBoxPlot(disributionBoxPlot: any, distributionDiv: HTMLElement,
    dm: DG.DataFrame, lb: DG.DataFrame, labValue: string, trArm: string) {
    distributionDiv.innerHTML = '';
    disributionBoxPlot = this.labValuesDistributionPlot(dm, lb, labValue, trArm);
    distributionDiv.append(disributionBoxPlot.root);
  }


  private hysLawScatterPlot(dm: DG.DataFrame, lb: DG.DataFrame) {
    let hysLawDataframe = createHysLawDataframe(lb, dm);
    grok.data.linkTables(lb, hysLawDataframe,
      [SUBJECT_ID], [SUBJECT_ID],
      [DG.SYNC_TYPE.CURRENT_ROW_TO_ROW, DG.SYNC_TYPE.CURRENT_ROW_TO_SELECTION]);
    grok.data.linkTables(hysLawDataframe, lb,
      [SUBJECT_ID], [SUBJECT_ID],
      [DG.SYNC_TYPE.SELECTION_TO_SELECTION, DG.SYNC_TYPE.SELECTION_TO_SELECTION]);
    let hysLawScatterPlot = createHysLawScatterPlot(hysLawDataframe, ALT, BILIRUBIN, TREATMENT_ARM);
    return hysLawScatterPlot;
  }

  private baselineEndpointPlot(dm: DG.DataFrame, lb: DG.DataFrame, value: any, bl: any, ep: any) {
    let visitCol = LAB_VISIT_NAME;
    let blVisit = bl;
    let epVisit = ep;
    let labValue = value;
    let blNumCol = `${labValue}_BL`;
    let epNumCol = `${labValue}_EP`;
    let baselineEndpointDataframe = createBaselineEndpointDataframe(lb, dm, [TREATMENT_ARM], labValue, blVisit, epVisit, visitCol, blNumCol, epNumCol);
    grok.data.linkTables(lb, baselineEndpointDataframe,
      [ SUBJECT_ID ], [ SUBJECT_ID ],
      [ DG.SYNC_TYPE.CURRENT_ROW_TO_ROW, DG.SYNC_TYPE.CURRENT_ROW_TO_SELECTION ]);
    let baselineEndpointPlot = createBaselineEndpointScatterPlot(baselineEndpointDataframe, blNumCol, epNumCol, TREATMENT_ARM,
      baselineEndpointDataframe.get(LAB_LO_LIM_N, 0), baselineEndpointDataframe.get(LAB_HI_LIM_N, 0));
    return baselineEndpointPlot;
  }

  private labValuesDistributionPlot(dm: DG.DataFrame, lb: DG.DataFrame, selectedlabValue: any, trArm: any) {
    let labValue = selectedlabValue;
    let labValueNumColumn = `${labValue} values`;
    let disributionDataframe = createLabValuesByVisitDataframe(lb, dm, labValue, trArm, labValueNumColumn, LAB_VISIT_DAY);
    let disributionBoxPlot = DG.Viewer.boxPlot(disributionDataframe, {
      category: LAB_VISIT_DAY,
      value: labValueNumColumn,
    });
    return disributionBoxPlot;
  }

}
