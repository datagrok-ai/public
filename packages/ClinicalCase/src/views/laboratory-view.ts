import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { ClinRow, study } from "../clinical-study";
import { createBaselineEndpointDataframe, createHysLawDataframe, createLabValuesByVisitDataframe } from '../data-preparation/data-preparation';
import { ALT, AP, BILIRUBIN, TREATMENT_ARM } from '../constants';
import { createBaselineEndpointScatterPlot, createHysLawScatterPlot } from '../custom-scatter-plots/custom-scatter-plots';
import { getUniqueValues } from '../data-preparation/utils';

export class LaboratoryView extends DG.ViewBase {

  constructor() {
    super();
    let lb = study.domains.lb;
    let dm = study.domains.dm;

    let grid = lb.plot.grid();
    lb.onCurrentRowChanged.subscribe((_) => {
      grok.shell.o = new ClinRow(lb.currentRow);
    });

    let hysLawScatterPlot = this.hysLawScatterPlot(dm, lb);

    let uniqueLabValues = Array.from(getUniqueValues(lb, 'LBTEST'));
    let uniqueVisits = Array.from(getUniqueValues(lb, 'VISIT'));
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


    let baselineEndpointDiv = ui.div(baselineEndpointPlot ? baselineEndpointPlot.root : null);
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

    let distributionDiv = ui.div(disributionBoxPlot ? disributionBoxPlot.root : null);
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

    this.root.appendChild(ui.div([
      ui.divH([
        ui.block50([ ui.h2('Hy\'s Law Chart'), hysLawScatterPlot ? hysLawScatterPlot.root : null ]),
        ui.block50([
          ui.divV([ ui.h2('Baseline/Endpoint with LLN/ULN'),
          labValueChoices.root,
          blVisitChoices.root,
          epVisitChoices.root
          ]),
          baselineEndpointDiv ]),
      ], { style: { width: '100%' } }),
      ui.divH([
        ui.block50([ ui.h2('Laboratory results'), grid.root ]),
        ui.block50([
          ui.divV([ ui.h2('Laboratory results distribution'),
          labValueChoicesBoxPlot.root,
          treatmentArmsChoices.root
          ]),
          distributionDiv ])
      ], { style: { width: '100%' } })
    ]));

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
    if (dm.columns.contains(TREATMENT_ARM) &&
      lb.columns.contains('LBTEST') &&
      lb.columns.contains('USUBJID') &&
      lb.columns.contains('LBORRES') &&
      lb.columns.contains('LBORNRHI')) {
      let hysLawDataframe = createHysLawDataframe(lb, dm);
      grok.data.linkTables(lb, hysLawDataframe,
        [ `USUBJID` ], [ `USUBJID` ],
        [ DG.SYNC_TYPE.CURRENT_ROW_TO_ROW, DG.SYNC_TYPE.CURRENT_ROW_TO_SELECTION ]);
      let hysLawScatterPlot = createHysLawScatterPlot(hysLawDataframe, ALT, BILIRUBIN, TREATMENT_ARM, AP);
      return hysLawScatterPlot;
    }
    return null;
  }

  private baselineEndpointPlot(dm: DG.DataFrame, lb: DG.DataFrame, value: any, bl: any, ep: any) {
    if (dm.columns.contains(TREATMENT_ARM) &&
      lb.columns.contains('LBTEST') &&
      lb.columns.contains('USUBJID') &&
      lb.columns.contains('LBORRES') &&
      lb.columns.contains('VISIT') &&
      lb.columns.contains('LBORNRLO') &&
      lb.columns.contains('LBORNRHI')) {
      let visitCol = 'VISIT';
      let blVisit = bl;
      let epVisit = ep;
      let labValue = value;
      let blNumCol = `${labValue}_BL`;
      let epNumCol = `${labValue}_EP`;
      let baselineEndpointDataframe = createBaselineEndpointDataframe(lb, dm, labValue, blVisit, epVisit, visitCol, blNumCol, epNumCol);
      grok.data.linkTables(lb, baselineEndpointDataframe,
        [ `USUBJID` ], [ `USUBJID` ],
        [ DG.SYNC_TYPE.CURRENT_ROW_TO_ROW, DG.SYNC_TYPE.CURRENT_ROW_TO_SELECTION ]);
      let baselineEndpointPlot = createBaselineEndpointScatterPlot(baselineEndpointDataframe, blNumCol, epNumCol, TREATMENT_ARM,
        parseFloat(baselineEndpointDataframe.get('LBORNRLO', 0)), parseFloat(baselineEndpointDataframe.get('LBORNRHI', 0)));
      return baselineEndpointPlot;
    }
    return null;
  }

  private labValuesDistributionPlot(dm: DG.DataFrame, lb: DG.DataFrame, selectedlabValue: any, trArm: any) {

    if (lb.columns.contains('LBTEST') &&
      lb.columns.contains('USUBJID') &&
      lb.columns.contains('LBORRES') &&
      lb.columns.contains('VISITDY') &&
      dm.columns.contains(TREATMENT_ARM)) {
      let labValue = selectedlabValue;
      let labValueNumColumn = `${labValue} values`;
      let disributionDataframe = createLabValuesByVisitDataframe(lb, dm, labValue, trArm, labValueNumColumn, 'VISITDY');
      let disributionBoxPlot = DG.Viewer.boxPlot(disributionDataframe, {
        category: 'VISITDY',
        value: labValueNumColumn,
      });
      return disributionBoxPlot;
    }
    return null;
  }

}
