import { ClinicalCaseView } from "../clinical-case-view";
import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";

export class StudySummaryView extends DG.ViewBase {

  constructor() {
    super();

    let subjsPerDay = study.domains.dm.groupBy(['RFSTDTC']).uniqueCount('USUBJID').aggregate();
    let refStartCol = subjsPerDay.col('RFSTDTC');
    for (let i = 0, rowCount = subjsPerDay.rowCount; i < rowCount; i++) {
      if (refStartCol.isNone(i)) subjsPerDay.rows.removeAt(i);
    }

    let lc = DG.Viewer.lineChart(subjsPerDay);
    let bc = DG.Viewer.barChart(study.domains.dm, { split: 'SEX', stack: 'ARMCD' });
    let bp = DG.Viewer.boxPlot(study.domains.dm, {
      labelOrientation: 'Horz',
      showStatistics: true,
      value: 'AGE',
      category: 'ARM'
    });

    this.root.appendChild(ui.div([
      ui.h1('Study Summary'),
      ui.block([ui.divText('Patient Enrollment'), lc.root]),
      ui.divH([
        ui.block25([ui.h2('Sex Distribution'), bc.root]),
        ui.block75([ui.h2('Age Statistics'), bp.root]),
      ], { style: { width: '100%' } }),
    ]));
  }
}