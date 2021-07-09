import { ClinicalCaseView } from "../clinical-case-view";
import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";

export class StudySummaryView extends DG.ViewBase {

  constructor() {
    super();

    this.name = study.name;

    let subjsPerDay = study.domains.dm.groupBy(['RFSTDTC']).uniqueCount('USUBJID').aggregate();
    let refStartCol = subjsPerDay.col('RFSTDTC');
    for (let i = 0, rowCount = subjsPerDay.rowCount; i < rowCount; i++) {
      if (refStartCol.isNone(i))
        subjsPerDay.rows.removeAt(i);
    }
    let lc = DG.Viewer.lineChart(subjsPerDay);

    let summary = ui.tableFromMap({
      'subjects': study.subjectsCount,
      'sites': study.sitesCount
    });

    this.root.appendChild(ui.div([
      ui.divH([
        ui.block25([ui.h2('Summary'), summary]),
        ui.block75([ui.h2('Enrollment'), lc.root])
      ]),
      ui.divH([
        ui.block25([ui.h2('ARM'), DG.Viewer.barChart(study.domains.dm, { split: 'arm', style: 'dashboard' }).root]),
        ui.block25([ui.h2('Sex'), DG.Viewer.barChart(study.domains.dm, { split: 'sex', style: 'dashboard' }).root]),
        ui.block25([ui.h2('Race'), DG.Viewer.barChart(study.domains.dm, { split: 'race', style: 'dashboard' }).root]),
        ui.block25([ui.h2('Age'), DG.Viewer.histogram(study.domains.dm, { value: 'age', style: 'dashboard' }).root]),
      ], { style: { width: '100%' } }),
    ]));
  }
}