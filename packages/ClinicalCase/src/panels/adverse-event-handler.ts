import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { study, ClinRow } from "../clinical-study";
import { AE_DECOD_TERM, AE_START_DATE, SUBJECT_ID } from '../constants/columns-constants';

export class AdverseEventHandler extends DG.ObjectHandler {
  get type() { return 'Adverse Event' }

  isApplicable(x) {
    return (x instanceof ClinRow) && x.row.table.name.toLowerCase() == 'ae';
  }

  renderProperties(x: ClinRow) {
    let day = x.row[AE_START_DATE];

    study.domains.cm.rows
      .match({
        USUBJID: x.row[SUBJECT_ID],
        CMSTDY: `${day - 10}-${day}`,
      })
      .filter();

    return ui.divV([
      ui.h2(x.row[SUBJECT_ID]),
      ui.h2(x.row[AE_DECOD_TERM]),
      study.domains.cm.plot.grid().root
    ]);
  }
}