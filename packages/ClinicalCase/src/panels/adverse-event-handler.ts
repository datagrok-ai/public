import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { study, ClinRow } from "../clinical-study";

export class AdverseEventHandler extends DG.ObjectHandler {
  get type() { return 'Adverse Event' }

  isApplicable(x) {
    return (x instanceof ClinRow) && x.row.table.name.toLowerCase() == 'ae';
  }

  renderProperties(x: ClinRow) {
    let day = x.row['AESTDY'];

    study.domains.cm.rows
      .match({
        USUBJID: x.row['USUBJID'],
        CMSTDY: `${day - 10}-${day}`,
      })
      .filter();

    return ui.divV([
      ui.h2(x.row['USUBJID']),
      ui.h2(x.row['AEDECOD']),
      study.domains.cm.plot.grid().root
    ]);
  }
}