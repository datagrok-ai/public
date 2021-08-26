import { ClinicalCaseView } from "../clinical-case-view";
import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";

let links = {
  ae: { key: 'USUBJID', start: 'AESTDY', end: 'AEENDY', event: 'AETERM' },
  cm: { key: 'USUBJID', start: 'CMSTDY', event: 'CMTRT' },
  ex: { key: 'USUBJID', start: 'EXSTDY', event: 'EXTRT' },
  //ex: { key: 'USUBJID', start: 'EXSTDY', end: 'EXENDY', event: 'EXTRT' },
  lb: { key: 'USUBJID', start: 'LBDY', event: 'LBTEST' }
};

export class TimelinesView extends DG.ViewBase {

  constructor() {
    super();

    let result = null;

    let prepare = function (domain: DG.DataFrame) {
      let info = links[domain.name];
      let t = study.domains[domain.name]
        .clone(null, Object.keys(info).map(e => info[e]));
      t.columns.addNew('domain', DG.TYPE.STRING).init(domain.name);
      for (let name in info)
        t.col(info[name]).name = name;
      return t;
    }

    for (let dt of study.domains.all().filter((t) => links[t.name] != null)) {
      let t = prepare(dt);
      if (result == null)
        result = t;
      else
        result.append(t, true);
    }

    result.plot.fromType(DG.VIEWER.TIMELINES).then((v: DG.Viewer) => {
      v.setOptions({
        subjectColumnName: 'key',
        startColumnName: 'start',
        endColumnName: 'end',
        colorByColumnName: 'event',
      });
      this.root.appendChild(v.root);
    });
  }
}