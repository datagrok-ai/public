import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study } from '../clinical-study';
import { createPropertyPanel } from '../panels/panels-service';
import { SUBJECT_ID } from '../constants';

let filters = ['USUBJID', 'AESEV', 'AEBODSYS', 'AESTDY']

export class AeBrowserView extends DG.ViewBase {
  
  domains = ['dm', 'ae', 'cm', 'ex'];
  aeToSelect: DG.DataFrame;
  ae: DG.DataFrame;
  dm: DG.DataFrame;
  cm: DG.DataFrame;
  ex: DG.DataFrame;
  selectedAe: DG.DataFrame;

  constructor(name) {
    super(name);
    this.name = name;
    this.domains.forEach(it => this[it] = study.domains[it].clone());
    this.aeToSelect = study.domains.ae.clone();

     grok.data.linkTables(this.aeToSelect, this.dm,
        [ `USUBJID`], [ `USUBJID` ],
        [ DG.SYNC_TYPE.CURRENT_ROW_TO_ROW, DG.SYNC_TYPE.CURRENT_ROW_TO_FILTER ]); 

   
    this.root.className = 'grok-view ui-box';

    this.root.append(ui.splitH([
        this.getFilters(),
        this.aeToSelect.plot.grid().root
    ]))

    this.subscribeToCurrentRow();

  }

  private getFilters() {
    let chart = DG.Viewer.fromType('Filters', this.ae, {
      'columnNames': filters,
      'showContextMenu': false,
    }).root;
    chart.style.overflowY = 'scroll';
    return chart
  }

  private subscribeToCurrentRow(){
    this.aeToSelect.onCurrentRowChanged.subscribe(() => {
      const currentSubjId = this.aeToSelect.get(SUBJECT_ID, this.aeToSelect.currentRowIdx);
      const currentAeDate = this.aeToSelect.get('AESTDY', this.aeToSelect.currentRowIdx);
      this.domains.forEach(domain => {
          if(domain !== 'dm'){
              this[domain] = study.domains[domain]
              .clone()
              .groupBy(study.domains[domain].columns.names())
              .where(`${SUBJECT_ID} = ${currentSubjId} and ${domain.toUpperCase()}STDY < ${currentAeDate}`)
              .aggregate();
          }
      })
      createPropertyPanel(this);
    })
  }
}