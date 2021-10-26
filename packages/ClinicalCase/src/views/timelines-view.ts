import { ClinicalCaseView } from "../clinical-case-view";
import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";
import { createFilteredTable, dataframeContentToRow } from "../data-preparation/utils";
import { updateDivInnerHTML } from "./utils";
import $ from "cash-dom";
import { createPropertyPanel } from "../panels/panels-service";
import { ILazyLoading } from "../lazy-loading/lazy-loading";

let links = {
  ae: { key: 'USUBJID', start: 'AESTDY', end: 'AEENDY', event: 'AETERM'},
  cm: { key: 'USUBJID', start: 'CMSTDY', end: 'CMENDY', event: 'CMTRT' },
  ex: { key: 'USUBJID', start: 'EXSTDY', end: 'EXENDY', event: 'EXTRT' },
  //ex: { key: 'USUBJID', start: 'EXSTDY', end: 'EXENDY', event: 'EXTRT' },
  //lb: { key: 'USUBJID', start: 'LBDY', event: 'LBTEST' }
};

let filters = {
  ae: {'AE severity': 'AESEV', 'AE body system': 'AEBODSYS' },
  cm: {'Concomitant medication': 'CMTRT'},
  ex: {'Treatment arm': 'EXTRT'}
}

let multichoiceTableOptions = { 'Adverse events': 'ae', 'Concomitant medication intake': 'cm', 'Drug exposure': 'ex' }

export class TimelinesView extends DG.ViewBase implements ILazyLoading {

  options = {
    subjectColumnName: 'key',
    startColumnName: 'start',
    endColumnName: 'end',
    colorByColumnName: 'domain',
    eventColumnName: 'event',
    showEventInTooltip: true
  }

  selectedOptions = [ Object.keys(multichoiceTableOptions)[ 0 ] ];
  selectedDataframes = [ Object.values(multichoiceTableOptions)[ 0 ] ];
  timelinesDiv = ui.box();
  filtersDiv = ui.box();
  resultTables: DG.DataFrame;

  constructor(name) {
    super({});
    this.name = name;
  }

  loaded: boolean;

  load(): void {
        let multiChoiceOptions = ui.multiChoiceInput('', [ this.selectedOptions[ 0 ] ] as any, Object.keys(multichoiceTableOptions))
    multiChoiceOptions.onChanged((v) => {
      this.selectedOptions = multiChoiceOptions.value;
      this.updateSelectedDataframes(this.selectedOptions);
      this.updateTimelinesPlot();
      this.subscribeToSelection();
    });

    let customTitle = {style:{
      'color':'var(--grey-6)',
      'margin-top':'8px',
      'font-size':'16px',
    }};

    let viewerTitle = {style:{
      'color':'var(--grey-6)',
      'margin':'12px 0px 6px 12px',
      'font-size':'16px',
    }};

    this.root.className = 'grok-view ui-box';
    this.append(ui.splitH([
      ui.splitV([
        ui.box(ui.panel([
          ui.divText('Events',customTitle),
          multiChoiceOptions.root
        ]), {style:{maxHeight:'130px'}}),
        ui.box(
          ui.divText('Filters',viewerTitle), {style:{maxHeight:'40px'}}),
        this.filtersDiv
      ], {style:{maxWidth:'250px',}}),
      this.timelinesDiv
    ]))

    this.updateTimelinesPlot();
    this.subscribeToSelection();
  }

  private subscribeToSelection(){
    this.resultTables.onSelectionChanged.subscribe(() => { 
      createPropertyPanel(this);
    })
  }

  private prepare(domain: DG.DataFrame){
    let info = links[ domain.name ];
    let df = study.domains[ domain.name ];
    let t = df.clone(null, Object.keys(info).map(e => info[ e ]));
    let filterCols = filters[domain.name]
    Object.keys(filterCols).forEach(key => {t.columns.addNewString(key).init((i) => df.get(filterCols[key], i));})
    t.columns.addNew('domain', DG.TYPE.STRING).init(domain.name);
    for (let name in info)
      t.col(info[ name ]).name = name;
    // t = createFilteredTable(t, t.columns.names(), 'key = 01-701-1015');
    return t;
  }

  private updateTimelinesTables(){
    this.resultTables = null;
    for (let dt of study.domains.all().filter((t) => this.selectedDataframes.includes(t.name))) {
      let t = this.prepare(dt);
      if (this.resultTables == null)
        this.resultTables = t;
      else
        this.resultTables.append(t, true);
    }
  }

  private updateTimelinesPlot(){
    this.updateTimelinesTables();
    if(this.resultTables){
      this.resultTables.plot.fromType(DG.VIEWER.TIMELINES, {
        paramOptions: JSON.stringify(this.options),
      }).then((v: any) => {
        v.setOptions({
          subjectColumnName: 'key',
          startColumnName: 'start',
          endColumnName: 'end',
          colorByColumnName: 'domain',
          eventColumnName: 'event',
          showEventInTooltip: true
        });
        $(v.root).css('position', 'relative')
        v.zoomState = [[0, 10], [0, 10], [90, 100], [90, 100]];
        v.render();
        this.updateTimelinesDivs(v.root, this.getFilters());
      });
    } else {
      this.updateTimelinesDivs('', '');
    }
    
  }

  private updateSelectedDataframes(options: string[]){
    this.selectedDataframes = [];
    options.forEach(item => {
      this.selectedDataframes.push(multichoiceTableOptions[item])
    })
  }

  private getFilters() {
    let filterColumns = [];
    Object.keys(filters).forEach(domain => filterColumns = filterColumns.concat(Object.keys(filters[ domain ])))
    let chart = DG.Viewer.fromType('Filters', this.resultTables, {
      'columnNames': filterColumns,
      'showContextMenu': false,
    }).root;
    chart.style.overflowY = 'scroll';
    return chart
  }

  private updateTimelinesDivs(timelinesContent: any, filtersContent: any) {
    updateDivInnerHTML(this.timelinesDiv, timelinesContent);
    updateDivInnerHTML(this.filtersDiv, filtersContent);
  }
}