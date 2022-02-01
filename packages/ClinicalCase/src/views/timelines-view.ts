import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import * as grok from 'datagrok-api/grok';
import { study } from "../clinical-study";
import { updateDivInnerHTML } from "./utils";
import $ from "cash-dom";
import { AEBrowserHelper } from "../helpers/ae-browser-helper";
import { _package } from "../package";
import { AE_BODY_SYSTEM, AE_END_DAY, AE_SEVERITY, AE_START_DAY, AE_TERM, CON_MED_DOSE, CON_MED_DOSE_FREQ, CON_MED_DOSE_UNITS, CON_MED_END_DAY, CON_MED_NAME, CON_MED_ROUTE, CON_MED_START_DAY, DOMAIN, INV_DRUG_DOSE, INV_DRUG_DOSE_FORM, INV_DRUG_DOSE_FREQ, INV_DRUG_DOSE_UNITS, INV_DRUG_END_DAY, INV_DRUG_NAME, INV_DRUG_ROUTE, INV_DRUG_START_DAY, SUBJECT_ID, TREATMENT_ARM } from "../columns-constants";
import { ClinicalCaseViewBase } from "../model/ClinicalCaseViewBase";
import { addDataFromDmDomain, getNullOrValue } from "../data-preparation/utils";

let links = {
  ae: { key: SUBJECT_ID, start: AE_START_DAY, end: AE_END_DAY, event: AE_TERM },
  cm: { key: SUBJECT_ID, start: CON_MED_START_DAY, end: CON_MED_END_DAY, event: CON_MED_NAME },
  ex: { key: SUBJECT_ID, start: INV_DRUG_START_DAY, end: INV_DRUG_END_DAY, event: INV_DRUG_NAME },
  //ex: { key: 'USUBJID', start: 'EXSTDY', end: 'EXENDY', event: 'EXTRT' },
  //lb: { key: 'USUBJID', start: 'LBDY', event: 'LBTEST' }
};

let filters = {
  ae: { 'AE severity': AE_SEVERITY, 'AE body system': AE_BODY_SYSTEM },
  cm: { 'Concomitant medication': CON_MED_NAME },
  ex: { 'Treatment arm': INV_DRUG_NAME }
}

let multichoiceTableDict = { 'Adverse events': 'ae', 'Concomitant medication intake': 'cm', 'Drug exposure': 'ex' }

export class TimelinesView extends ClinicalCaseViewBase {

  options = {
    splitByColumnName: 'key',
    startColumnName: 'start',
    endColumnName: 'end',
    colorByColumnName: 'domain',
    eventColumnName: 'event',
    showEventInTooltip: true
  }

  multichoiceTableOptions: any;
  selectedOptions: any;
  selectedDataframes: any;
  timelinesDiv = ui.box();
  filtersDiv = ui.box();
  resultTables: DG.DataFrame;
  aeBrowserHelper: AEBrowserHelper;
  filterColumns = [];

  constructor(name) {
    super({});
    this.name = name;
    this.helpUrl = `${_package.webRoot}/views_help/timelines.md`;
    //@ts-ignore
    this.basePath = '/timelines';
  }

  createView(): void {
    let existingTables = study.domains.all().map(it => it.name);
    Object.keys(filters).forEach(domain => this.filterColumns = this.filterColumns.concat(Object.keys(filters[domain])));
    this.multichoiceTableOptions = {};
    this.multichoiceTableOptions = Object.fromEntries(Object.entries(multichoiceTableDict).filter(([k, v]) => existingTables.includes(v)));
    this.selectedOptions = [Object.keys(this.multichoiceTableOptions)[0]];
    this.selectedDataframes = [Object.values(this.multichoiceTableOptions)[0]];
    let multiChoiceOptions = ui.multiChoiceInput('', [this.selectedOptions[0]] as any, Object.keys(this.multichoiceTableOptions));
    multiChoiceOptions.onChanged((v) => {
      this.selectedOptions = multiChoiceOptions.value;
      this.updateSelectedDataframes(this.selectedOptions);
      this.updateTimelinesPlot();
      this.subscribeToSelection();
    });

    let customTitle = {
      style: {
        'color': 'var(--grey-6)',
        'margin-top': '8px',
        'font-size': '16px',
      }
    };

    let viewerTitle = {
      style: {
        'color': 'var(--grey-6)',
        'margin': '12px 0px 6px 12px',
        'font-size': '16px',
      }
    };

    this.root.className = 'grok-view ui-box';
    this.append(ui.splitH([
      ui.splitV([
        ui.box(ui.panel([
          ui.divText('Events', customTitle),
          multiChoiceOptions.root
        ]), { style: { maxHeight: '130px' } }),
        ui.box(
          ui.divText('Filters', viewerTitle), { style: { maxHeight: '40px' } }),
        this.filtersDiv
      ], { style: { maxWidth: '250px', } }),
      this.timelinesDiv
    ]))

    this.updateTimelinesPlot();
    this.subscribeToSelection();
  }

  private subscribeToSelection() {
    let setPropPanel = async () => {
      const panel = await this.propertyPanel();
      if (panel) {
        grok.shell.o = panel;
      }
    }
    this.resultTables.onSelectionChanged.subscribe(() => {
      setPropPanel();
    })

  }

  private prepare(domain: DG.DataFrame) {
    let info = links[domain.name];
    let df = study.domains[domain.name];
    let t = df.clone(null, Object.keys(info).map(e => info[e]));
    let filterCols = Object.keys(filters[domain.name]).filter(it => df.columns.names().includes(it));
    filterCols.forEach(key => { t.columns.addNewString(key).init((i) => df.get(filterCols[key], i)); })
    t.columns.addNew('domain', DG.TYPE.STRING).init(domain.name.toLocaleLowerCase());
    t.columns.addNewFloat('rowNum').init((i) => i);
    for (let name in info)
      t.col(info[name]).name = name;
    return t;
  }

  private updateTimelinesTables() {
    this.resultTables = null;
    for (let dt of study.domains.all().filter((t) => this.selectedDataframes.includes(t.name))) {
      let t = this.prepare(dt);
      if (this.resultTables == null)
        this.resultTables = t;
      else
        this.resultTables.append(t, true);
    }
  }

  private updateTimelinesPlot() {
    this.updateTimelinesTables();
    if (this.resultTables) {
      this.resultTables.plot.fromType(DG.VIEWER.TIMELINES, {
        paramOptions: JSON.stringify(this.options),
      }).then((v: any) => {
        v.setOptions({
          splitByColumnName: 'key',
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

  private updateSelectedDataframes(options: string[]) {
    this.selectedDataframes = [];
    options.forEach(item => {
      this.selectedDataframes.push(this.multichoiceTableOptions[item])
    })
  }

  private getFilters() {
    let chart = DG.Viewer.fromType('Filters', this.resultTables, {
      'columnNames': this.filterColumns,
      'showContextMenu': true,
    }).root;
    chart.style.overflowY = 'scroll';
    return chart
  }

  private updateTimelinesDivs(timelinesContent: any, filtersContent: any) {
    updateDivInnerHTML(this.timelinesDiv, timelinesContent);
    updateDivInnerHTML(this.filtersDiv, filtersContent);
  }


  override async propertyPanel() {
    const selectedInd = this.resultTables.selection.getSelectedIndexes();

    if (!selectedInd.length) {
      if (this.selectedDataframes.includes('ae')) {
        const aeWithArm = addDataFromDmDomain(this.resultTables.clone(), study.domains.dm, this.resultTables.columns.names(), [TREATMENT_ARM], 'key');
        const aeNumberByArm = aeWithArm.groupBy([TREATMENT_ARM]).where(`${DOMAIN} = ae`).count().aggregate();
        const subjNumberByArm = study.domains.dm.groupBy([TREATMENT_ARM]).uniqueCount(SUBJECT_ID).aggregate();
        const aeNumByArmDict = {};
        for (let i = 0; i < aeNumberByArm.rowCount; i++) {
          aeNumByArmDict[aeNumberByArm.get(TREATMENT_ARM, i)] = aeNumberByArm.get('count', i) /
            subjNumberByArm.
              groupBy([TREATMENT_ARM, `unique(${SUBJECT_ID})`])
              .where(`${TREATMENT_ARM} = ${aeNumberByArm.get(TREATMENT_ARM, i)}`)
              .aggregate()
              .get(`unique(${SUBJECT_ID})`, 0);
        }

        const aeTop5Dict = {};
        const aeTop5 = aeWithArm.groupBy([TREATMENT_ARM, 'event']).where(`${DOMAIN} = ae`).count().aggregate();
        const order = aeTop5.getSortedOrder([TREATMENT_ARM, 'event'] as any);
        order.forEach(item => {
          const arm = aeTop5.get(TREATMENT_ARM, item);
          if (Object.keys(aeTop5Dict).includes(arm) && aeTop5Dict[arm].length < 5) {
            aeTop5Dict[arm].push(aeTop5.get('event', item));
          } else {
            aeTop5Dict[arm] = [aeTop5.get('event', item)]
          }
        })

        Object.keys(aeTop5Dict).forEach(key => {
          aeTop5Dict[key] = aeTop5Dict[key].join(', ')
        })

        let acc = this.createAccWithTitle(this.name);
        let avarageae = ui.tableFromMap(aeNumByArmDict);
        let mostae = ui.tableFromMap(aeTop5Dict);

        acc.addPane('Average AE per patient', () => {
          $(avarageae).find('tr').css('vertical-align', 'top');
          $(avarageae).find('td').css('padding-bottom', '10px');
          $(avarageae).find('.d4-entity-list>span').css('margin', '0px');
          return avarageae
        }, true);

        acc.addPane('Frequent AEs', () => {
          $(mostae).find('tr').css('vertical-align', 'top');
          $(mostae).find('td').css('padding-bottom', '10px');
          $(mostae).find('.d4-entity-list>span').css('margin', '0px');
          return mostae
        });
        return acc.root;
      }
    } else {
      const eventArray = [];
      const eventIndexesArray = [];

      let domainAdditionalFields = {
        ae: { fields: [AE_SEVERITY], pos: 'before' },
        cm: { fields: [CON_MED_DOSE, CON_MED_DOSE_UNITS, CON_MED_DOSE_FREQ, CON_MED_ROUTE], pos: 'after' },
        ex: { fields: [INV_DRUG_DOSE, INV_DRUG_DOSE_UNITS, INV_DRUG_DOSE_FORM, INV_DRUG_DOSE_FREQ, INV_DRUG_ROUTE], pos: 'after' },
      }

      let switchToAEBrowserPanel = (aeRowNum) => {
        if (this.aeBrowserHelper.aeToSelect.currentRowIdx === aeRowNum) {
          this.aeBrowserHelper.propertyPanel();
          return null;
        } else {
          //@ts-ignore
          this.aeBrowserHelper.aeToSelect.currentRow = aeRowNum;
        }
      }

      if (selectedInd.length === 1 && this.resultTables.get('domain', selectedInd[0]) === 'ae') {
        const aeRowNum = this.resultTables.get('rowNum', selectedInd[0]);
        switchToAEBrowserPanel(aeRowNum);
      } else {
        selectedInd.forEach((item) => {
          const domain = this.resultTables.get('domain', item);
          const addDomainInfo = domainAdditionalFields[domain];
          let addInfoString = '';
          const index = this.resultTables.get('rowNum', item);
          addDomainInfo.fields.forEach(it => {
            if (study.domains[domain].columns.names().includes(it)) {
              addInfoString += `${study.domains[domain].get(it, index)} `;
            }
          });
          const eventName = String(getNullOrValue(this.resultTables, 'event', item)).toLowerCase();
          const fullEventName = addDomainInfo.pos === 'before' ? `${addInfoString}${eventName}` : `${eventName} ${addInfoString}`;
          const eventStart = getNullOrValue(this.resultTables, 'start', item);
          const eventEnd = getNullOrValue(this.resultTables, 'end', item);
          eventArray.push({
            Domain: domain,
            Event: fullEventName,
            Days: `${eventStart} - ${eventEnd}`
          })
          eventIndexesArray.push(index);
        })

        let acc2 = this.createAccWithTitle(`${this.name} patient`, `${this.resultTables.get('key', selectedInd[0])}`);
  
        const eventTable = DG.DataFrame.fromObjects(eventArray);
        const eventGrid = eventTable.plot.grid();
        eventGrid.columns.byName('domain').width = 55;
        let col = eventGrid.columns.byName('event');
        col.width = 170;
        col.cellType = 'html';
  
        eventGrid.onCellPrepare(function (gc) {
          if (gc.isTableCell && eventTable.get('Domain', gc.gridRow) === 'ae' && gc.gridColumn.name === 'Event') {
            let eventElement = ui.link(gc.cell.value, {}, '', { id: `${eventIndexesArray[gc.gridRow]}` });
            eventElement.addEventListener('click', (event) => {
              switchToAEBrowserPanel(parseInt(eventElement.id));
              event.stopPropagation();
            });
            gc.style.element = ui.div(eventElement, { style: { 'white-space': 'nowrap' } });
          } else {
            gc.style.element = ui.divText(gc.cell.value, { style: { 'white-space': 'nowrap' } });
          }
          gc.style.element.style.paddingTop = '7px';
          gc.style.element.style.paddingLeft = '7px';
          ui.tooltip.bind(gc.style.element, gc.cell.value);
        });
  
        acc2.addPane('Events', () => {
          return ui.div(eventGrid.root);
        });
        const accPane = acc2.getPane('Events').root.getElementsByClassName('d4-accordion-pane-content')[0] as HTMLElement;
        accPane.style.margin = '0px';
        accPane.style.paddingLeft = '0px';
  
        return acc2.root;

      }
    }

  }

}