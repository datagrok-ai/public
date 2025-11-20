import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {addDataFromDmDomain, getUniqueValues} from '../data-preparation/utils';
import {createBaselineEndpointDataframe, createVisitDayStrCol} from '../data-preparation/data-preparation';
import {ETHNIC, LAB_RES_N, LAB_TEST, VISIT_DAY, RACE, SEX, SUBJECT_ID, VS_TEST, VS_RES_N,
  VISIT_DAY_STR, BW_TEST, BW_RES_N, BG_TEST, BG_RES_N} from '../constants/columns-constants';
import {updateDivInnerHTML} from '../utils/utils';
import {_package} from '../package';
import {ClinicalCaseViewBase} from '../model/ClinicalCaseViewBase';
import {TRT_ARM_FIELD, VISIT_FIELD} from '../views-config';
import {DISTRIBUTIONS_VIEW_NAME} from '../constants/view-names-constants';
import {tTest} from '@datagrok-libraries/statistics/src/tests';
import {studies} from '../package';
const {jStat} = require('jstat');


export class BoxPlotsView extends ClinicalCaseViewBase {
  domains = ['lb', 'vs', 'bw', 'bg'];
  domainFields = {'lb': {'test': LAB_TEST, 'res': LAB_RES_N}, 'vs': {'test': VS_TEST, 'res': VS_RES_N},
    'bw': {'test': BW_TEST, 'res': BW_RES_N}, 'bg': {'test': BG_TEST, 'res': BG_RES_N}};
  distrDataframe: DG.DataFrame;

  boxPlotDiv = ui.div();
  boxPlots = [];

  uniqueValues = {};
  selectedValuesByDomain = {};
  uniqueVisits: any;
  splitBy: any;

  bl = '';
  selectedSplitBy = studies[this.studyId].viewsConfig.config[DISTRIBUTIONS_VIEW_NAME][TRT_ARM_FIELD];
  visitDayColumnsDict: {[key: string]: string} = {'lb': VISIT_DAY, 'vs': VISIT_DAY, 'bw': VISIT_DAY, 'bg': VISIT_DAY};

  distrWithDmData: DG.DataFrame;

  pValuesArray: any;

  viewerTitle = {
    style: {
      'color': 'var(--grey-6)',
      'margin': '12px 0px 6px 12px',
      'font-size': '14px',
      'font-weight': 'bold',
    },
  };

  viewerTitlePValue = {
    style: {
      'color': 'var(--grey-6)',
      'margin': '12px 0px 6px 100px',
      'font-size': '12px',
    },
  };


  constructor(name, studyId) {
    super(name, studyId);
    this.name = name;
    this.helpUrl = `${_package.webRoot}/views_help/biomarkers_distribution.md`;
  }

  createView(): void {
    this.selectedValuesByDomain = {};
    this.domains = this.domains.filter((it) => studies[this.studyId].domains[it] !== null &&
      !this.optDomainsWithMissingCols.includes(it));
    this.splitBy = [studies[this.studyId].viewsConfig.config[DISTRIBUTIONS_VIEW_NAME][TRT_ARM_FIELD], SEX, RACE, ETHNIC]
      .filter((it) => studies[this.studyId].domains.dm.columns.names().includes(it));
    this.selectedSplitBy = this.splitBy[0];

    if (studies[this.studyId].viewsConfig.config[this.name][VISIT_FIELD] === VISIT_DAY_STR)
      this.domains.forEach((it) => createVisitDayStrCol(studies[this.studyId].domains[it], this.visitDayColumnsDict));
    this.domains = this.domains.filter((it) =>
      studies[this.studyId].domains[it] !== null && !this.optDomainsWithMissingCols.includes(it) &&
        Object.keys(this.visitDayColumnsDict).includes(it));

    this.domains.forEach((it) => {
      const df = (studies[this.studyId].domains[it] as DG.DataFrame).clone(null, [SUBJECT_ID,
        studies[this.studyId].viewsConfig.config[this.name][VISIT_FIELD], this.visitDayColumnsDict[it],
        this.domainFields[it]['test'], this.domainFields[it]['res']]);
      df.col(this.visitDayColumnsDict[it]).name = VISIT_DAY;
      df.getCol(this.domainFields[it]['test']).name = 'test';
      df.getCol(this.domainFields[it]['res']).name = 'res';
      if (!this.distrDataframe)
        this.distrDataframe = df;
      else
        this.distrDataframe.append(df, true);
    });

    this.domains.forEach((it) => {
      this.uniqueValues[it] = Array.from(getUniqueValues(studies[this.studyId].domains[it],
        this.domainFields[it]['test']));
    });

    /* let minLabVisit = this.distrDataframe.getCol(VISIT_DAY).stats[ 'min' ];
    let minVisitName = this.distrDataframe
      .groupBy([ VISIT_DAY, VISIT ])
      .where(`${VISIT_DAY} = ${minLabVisit}`)
      .aggregate()
      .get(VISIT, 0);
    this.bl = minVisitName; */

    this.uniqueVisits = Array.from(getUniqueValues(this.distrDataframe,
      studies[this.studyId].viewsConfig.config[this.name][VISIT_FIELD]));
    this.bl = this.uniqueVisits[0];
    this.distrWithDmData = addDataFromDmDomain(this.distrDataframe, studies[this.studyId].domains.dm,
      [SUBJECT_ID, VISIT_DAY, studies[this.studyId].viewsConfig.config[this.name][VISIT_FIELD],
        'test', 'res'], this.splitBy);
    this.distrWithDmData = this.distrWithDmData
      .groupBy(this.distrWithDmData.columns.names())
      .where(`${studies[this.studyId].viewsConfig.config[this.name][VISIT_FIELD]} = ${this.bl}`)
      .aggregate();
    this.getTopPValues(4);

    this.updateBoxPlots(this.viewerTitle, this.viewerTitlePValue, this.selectedSplitBy);

    const blVisitChoices = ui.input.choice('Baseline', {value: this.bl, items: this.uniqueVisits});
    blVisitChoices.onChanged.subscribe((value) => {
      this.bl = value;
      this.updateBoxPlots(this.viewerTitle, this.viewerTitlePValue, this.selectedSplitBy);
    });

    const splitByChoices = ui.input.choice('Split by', {value: this.selectedSplitBy, items: this.splitBy});
    splitByChoices.onChanged.subscribe((value) => {
      this.selectedSplitBy = value;
      this.updateBoxPlots(this.viewerTitle, this.viewerTitlePValue, this.selectedSplitBy);
    });

    const selectBiomarkers = ui.iconFA('cog', () => {
      const multichoices = {};
      this.domains.forEach((domain) => {
        const valuesMultiChoices = ui.input.multiChoice('', {value: this.selectedValuesByDomain[domain],
          items: this.uniqueValues[domain]});
        valuesMultiChoices.onChanged.subscribe((value) => {
          this.selectedValuesByDomain[domain] = value;
        });
        //@ts-ignore
        valuesMultiChoices.input.style.maxWidth = '100%';
        //@ts-ignore
        valuesMultiChoices.input.style.maxHeight = '100%';
        multichoices[domain] = valuesMultiChoices;
      });

      const acc = ui.accordion();
      this.domains.forEach((domain) => {
        acc.addCountPane(`${domain}`, () => multichoices[domain].root,
          () => this.selectedValuesByDomain[domain].length, false);
        const panel = acc.getPane(`${domain}`);
        //@ts-ignore
        $(panel.root).css('display', 'flex');
        //@ts-ignore
        $(panel.root).css('opacity', '1');
      });

      ui.dialog({title: 'Select values'})
        .add(ui.div(acc.root, {style: {width: '400px', height: '300px'}}))
        .onOK(() => {
          this.updateBoxPlots(this.viewerTitle, this.viewerTitlePValue, this.selectedSplitBy);
        })
        .show();
    });

    this.setRibbonPanels(
      [[blVisitChoices.root], [splitByChoices.root], [selectBiomarkers]],
    );

    this.root.append(ui.div([
      ui.block([this.boxPlotDiv]),
    ]));

    grok.data.linkTables(studies[this.studyId].domains.dm, this.distrDataframe,
      [SUBJECT_ID], [SUBJECT_ID],
      [DG.SYNC_TYPE.FILTER_TO_FILTER]);
  }

  updateGlobalFilter(): void {
    this.updateBoxPlots(this.viewerTitle, this.viewerTitlePValue, this.selectedSplitBy);
  }

  private updateBoxPlots(viewerTitle: any, viewerTitlePValue: any, category: string) {
    if (Object.keys(this.selectedValuesByDomain).length && this.bl) {
      this.boxPlots = [];
      this.pValuesArray = [];
      Object.keys(this.selectedValuesByDomain).forEach((domain) => {
        this.selectedValuesByDomain[domain].forEach((it) => {
          const df = createBaselineEndpointDataframe(
            this.distrDataframe.clone(this.distrDataframe.filter), studies[this.studyId].domains.dm, [category],
            'test', 'res', [], it, this.bl, '',
            studies[this.studyId].viewsConfig.config[this.name][VISIT_FIELD], `${it}_BL`);
          this.getPValues(df, domain, it, category, `${it}_BL`);
          const plot = DG.Viewer.boxPlot(df, {
            categoryColumnNames: [category],
            value: `${it}_BL`,
            labelOrientation: DG.TextOrientation.Horz,
            markerColor: category,
            showCategorySelector: false,
            showValueSelector: false,
            showPValue: true,
          });
          plot.root.prepend(ui.splitH([
            ui.divText(it, viewerTitle),
            ui.divText(`p-value: ${this.pValuesArray.find((val) => val.value === it).pValue.toPrecision(5)}`,
              viewerTitlePValue),
          ], {style: {maxHeight: '35px'}}));
          const boxPlot = Array.from(getUniqueValues(df, category)).length > 3 ? ui.block([plot.root]) :
            ui.block50([plot.root]);
          this.boxPlots.push(boxPlot);
        });
      });
      updateDivInnerHTML(this.boxPlotDiv, ui.block(this.boxPlots));
    }
  }

  private getTopPValues(topNum: number) {
    this.pValuesArray = [];
    Object.keys(this.uniqueValues).forEach((domain) => {
      this.uniqueValues[domain].forEach((val) =>
        this.getPValues(this.distrWithDmData, domain, val, this.selectedSplitBy, 'res'));
    });
    //@ts-ignore
    this.pValuesArray.sort((a, b) => a.pValue-b.pValue || isNaN(a.pValue)-isNaN(b.pValue));
    for (let i = 0; i < topNum; i++) {
      const domain = this.pValuesArray[i].domain;
      const value = this.pValuesArray[i].value;
      if (!this.selectedValuesByDomain[domain])
        this.selectedValuesByDomain[domain] = [value];
      else
        this.selectedValuesByDomain[domain] = this.selectedValuesByDomain[domain].concat(value);
    }
  }

  getPValues(df: DG.DataFrame, domain: string, labVal: any, category: any, resColName: string) {
    const valueData = df
      .groupBy([SUBJECT_ID, 'test', resColName, category])
      .where(`test = ${labVal}`)
      .aggregate();
    const valuesByArm = valueData.groupBy([category]).getGroups();
    const dataForAnova = [];
    Object.values(valuesByArm).forEach((it) => {
      const labResults = it.getCol(resColName).getRawData();
      dataForAnova.push(Array.from(labResults));
    });
    const pValue = dataForAnova.length === 2 ? tTest(dataForAnova[0],
      dataForAnova[1])['p-value'] : jStat.anovaftest(...dataForAnova);
    this.pValuesArray.push({value: labVal, pValue: pValue, domain: domain});
  }
}
