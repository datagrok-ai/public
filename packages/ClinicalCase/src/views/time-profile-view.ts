import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {addDataFromDmDomain, createPivotedDataframeAvg,
  getUniqueValues, getVisitNamesAndDays} from '../data-preparation/utils';
import {ETHNIC, LAB_RES_N, LAB_TEST, VISIT_DAY, RACE, SEX,
  SUBJECT_ID, VS_TEST, VS_RES_N, VISIT_DAY_STR, BW_TEST, BW_RES_N,
  BG_TEST, BG_RES_N,
  VISIT} from '../constants/columns-constants';
import {createVisitDayStrCol, dynamicComparedToBaseline} from '../data-preparation/data-preparation';
import {updateDivInnerHTML} from '../utils/utils';
import {_package} from '../package';
import {ClinicalCaseViewBase} from '../model/ClinicalCaseViewBase';
import {TRT_ARM_FIELD} from '../views-config';
import {TIME_PROFILE_VIEW_NAME} from '../constants/view-names-constants';
import {studies} from '../package';
import {CDISC_STANDARD} from '../utils/types';


export class TimeProfileView extends ClinicalCaseViewBase {
  blVisitChoices: DG.InputBase;
  epVisitChoices: DG.InputBase;
  labChoices: DG.InputBase;
  blVisitDiv = ui.div();
  epVisitDiv = ui.div();
  labChoicesDiv = ui.div();
  laboratoryDataFrame: DG.DataFrame;
  relativeChangeFromBlDataFrame: DG.DataFrame;
  uniqueLabValues: any;
  uniqueVisits: any;
  splitBy: any;
  types = ['Values', 'Changes'];
  domains = ['vs', 'lb', 'bw', 'bg'];
  domainFields = {'lb': {'test': LAB_TEST, 'res': LAB_RES_N}, 'vs': {'test': VS_TEST, 'res': VS_RES_N},
    'bw': {'test': BW_TEST, 'res': BW_RES_N}, 'bg': {'test': BG_TEST, 'res': BG_RES_N}};
  visitDayColumnsDict: {[key: string]: string} = {'lb': VISIT_DAY, 'vs': VISIT_DAY, 'bw': VISIT_DAY, 'bg': VISIT_DAY};
  selectedLabValue: string;
  selectedType: string;
  bl: string;
  ep: string;
  visitNamesAndDays: any [];
  linechart: any;
  selectedDomain: string;
  isSend = false;

  constructor(name, studyId) {
    super(name, studyId);
    this.name = name;
    this.helpUrl = `${_package.webRoot}/views_help/time_profile.md`;
  }

  createView(): void {
    this.isSend = studies[this.studyId].config.standard === CDISC_STANDARD.SEND;
    this.splitBy = [studies[this.studyId].viewsConfig.config[TIME_PROFILE_VIEW_NAME][TRT_ARM_FIELD], SEX, RACE, ETHNIC]
      .filter((it) => studies[this.studyId].domains.dm &&
      studies[this.studyId].domains.dm.columns.names().includes(it));
    if (this.isSend)
      this.domains.forEach((it) => createVisitDayStrCol(studies[this.studyId].domains[it], this.visitDayColumnsDict));
    this.domains = this.domains.filter((it) =>
      studies[this.studyId].domains[it] !== null && !this.optDomainsWithMissingCols.includes(it) &&
        Object.keys(this.visitDayColumnsDict).includes(it));
    this.selectedDomain = this.domains[0];
    this.uniqueLabValues = Array.from(getUniqueValues(studies[this.studyId].domains[this.selectedDomain],
      this.domainFields[this.selectedDomain]['test']));
    this.uniqueVisits = Array.from(getUniqueValues(studies[this.studyId].domains[this.selectedDomain],
      this.isSend ? VISIT_DAY_STR : VISIT));
    this.selectedLabValue = this.uniqueLabValues[0] as string;
    this.selectedType = this.types[0];
    this.visitNamesAndDays = getVisitNamesAndDays(studies[this.studyId].domains[this.selectedDomain],
      this.isSend ? VISIT_DAY_STR : VISIT, this.visitDayColumnsDict[this.selectedDomain]);
    this.bl = this.visitNamesAndDays[0].name;
    this.ep = this.visitNamesAndDays[this.visitNamesAndDays.length-1].name;
    this.createLaboratoryDataframe();

    const domainChoices = ui.input.choice('', {value: this.selectedDomain, items: this.domains});
    domainChoices.onChanged.subscribe((value) => {
      this.selectedDomain = value;
      this.uniqueLabValues = Array.from(getUniqueValues(studies[this.studyId].domains[this.selectedDomain],
        this.domainFields[this.selectedDomain]['test']));
      this.uniqueVisits = Array.from(getUniqueValues(studies[this.studyId].domains[this.selectedDomain],
        this.isSend ? VISIT_DAY_STR : VISIT));
      this.selectedLabValue = this.uniqueLabValues[0] as string;
      this.visitNamesAndDays = getVisitNamesAndDays(studies[this.studyId].domains[this.selectedDomain],
        this.isSend ? VISIT_DAY_STR : VISIT,
        this.visitDayColumnsDict[this.selectedDomain]);
      if (this.visitNamesAndDays.findIndex((it) => it.name === this.bl) === -1)
        this.bl = this.visitNamesAndDays[0].name;

      if (this.visitNamesAndDays.findIndex((it) => it.name === this.ep) === -1)
        this.ep = this.visitNamesAndDays[this.visitNamesAndDays.length-1].name;

      this.createDropdownLists();
      this.updateTimeProfile();
    });

    const typeChoices = ui.input.choice('', {value: this.selectedType, items: this.types});
    typeChoices.onChanged.subscribe((value) => {
      this.selectedType = value;
      this.updateTimeProfile();
    });

    this.createDropdownLists();

    this.root.className = 'grok-view ui-box';
    this.linechart = DG.Viewer.lineChart(this.laboratoryDataFrame, {
      splitColumnName: this.splitBy.length ? this.splitBy[0] : '',
      xColumnName: this.visitDayColumnsDict[this.selectedDomain],
      yColumnNames: [`${this.selectedLabValue} avg(${this.domainFields[this.selectedDomain]['res']})`],
      whiskersType: 'Med | Q1, Q3',
    });
    // if (studies[this.studyId].domains.dm) {
    //   grok.data.linkTables(studies[this.studyId].domains.dm, this.linechart.dataFrame,
    //     [SUBJECT_ID], [SUBJECT_ID],
    //     [DG.SYNC_TYPE.FILTER_TO_FILTER]);
    // }

    this.root.append(this.linechart.root);
    this.setRibbonPanels([
      [
        ui.span(['Plot ']),
        domainChoices.root,
        this.labChoicesDiv,
        typeChoices.root,
        ui.span([' from ']),
        this.blVisitDiv,
        ui.span([' to ']),
        this.epVisitDiv,
      ],
    ]);
  }

  private updateTimeProfile() {
    switch (this.selectedType) {
    case 'Values': {
      this.createLaboratoryDataframe();
      this.linechart.dataFrame = this.laboratoryDataFrame;
      this.linechart.setOptions({
        yColumnNames: [`${this.selectedLabValue} avg(${this.domainFields[this.selectedDomain]['res']})`],
      });
      break;
    }
    case 'Changes': {
      this.createrelativeChangeFromBlDataframe();
      this.linechart.dataFrame = this.relativeChangeFromBlDataFrame;
      this.linechart.setOptions({
        yColumnNames: [`${this.selectedLabValue} avg(${this.domainFields[this.selectedDomain]['res']})`],
      });
      break;
    }
    default: {
      break;
    }
    }
    // if (studies[this.studyId].domains.dm) {
    //   grok.data.linkTables(studies[this.studyId].domains.dm, this.linechart.dataFrame,
    //     [SUBJECT_ID], [SUBJECT_ID],
    //     [DG.SYNC_TYPE.FILTER_TO_FILTER]);
    // }
  }

  private createDropdownLists() {
    this.createValuesChoices();
    this.createBlChoices();
    this.createEpChoices();
  }

  private createValuesChoices() {
    this.labChoices = ui.input.choice('', {value: this.selectedLabValue, items: this.uniqueLabValues});
    this.labChoices.onChanged.subscribe((value) => {
      this.selectedLabValue = value;
      this.updateTimeProfile();
    });
    //@ts-ignore
    this.labChoices.input.style.width = '150px';
    updateDivInnerHTML(this.labChoicesDiv, this.labChoices.root);
  }

  private createBlChoices() {
    this.blVisitChoices = ui.input.choice('', {value: this.bl, items: this.uniqueVisits});
    this.blVisitChoices.onChanged.subscribe((value) => {
      this.bl = value;
      this.updateTimeProfile();
    });
    //@ts-ignore
    this.blVisitChoices.input.style.width = '150px';
    updateDivInnerHTML(this.blVisitDiv, this.blVisitChoices.root);
  }

  private createEpChoices() {
    this.epVisitChoices = ui.input.choice('', {value: this.ep, items: this.uniqueVisits});
    this.epVisitChoices.onChanged.subscribe((value) => {
      this.ep = value;
      this.updateTimeProfile();
    });
    //@ts-ignore
    this.epVisitChoices.input.style.width = '150px';
    updateDivInnerHTML(this.epVisitDiv, this.epVisitChoices.root);
  }

  private createLaboratoryDataframe() {
    let df = this.filterDataFrameByDays(studies[this.studyId].domains[this.selectedDomain].clone());
    if (this.splitBy.length) {
      df = addDataFromDmDomain(df, studies[this.studyId].domains.dm,
        [SUBJECT_ID, this.visitDayColumnsDict[this.selectedDomain],
          this.isSend ? VISIT_DAY_STR : VISIT]
          .concat(Object.values(this.domainFields[this.selectedDomain])),
        this.splitBy);
    }

    this.laboratoryDataFrame = createPivotedDataframeAvg(df, [SUBJECT_ID,
      this.visitDayColumnsDict[this.selectedDomain]],
    this.domainFields[this.selectedDomain]['test'], this.domainFields[this.selectedDomain]['res'], this.splitBy);
  }

  private createrelativeChangeFromBlDataframe() {
    let df = this.filterDataFrameByDays(studies[this.studyId].domains[this.selectedDomain].clone());
    dynamicComparedToBaseline(df, this.domainFields[this.selectedDomain]['test'],
      this.domainFields[this.selectedDomain]['res'], this.bl,
      this.isSend ? VISIT_DAY_STR : VISIT, 'LAB_DYNAMIC_BL', true);
    if (this.splitBy.length) {
      df = addDataFromDmDomain(df, studies[this.studyId].domains.dm, [SUBJECT_ID,
        this.visitDayColumnsDict[this.selectedDomain],
        this.isSend ? VISIT_DAY_STR : VISIT,
        this.domainFields[this.selectedDomain]['test'], this.domainFields[this.selectedDomain]['res']], this.splitBy);
    }

    this.relativeChangeFromBlDataFrame = createPivotedDataframeAvg(df, [SUBJECT_ID,
      this.visitDayColumnsDict[this.selectedDomain]],
    this.domainFields[this.selectedDomain]['test'], this.domainFields[this.selectedDomain]['res'], this.splitBy);
  }


  private filterDataFrameByDays(df: DG.DataFrame) {
    const blDay = this.visitNamesAndDays.find((it) => it.name === this.bl).day;
    const epDay = this.visitNamesAndDays.find((it) => it.name === this.ep).day;
    const filteredDf = df.groupBy(df.columns.names())
      // eslint-disable-next-line max-len
      .where(`${this.visitDayColumnsDict[df.name]} >= ${blDay} and ${this.visitDayColumnsDict[df.name]} <= ${epDay} and ${this.domainFields[this.selectedDomain]['test']} = ${this.selectedLabValue}`)
      .aggregate();
    return filteredDf;
  }
}
