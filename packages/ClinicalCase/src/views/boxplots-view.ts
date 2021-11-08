import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";
import { addDataFromDmDomain, getUniqueValues } from '../data-preparation/utils';
import { createBaselineEndpointDataframe } from '../data-preparation/data-preparation';
import { ETHNIC, LAB_RES_N, LAB_TEST, LAB_VISIT_DAY, LAB_VISIT_NAME, RACE, SEX, SUBJECT_ID, TREATMENT_ARM } from '../columns-constants';
import { checkMissingDomains, updateDivInnerHTML } from './utils';
import { ILazyLoading } from '../lazy-loading/lazy-loading';
import { _package } from '../package';
import { requiredColumnsByView } from '../constants';
var { jStat } = require('jstat')


export class BoxPlotsView extends DG.ViewBase implements ILazyLoading {

  boxPlotDiv = ui.div();
  boxPlots = [];

  uniqueLabValues: any;
  uniqueVisits: any;
  splitBy =  [TREATMENT_ARM, SEX, RACE, ETHNIC];

  selectedLabValues = null;
  bl = '';
  selectedSplitBy = TREATMENT_ARM;

  labWithDmData: DG.DataFrame;

  pValuesArray: any;


  constructor(name) {
    super({});
    this.name = name;
    this.helpUrl = `${_package.webRoot}/views_help/biomarkers_distribution.md`;
  }

  loaded: boolean;

  load(): void {
    checkMissingDomains(requiredColumnsByView[this.name], false, this);
 }

  createView(): void {

    let viewerTitle = {
      style: {
        'color': 'var(--grey-6)',
        'margin': '12px 0px 6px 12px',
        'font-size': '14px',
        'font-weight': 'bold'
      }
    };

    let viewerTitlePValue = {
      style: {
        'color': 'var(--grey-6)',
        'margin': '12px 0px 6px 100px',
        'font-size': '12px',
      }
    };

    let minLabVisit = study.domains.lb.getCol(LAB_VISIT_DAY).stats[ 'min' ];
    let minVisitName = study.domains.lb
      .groupBy([ LAB_VISIT_DAY, LAB_VISIT_NAME ])
      .where(`${LAB_VISIT_DAY} = ${minLabVisit}`)
      .aggregate()
      .get(LAB_VISIT_NAME, 0);
    this.bl = minVisitName;

    this.uniqueVisits = Array.from(getUniqueValues(study.domains.lb, LAB_VISIT_NAME));
    this.labWithDmData = addDataFromDmDomain(study.domains.lb, study.domains.dm, [ SUBJECT_ID, LAB_VISIT_DAY, LAB_VISIT_NAME, LAB_TEST, LAB_RES_N ], [TREATMENT_ARM, SEX, RACE, ETHNIC]);
    this.uniqueLabValues = Array.from(getUniqueValues(this.labWithDmData, LAB_TEST));
    this.labWithDmData = this.labWithDmData
    .groupBy(this.labWithDmData.columns.names())
    .where(`${LAB_VISIT_DAY} = ${minLabVisit}`)
    .aggregate();
    this.getTopPValues(4);
    this.updateBoxPlots(viewerTitle, viewerTitlePValue, this.selectedSplitBy);

    let blVisitChoices = ui.choiceInput('Baseline', this.bl, this.uniqueVisits);
    blVisitChoices.onChanged((v) => {
      this.bl = blVisitChoices.value;
      this.updateBoxPlots(viewerTitle, viewerTitlePValue, this.selectedSplitBy);
    });

    let splitByChoices = ui.choiceInput('Split by', this.selectedSplitBy, this.splitBy);
    splitByChoices.onChanged((v) => {
      this.selectedSplitBy = splitByChoices.value;
      this.updateBoxPlots(viewerTitle, viewerTitlePValue, this.selectedSplitBy);
    });

    let selectBiomarkers = ui.iconFA('cog', () => {
      let labValuesMultiChoices = ui.multiChoiceInput('Select values', this.selectedLabValues, this.uniqueLabValues)
      labValuesMultiChoices.onChanged((v) => {
        this.selectedLabValues = labValuesMultiChoices.value;
      });
      //@ts-ignore
      labValuesMultiChoices.input.style.maxWidth = '100%';
      //@ts-ignore
      labValuesMultiChoices.input.style.maxHeight = '100%';
      ui.dialog({ title: 'Select values' })
        .add(ui.div([ labValuesMultiChoices ], { style: { width: '400px', height: '300px' } }))
        .onOK(() => {
          this.updateBoxPlots(viewerTitle, viewerTitlePValue, this.selectedSplitBy);
        })
        .show();
    });

    this.setRibbonPanels(
      [ [blVisitChoices.root], [splitByChoices.root], [selectBiomarkers] ] ,
 );

    this.root.append(ui.div([
      ui.block([this.boxPlotDiv])
    ]));
    
  }

  private updateBoxPlots(viewerTitle: any, viewerTitlePValue: any, category: string) {
    if (this.selectedLabValues && this.bl) {
      this.boxPlots = [];
      this.pValuesArray = [];
      this.selectedLabValues.forEach(it => {
        let df = createBaselineEndpointDataframe(study.domains.lb, study.domains.dm, [category], it, this.bl, '', LAB_VISIT_NAME, `${it}_BL`);
        this.getPValues(df, it, category, `${it}_BL`);
        const plot = DG.Viewer.boxPlot(df, {
          category: category,
          value: `${it}_BL`,
          labelOrientation: 'Horz',
          markerColor: category,
          showCategorySelector: false,
          showValueSelector: false,
          showPValue: true
        });
        plot.root.prepend(ui.splitH([
          ui.divText(it, viewerTitle), 
          ui.divText(`p-value: ${this.pValuesArray.find(val => val.labValue === it).pValue.toPrecision(5)}`, viewerTitlePValue)
        ], { style: { maxHeight: '35px' } }));
        const boxPlot = Array.from(getUniqueValues(df, category)).length > 3 ? ui.block([plot.root]) : ui.block50([plot.root]);
        this.boxPlots.push(boxPlot);
      })
      updateDivInnerHTML(this.boxPlotDiv, ui.block(this.boxPlots));
    }
  }

  private getTopPValues(topNum: number){
    this.pValuesArray = [];
    this.uniqueLabValues.forEach(val => this.getPValues(this.labWithDmData, val, TREATMENT_ARM, LAB_RES_N));
    this.pValuesArray.sort((a, b) => (a.pValue - b.pValue));
    this.selectedLabValues = this.pValuesArray.slice(0, topNum).map(it => it.labValue);
  }

  getPValues(df: DG.DataFrame, labVal: any, category: any, resColName: string){
      const valueData = df
        .groupBy([ SUBJECT_ID, LAB_TEST, resColName, category ])
        .where(`${LAB_TEST} = ${labVal}`)
        .aggregate();
      const valuesByArm = valueData.groupBy([ category ]).getGroups();
      const dataForAnova = [];
      Object.values(valuesByArm).forEach(it => {
        const labResults = it.getCol(resColName).getRawData();
        dataForAnova.push(Array.from(labResults));
      })
      const pValue = jStat.anovaftest(...dataForAnova);
      this.pValuesArray.push({labValue: labVal, pValue: pValue});
  }
}