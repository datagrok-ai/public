import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";
import { addDataFromDmDomain, getUniqueValues } from '../data-preparation/utils';
import { createBaselineEndpointDataframe } from '../data-preparation/data-preparation';
import { ETHNIC, RACE, SEX, TREATMENT_ARM } from '../constants';
import { checkDomainExists, updateDivInnerHTML } from './utils';
import { ILazyLoading } from '../lazy-loading/lazy-loading';
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
    this.helpUrl = 'https://raw.githubusercontent.com/datagrok-ai/public/master/packages/ClinicalCase/views_help/biomarkers_distribution.md';
  }

  loaded: boolean;

  load(): void {
    checkDomainExists(['dm', 'lb'], false, this);
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

    let minLabVisit = study.domains.lb.getCol('VISITDY').stats[ 'min' ];
    let minVisitName = study.domains.lb
      .groupBy([ 'VISITDY', 'VISIT' ])
      .where(`VISITDY = ${minLabVisit}`)
      .aggregate()
      .get('VISIT', 0);
    this.bl = minVisitName;

    this.uniqueVisits = Array.from(getUniqueValues(study.domains.lb, 'VISIT'));
    this.labWithDmData = addDataFromDmDomain(study.domains.lb, study.domains.dm, [ 'USUBJID', 'VISITDY', 'VISIT', 'LBTEST', 'LBSTRESN' ], [TREATMENT_ARM, SEX, RACE, ETHNIC]);
    this.uniqueLabValues = Array.from(getUniqueValues(this.labWithDmData, 'LBTEST'));
    this.labWithDmData = this.labWithDmData
    .groupBy(this.labWithDmData.columns.names())
    .where(`VISITDY = ${minLabVisit}`)
    .aggregate();
    this.getTopPValues(4);
    this.updateBoxPlots(viewerTitle, viewerTitlePValue, this.selectedSplitBy);

    let blVisitChoices = ui.choiceInput('Baseleine', this.bl, this.uniqueVisits);
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
        let df = createBaselineEndpointDataframe(study.domains.lb, study.domains.dm, [category], it, this.bl, '', 'VISIT', `${it}_BL`);
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
    this.uniqueLabValues.forEach(val => this.getPValues(this.labWithDmData, val, TREATMENT_ARM, 'LBSTRESN'));
    this.pValuesArray.sort((a, b) => (a.pValue - b.pValue));
    this.selectedLabValues = this.pValuesArray.slice(0, topNum).map(it => it.labValue);
  }

  getPValues(df: DG.DataFrame, labVal: any, category: any, resColName: string){
      const valueData = df
        .groupBy([ 'USUBJID', 'LBTEST', resColName, category ])
        .where(`LBTEST = ${labVal}`)
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