import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";
import { addDataFromDmDomain, getUniqueValues } from '../data-preparation/utils';
import { createBaselineEndpointDataframe } from '../data-preparation/data-preparation';
import { ETHNIC, LAB_RES_N, LAB_TEST, VISIT_DAY, VISIT_NAME, RACE, SEX, SUBJECT_ID, TREATMENT_ARM, LAB_LO_LIM_N, LAB_HI_LIM_N, VS_TEST, VS_RES_N } from '../columns-constants';
import { updateDivInnerHTML } from './utils';
import { _package } from '../package';
import { ClinicalCaseViewBase } from '../model/ClinicalCaseViewBase';
var { jStat } = require('jstat')


export class BoxPlotsView extends ClinicalCaseViewBase {

  domains = ['lb', 'vs'];
  domainFields = {'lb': {'test': LAB_TEST, 'res': LAB_RES_N}, 'vs': {'test': VS_TEST, 'res': VS_RES_N}};
  distrDataframe: DG.DataFrame;

  boxPlotDiv = ui.div();
  boxPlots = [];

  uniqueValues = {};
  selectedValuesByDomain = {};
  uniqueVisits: any;
  splitBy =  [TREATMENT_ARM, SEX, RACE, ETHNIC];

  bl = '';
  selectedSplitBy = TREATMENT_ARM;

  distrWithDmData: DG.DataFrame;

  pValuesArray: any;


  constructor(name) {
    super({});
    this.name = name;
    this.helpUrl = `${_package.webRoot}/views_help/biomarkers_distribution.md`;
  }

  createView(): void {

    this.domains = this.domains.filter(it => study.domains[it] !== null);
    this.splitBy = this.splitBy.filter(it => study.domains.dm.columns.names().includes(it));
    this.selectedSplitBy = this.splitBy[0];
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

    this.domains.forEach(it => {
      let df = study.domains[it].clone(null, [SUBJECT_ID, VISIT_NAME, VISIT_DAY, this.domainFields[it]['test'], this.domainFields[it]['res']]);
      df.getCol(this.domainFields[it]['test']).name = 'test';
      df.getCol(this.domainFields[it]['res']).name = 'res';
      if (!this.distrDataframe) {
        this.distrDataframe = df;
      } else {
        this.distrDataframe.append(df, true);
      }
    });

    this.domains.forEach(it => {
      this.uniqueValues[it] = Array.from(getUniqueValues(study.domains[it], this.domainFields[it]['test']));
    });

    /* let minLabVisit = this.distrDataframe.getCol(VISIT_DAY).stats[ 'min' ];
    let minVisitName = this.distrDataframe
      .groupBy([ VISIT_DAY, VISIT_NAME ])
      .where(`${VISIT_DAY} = ${minLabVisit}`)
      .aggregate()
      .get(VISIT_NAME, 0);
    this.bl = minVisitName; */

    this.uniqueVisits = Array.from(getUniqueValues(this.distrDataframe, VISIT_NAME));
    this.bl = this.uniqueVisits[0];
    this.distrWithDmData = addDataFromDmDomain(this.distrDataframe, study.domains.dm, [ SUBJECT_ID, VISIT_DAY, VISIT_NAME, 'test', 'res' ], this.splitBy);
    this.distrWithDmData = this.distrWithDmData
      .groupBy(this.distrWithDmData.columns.names())
      .where(`${VISIT_NAME} = ${this.bl}`)
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
      let multichoices = {};
      this.domains.forEach(domain => {
        let valuesMultiChoices = ui.multiChoiceInput('', this.selectedValuesByDomain[domain], this.uniqueValues[domain])
        valuesMultiChoices.onChanged((v) => {
          this.selectedValuesByDomain[domain] = valuesMultiChoices.value;
        });
        //@ts-ignore
        valuesMultiChoices.input.style.maxWidth = '100%';
        //@ts-ignore
        valuesMultiChoices.input.style.maxHeight = '100%';
        multichoices[domain] = valuesMultiChoices;
      })

      let acc = ui.accordion();
      this.domains.forEach(domain => {
        acc.addCountPane(`${domain}`, () => multichoices[domain].root, () => this.selectedValuesByDomain[domain].length, false);
        let panel = acc.getPane(`${domain}`);
        //@ts-ignore
        $(panel.root).css('display', 'flex');
        //@ts-ignore
        $(panel.root).css('opacity', '1');
      })

      ui.dialog({ title: 'Select values' })
        .add(ui.div(acc.root, { style: { width: '400px', height: '300px' } }))
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
    if (Object.keys(this.selectedValuesByDomain).length && this.bl) {
      this.boxPlots = [];
      this.pValuesArray = [];
      Object.keys(this.selectedValuesByDomain).forEach(domain => {
        this.selectedValuesByDomain[domain].forEach(it => {
          let df = createBaselineEndpointDataframe(this.distrDataframe, study.domains.dm, [category], 'test', 'res', [], it, this.bl, '', VISIT_NAME, `${it}_BL`);
          this.getPValues(df, domain, it, category, `${it}_BL`);
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
            ui.divText(`p-value: ${this.pValuesArray.find(val => val.value === it).pValue.toPrecision(5)}`, viewerTitlePValue)
          ], { style: { maxHeight: '35px' } }));
          const boxPlot = Array.from(getUniqueValues(df, category)).length > 3 ? ui.block([plot.root]) : ui.block50([plot.root]);
          this.boxPlots.push(boxPlot);
        })       
      })
      updateDivInnerHTML(this.boxPlotDiv, ui.block(this.boxPlots));
    }
  }

  private getTopPValues(topNum: number){
    this.pValuesArray = [];
    Object.keys(this.uniqueValues).forEach(domain => {
      this.uniqueValues[domain].forEach(val => this.getPValues(this.distrWithDmData, domain, val, this.selectedSplitBy, 'res'))
    });   
    //@ts-ignore
    this.pValuesArray.sort((a, b) => a.pValue-b.pValue || isNaN(a.pValue)-isNaN(b.pValue));
    for (let i = 0; i < topNum; i++) {
      let domain = this.pValuesArray[i].domain;
      let value = this.pValuesArray[i].value;
      if(!this.selectedValuesByDomain[domain]){
        this.selectedValuesByDomain[domain] = [value];
      } else {
        this.selectedValuesByDomain[domain] = this.selectedValuesByDomain[domain].concat(value);
      }      
    }
  }

  getPValues(df: DG.DataFrame, domain: string, labVal: any, category: any, resColName: string){
      const valueData = df
        .groupBy([ SUBJECT_ID, 'test', resColName, category ])
        .where(`test = ${labVal}`)
        .aggregate();
      const valuesByArm = valueData.groupBy([ category ]).getGroups();
      const dataForAnova = [];
      Object.values(valuesByArm).forEach(it => {
        const labResults = it.getCol(resColName).getRawData();
        dataForAnova.push(Array.from(labResults));
      })
      const pValue = jStat.anovaftest(...dataForAnova);
      this.pValuesArray.push({value: labVal, pValue: pValue, domain: domain});
  }
}