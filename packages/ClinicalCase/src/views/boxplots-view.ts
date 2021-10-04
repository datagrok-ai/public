import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study, ClinRow } from "../clinical-study";
import { addDataFromDmDomain, getUniqueValues } from '../data-preparation/utils';
import { createBaselineEndpointDataframe } from '../data-preparation/data-preparation';
import { ETHNIC, RACE, SEX, TREATMENT_ARM } from '../constants';
import { updateDivInnerHTML } from './utils';
var { jStat } = require('jstat')


export class BoxPlotsView extends DG.ViewBase {

  boxPlotDiv = ui.div();
  boxPlots = [];

  uniqueLabValues = Array.from(getUniqueValues(study.domains.lb, 'LBTEST'));
  uniqueVisits = Array.from(getUniqueValues(study.domains.lb, 'VISIT'));
  splitBy =  [TREATMENT_ARM, SEX, RACE, ETHNIC];

  selectedLabValues = null;
  bl = '';
  selectedSplitBy = TREATMENT_ARM;

  labWithDmData: DG.DataFrame;


  constructor(name) {
    super(name);
    this.name = name;
    this.labWithDmData = addDataFromDmDomain(study.domains.lb, study.domains.dm, [ 'USUBJID', 'VISITDY', 'VISIT', 'LBTEST', 'LBORRES' ], [TREATMENT_ARM, SEX, RACE, ETHNIC]);

    let viewerTitle = {
      style: {
        'color': 'var(--grey-6)',
        'margin': '12px 0px 6px 12px',
        'font-size': '16px',
      }
    };

/*     const plot = DG.Viewer.boxPlot(study.domains.dm, {
      category: 'ARM',
      value: 'AGE',
      labelOrientation: 'Horz'
    });
    plot.root.prepend(ui.divText('Age distribution by treatment arm', viewerTitle)); */

    let minLabVisit = this.labWithDmData.getCol('VISITDY').stats[ 'min' ];
    let minVisitName = this.labWithDmData
      .groupBy([ 'VISITDY', 'VISIT' ])
      .where(`VISITDY = ${minLabVisit}`)
      .aggregate()
      .get('VISIT', 0);
    this.bl = minVisitName;
    this.uniqueLabValues = Array.from(getUniqueValues(this.labWithDmData, 'LBTEST'));
    this.uniqueVisits = Array.from(getUniqueValues(this.labWithDmData, 'VISIT'));

    this.getTop5PValues(minLabVisit);
    this.updateBoxPlots(viewerTitle, this.selectedSplitBy);

    let blVisitChoices = ui.choiceInput('BL', this.bl, this.uniqueVisits);
    blVisitChoices.onChanged((v) => {
      this.bl = blVisitChoices.value;
      this.updateBoxPlots(viewerTitle, this.selectedSplitBy);
    });

    let splitByChoices = ui.choiceInput('Split by', this.selectedSplitBy, this.splitBy);
    splitByChoices.onChanged((v) => {
      this.selectedSplitBy = splitByChoices.value;
      this.updateBoxPlots(viewerTitle, this.selectedSplitBy);
    });

    let selectBiomarkers = ui.bigButton('Select biomarkers', () => {
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
          this.updateBoxPlots(viewerTitle, this.selectedSplitBy);
        })
        .show();
    });


    //this.root.className = 'grok-view ui-box';

    this.root.append(ui.div([
      ui.divH([ blVisitChoices.root, splitByChoices.root, selectBiomarkers ]),
      ui.block([this.boxPlotDiv])
    ]));
    
    /*this.root.append(ui.splitV([
     // ui.box(plot.root, { style: { maxHeight: '200px' } }),
      ui.box(ui.divH([ blVisitChoices.root, splitByChoices.root, selectBiomarkers ]), { style: { maxHeight: '100px' } }),
      this.boxPlotDiv
    ]))*/
  }

  private updateBoxPlots(viewerTitle: any, category: string) {
    if (this.selectedLabValues && this.bl) {
      this.boxPlots = [];
      this.selectedLabValues.forEach(it => {
        let df = createBaselineEndpointDataframe(study.domains.lb, study.domains.dm, [category], it, this.bl, '', 'VISIT', `${it}_BL`);
        const plot = DG.Viewer.boxPlot(df, {
          category: category,
          value: `${it}_BL`,
          labelOrientation: 'Horz'
        });
        plot.root.prepend(ui.divText(it, viewerTitle));
        const boxPlot = Array.from(getUniqueValues(df, category)).length > 3 ? ui.block([plot.root]) : ui.block50([plot.root]);
        this.boxPlots.push(boxPlot);
      })
      updateDivInnerHTML(this.boxPlotDiv, ui.block(this.boxPlots));
    }
  }

  private getTop5PValues(visit: number){
    let pValuesArray = [];
    this.uniqueLabValues.forEach(item => {
      const valueData = this.labWithDmData
        .groupBy([ 'USUBJID', 'VISITDY', 'LBTEST', 'LBORRES', TREATMENT_ARM ])
        .where(`LBTEST = ${item} and VISITDY = ${visit}`)
        .aggregate();
      const valuesByArm = valueData.groupBy([ TREATMENT_ARM ]).getGroups();
      const dataForAnova = [];
      Object.values(valuesByArm).forEach(it => {
        const labResults = it.getCol('LBORRES').getRawData();
        dataForAnova.push(Array.from(labResults));
      })
      const pValue = jStat.anovaftest(...dataForAnova);
      pValuesArray.push({labValue: item, pValue: pValue});

    });
    pValuesArray.sort((a, b) => (b.pValue - a.pValue));
    this.selectedLabValues = pValuesArray.slice(0, 5).map(it => it.labValue);
  }
}