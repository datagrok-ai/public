import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study, ClinRow } from "../clinical-study";
import { getUniqueValues } from '../data-preparation/utils';
import { createBaselineEndpointDataframe } from '../data-preparation/data-preparation';
import { TREATMENT_ARM } from '../constants';
import { updateDivInnerHTML } from './utils';


export class MatrixesView extends DG.ViewBase {

  matrixDiv = ui.box();
  dataframeDiv = ui.box();

  uniqueLabValues = Array.from(getUniqueValues(study.domains.lb, 'LBTEST'));
  uniqueVisits = Array.from(getUniqueValues(study.domains.lb, 'VISIT'));

  selectedLabValues: any;
  bl: any;
  matrixDataframe: DG.DataFrame;


  constructor(name) {
    super(name);
    this.name = name;

    this.createCorrelationMatrixDataframe();
    this.uniqueLabValues = Array.from(getUniqueValues(study.domains.lb, 'LBTEST'));
    this.uniqueVisits = Array.from(getUniqueValues(study.domains.lb, 'VISIT'));

    this.selectedLabValues = this.uniqueLabValues;
    this.bl = this.uniqueVisits[0];

    let blVisitChoices = ui.choiceInput('BL', this.bl, this.uniqueVisits);
    blVisitChoices.onChanged((v) => {
      this.bl = blVisitChoices.value;
      this.updateMarix();
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
                this.updateMarix();
              })
              .show();
     });

    
    this.root.className = 'grok-view ui-box';
    this.root.append(ui.splitV([
      ui.box(ui.divH([blVisitChoices.root, selectBiomarkers]), {style:{maxHeight:'100px'}}),
      this.matrixDiv
    ]))
    this.updateMarix();

  }

  private updateMarix(){
    if(this.selectedLabValues && this.bl) {
      this.matrixDataframe.rows.match(`VISIT = ${this.bl}`).filter();
      let filteredMatrixDataframe = this.matrixDataframe.clone(this.matrixDataframe.filter, this.selectedLabValues.map(it => `${it} avg(LBSTRESN)`)); 
    
      filteredMatrixDataframe.plot.fromType(DG.VIEWER.CORR_PLOT).then((v: any) => {
        let container = ui.splitV([
            v.root
        ])
        updateDivInnerHTML(this.matrixDiv, container);
      }); 
    
    }
  }

  private createCorrelationMatrixDataframe() {
    let df = study.domains.lb.clone();
    this.matrixDataframe = df
      .groupBy([ 'USUBJID', 'VISIT' ])
      .pivot('LBTEST')
      .avg('LBSTRESN')
      .aggregate();
  }
}