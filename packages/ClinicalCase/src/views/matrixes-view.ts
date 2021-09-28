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

  selectedLabValues = null;
  bl = '';


  constructor(name) {
    super(name);
    this.name = name;
    
    let viewerTitle = {style:{
      'color':'var(--grey-6)',
      'margin':'12px 0px 6px 12px',
      'font-size':'16px',
    }};

    this.uniqueLabValues = Array.from(getUniqueValues(study.domains.lb, 'LBTEST'));
    this.uniqueVisits = Array.from(getUniqueValues(study.domains.lb, 'VISIT'));

    this.selectedLabValues = null;
    this.bl = '';

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

  }

  private updateMarix(){
    if(this.selectedLabValues && this.bl) {
    study.labDataForCorelationMatrix.rows.match(`VISIT = ${this.bl}`).filter();
    let matrixDataframe = study.labDataForCorelationMatrix.clone(study.labDataForCorelationMatrix.filter, this.selectedLabValues.concat(['USUBJID'])) 
    updateDivInnerHTML(this.matrixDiv, grok.shell.addTableView(matrixDataframe).matrixPlot());
    
    /* matrixDataframe.plot.fromType(DG.VIEWER.MATRIX_PLOT).then((v: any) => {
        let container = ui.splitV([
            ui.box(ui.panel([matrixDataframe.plot.grid().root]), {style:{maxHeight:'200px'}}),
            v.root
        ])
        updateDivInnerHTML(this.matrixDiv, container);
      }); */
    
    }
  }
}