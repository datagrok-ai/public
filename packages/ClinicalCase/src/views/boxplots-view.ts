import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study, ClinRow } from "../clinical-study";
import { getUniqueValues } from '../data-preparation/utils';
import { createBaselineEndpointDataframe } from '../data-preparation/data-preparation';
import { TREATMENT_ARM } from '../constants';
import { updateDivInnerHTML } from './utils';


export class BoxPlotsView extends DG.ViewBase {

  boxPlotDiv = ui.box();
  boxPlots = [];

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

    const plot = DG.Viewer.boxPlot(study.domains.dm, {
      category: 'ARM',
      value: 'AGE',
    });
    plot.root.prepend(ui.divText('Age distribution by treatment arm', viewerTitle));
  

    this.uniqueLabValues = Array.from(getUniqueValues(study.domains.lb, 'LBTEST'));
    this.uniqueVisits = Array.from(getUniqueValues(study.domains.lb, 'VISIT'));

    this.selectedLabValues = null;
    this.bl = '';

    let blVisitChoices = ui.choiceInput('BL', this.bl, this.uniqueVisits);
    blVisitChoices.onChanged((v) => {
      this.bl = blVisitChoices.value;
      this.updateBoxPlots(viewerTitle);
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
                this.updateBoxPlots(viewerTitle);
              })
              .show();
     });

    
    this.root.className = 'grok-view ui-box';
    this.root.append(ui.splitV([
      ui.box(plot.root, {style:{maxHeight:'200px'}}),
      ui.box(ui.divH([blVisitChoices.root, selectBiomarkers]),{style:{maxHeight:'100px'}}),
      this.boxPlotDiv
    ]))
  }

  private updateBoxPlots(viewerTitle: any){
    if(this.selectedLabValues && this.bl) {
      this.boxPlots = [];
      this.selectedLabValues.forEach(it => {
        let df = createBaselineEndpointDataframe(study.domains.lb,study.domains.dm,it, this.bl, '', 'VISIT',`${it}_BL`);
        const plot = DG.Viewer.boxPlot(df, {
          category: TREATMENT_ARM,
          value: `${it}_BL`,
        });
        plot.root.prepend(ui.divText(it, viewerTitle));
        this.boxPlots.push(plot.root);
      })
      updateDivInnerHTML(this.boxPlotDiv, ui.splitV(this.boxPlots));
    }
  }
}