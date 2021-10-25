import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";
import { getUniqueValues } from '../data-preparation/utils';
import { ILazyLoading } from '../lazy-loading/lazy-loading';


export class MatrixesView extends DG.ViewBase implements ILazyLoading {

  matrixPlot: any;
  uniqueLabValues = Array.from(getUniqueValues(study.domains.lb, 'LBTEST'));
  uniqueVisits = Array.from(getUniqueValues(study.domains.lb, 'VISIT'));

  selectedLabValues: any;
  bl: any;
  matrixDataframe: DG.DataFrame;


  constructor(name) {
    super(name);
    this.name = name;
  }
  loaded: boolean;

  load(): void {
    
    this.createCorrelationMatrixDataframe();
    this.uniqueLabValues = Array.from(getUniqueValues(study.domains.lb, 'LBTEST'));
    this.uniqueVisits = Array.from(getUniqueValues(study.domains.lb, 'VISIT'));

    this.selectedLabValues = this.uniqueLabValues;
    this.bl = this.uniqueVisits[0];

    let blVisitChoices = ui.choiceInput('Baseline', this.bl, this.uniqueVisits);
    blVisitChoices.onChanged((v) => {
      this.bl = blVisitChoices.value;
      this.updateMarix();
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
        .add(ui.div([labValuesMultiChoices], { style: { width: '400px', height: '300px' } }))
        .onOK(() => {
          this.updateMarix();
        })
        .show();
    });

    this.matrixDataframe.plot.fromType(DG.VIEWER.CORR_PLOT).then((v: any) => {
      this.matrixPlot = v;
      this.root.className = 'grok-view ui-box';
      this.root.append(ui.box(this.matrixPlot.root, { style: { 'margin-top': '15px' } }));
      this.setRibbonPanels([
        [
          blVisitChoices.root
        ],
        [
          selectBiomarkers
        ]
      ]);
    });
  }

  private updateMarix() {
    if (this.selectedLabValues && this.bl) {
      this.matrixDataframe.rows.match(`VISIT = ${this.bl}`).filter();
      let filteredMatrixDataframe = this.matrixDataframe.clone(this.matrixDataframe.filter, this.selectedLabValues.map(it => `${it} avg(LBSTRESN)`));
      this.matrixPlot.dataFrame = filteredMatrixDataframe;
    }
  }

  private createCorrelationMatrixDataframe() {
    let df = study.domains.lb.clone();
    this.matrixDataframe = df
      .groupBy(['USUBJID', 'VISIT'])
      .pivot('LBTEST')
      .avg('LBSTRESN')
      .aggregate();
  }
}