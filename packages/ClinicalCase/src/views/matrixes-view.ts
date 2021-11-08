import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";
import { ILazyLoading } from '../lazy-loading/lazy-loading';
import { checkMissingDomains } from './utils';
import { _package } from '../package';
import { getUniqueValues } from '../data-preparation/utils';
import { LAB_RES_N, LAB_TEST, LAB_VISIT_NAME, SUBJECT_ID } from '../columns-constants';
import { requiredColumnsByView } from '../constants';


export class MatrixesView extends DG.ViewBase implements ILazyLoading {

  matrixPlot: any;
  uniqueLabValues: any;
  uniqueVisits: any;

  selectedLabValues: any;
  bl: any;
  matrixDataframe: DG.DataFrame;


  constructor(name) {
    super({});
    this.name = name;
    this.helpUrl = `${_package.webRoot}/views_help/correlation_matrix.md`;
  }

  loaded: boolean;

  load(): void {
    checkMissingDomains(requiredColumnsByView[this.name], false, this);
  }

  createView(): void {
    this.createCorrelationMatrixDataframe();
    this.uniqueLabValues = Array.from(getUniqueValues(study.domains.lb, LAB_TEST));
    this.uniqueVisits = Array.from(getUniqueValues(study.domains.lb, LAB_VISIT_NAME));

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
      this.root.append(ui.box(this.matrixPlot.root));
     // this.root.style.marginTop = '15px';
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
      this.matrixDataframe.rows.match(`${LAB_VISIT_NAME} = ${this.bl}`).filter();
      let filteredMatrixDataframe = this.matrixDataframe.clone(this.matrixDataframe.filter, this.selectedLabValues.map(it => `${it} avg(${LAB_RES_N})`));
      this.matrixPlot.dataFrame = filteredMatrixDataframe;
    }
  }

  private createCorrelationMatrixDataframe() {
    let df = study.domains.lb.clone();
    this.matrixDataframe = df
      .groupBy([SUBJECT_ID, LAB_VISIT_NAME])
      .pivot(LAB_TEST)
      .avg(LAB_RES_N)
      .aggregate();
  }
}