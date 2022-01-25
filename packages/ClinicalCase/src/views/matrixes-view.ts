import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";
import { updateDivInnerHTML } from './utils';
import { _package } from '../package';
import { getUniqueValues } from '../data-preparation/utils';
import { LAB_RES_N, LAB_TEST, VISIT_NAME, SUBJECT_ID, VS_TEST, VS_RES_N } from '../columns-constants';
import { ClinicalCaseViewBase } from '../model/ClinicalCaseViewBase';


export class MatrixesView extends ClinicalCaseViewBase {

  matrixPlot: any;
  martixPlotDiv = ui.box();
  uniqueValues = {};
  uniqueVisits: any;

  selectedValuesByDomain = {};
  selectedValues = [];
  bl: any;
  matrixDataframe: DG.DataFrame;
  domains = ['lb', 'vs'];
  domainFields = {'lb': {'test': LAB_TEST, 'res': LAB_RES_N}, 'vs': {'test': VS_TEST, 'res': VS_RES_N}};
  initialDataframe: DG.DataFrame;

  constructor(name) {
    super({});
    this.name = name;
    this.helpUrl = `${_package.webRoot}/views_help/correlation_matrix.md`;
  }

  loaded = false;

  createView(): void {
    this.domains = this.domains.filter(it => study.domains[it] !== null);
    this.domains.forEach(it => {
      let df = study.domains[it].clone(null, [SUBJECT_ID, VISIT_NAME, this.domainFields[it]['test'], this.domainFields[it]['res']]);
      df.getCol(this.domainFields[it]['test']).name = 'test';
      df.getCol(this.domainFields[it]['res']).name = 'res';
      if (!this.initialDataframe) {
        this.initialDataframe = df;
      } else {
        this.initialDataframe.append(df, true);
      }
    });
    this.createCorrelationMatrixDataframe(this.initialDataframe);
    this.domains.forEach(it => {
     this.uniqueValues[it] = Array.from(getUniqueValues(study.domains[it], this.domainFields[it]['test']));
    });
    this.uniqueVisits = Array.from(getUniqueValues(this.initialDataframe, VISIT_NAME));

    let topNum = 20;
    Object.keys(this.uniqueValues).forEach(key => {
      this.selectedValuesByDomain[key] = [];
      if (topNum > 0){
        if(this.uniqueValues[key].length > topNum){
          this.selectedValuesByDomain[key] = this.uniqueValues[key].slice(0, topNum);
          topNum = 0;
        } else {
          this.selectedValuesByDomain[key] = this.uniqueValues[key].slice(0, this.uniqueValues[key].length);
          topNum = topNum - this.uniqueValues[key].length;
        }
      }
    });
    Object.keys(this.selectedValuesByDomain).forEach(key => this.selectedValues = this.selectedValues.concat(this.selectedValuesByDomain[key]));

    this.bl = this.uniqueVisits[0];

    let blVisitChoices = ui.choiceInput('Baseline', this.bl, this.uniqueVisits);
    blVisitChoices.onChanged((v) => {
      this.bl = blVisitChoices.value;
      this.updateMarixPlot();
    });

    let selectBiomarkers = ui.iconFA('cog', () => {
      let multichoices = {};
      this.domains.forEach(domain => {
        let valuesMultiChoices = ui.multiChoiceInput('', this.selectedValuesByDomain[domain], this.uniqueValues[domain])
        valuesMultiChoices.onChanged((v) => {
          this.selectedValues = [];
          this.selectedValuesByDomain[domain] = valuesMultiChoices.value;
          Object.keys(this.selectedValuesByDomain).forEach(key => this.selectedValues = this.selectedValues.concat(this.selectedValuesByDomain[key]));
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
          this.updateMarixPlot();
        })
        .show();
    });

    this.root.className = 'grok-view ui-box';
    this.root.append(this.martixPlotDiv);
   // this.root.style.marginTop = '15px';
    this.setRibbonPanels([
      [
        blVisitChoices.root
      ],
      [
        selectBiomarkers
      ]
    ]);
    this.updateMarixPlot();
  }


  private updateMarixPlot() {
    if (this.selectedValues && this.bl) {
      let filteredDataframe = this.matrixDataframe.clone(null, this.selectedValues.map(it => `${it} avg(res)`).concat([SUBJECT_ID, VISIT_NAME]));
      filteredDataframe = filteredDataframe
      .groupBy(filteredDataframe.columns.names())
      .where(`${VISIT_NAME} = ${this.bl}`)
      .aggregate();
      filteredDataframe.plot.fromType(DG.VIEWER.CORR_PLOT).then((v: any) => {
        this.matrixPlot = v;
        this.root.className = 'grok-view ui-box';
        updateDivInnerHTML(this.martixPlotDiv, this.matrixPlot.root);
      });
    }
  }

  private createCorrelationMatrixDataframe(df: DG.DataFrame) {
    let dfForPivot = df.clone();
    this.matrixDataframe = dfForPivot
      .groupBy([SUBJECT_ID, VISIT_NAME])
      .pivot('test')
      .avg('res')
      .aggregate();
  }
}