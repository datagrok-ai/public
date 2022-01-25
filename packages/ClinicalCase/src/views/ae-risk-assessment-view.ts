import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study, ClinRow } from "../clinical-study";
import { createAERiskAssessmentDataframe } from '../data-preparation/data-preparation';
import { checkMissingDomains, updateDivInnerHTML } from './utils';
import { TREATMENT_ARM } from '../columns-constants';
import { ClinicalCaseViewBase } from '../model/ClinicalCaseViewBase';

export class AERiskAssessmentView extends ClinicalCaseViewBase {

  riskAssessmentDataframe: DG.DataFrame;
  volcanoPlotXAxis = '';
  volcanoPlotYAxis = '-log10(p-value)';
  placeboArm = [];
  activeArm = [];
  treatmentArmOptions: string[];
  riskAssessmentDiv = ui.box();
  initialGuide: any;
  xScaleOptions = { 'Linear': 'RD', 'Log2': 'log2(RD)' };
  volcanoPlot: any;
  pValueLimit = 0.05;
  incorrectTreatmentArms = false;


  constructor(name) {
    super({});
    this.name = name;
  }

  createView(): void {
    this.volcanoPlotXAxis = Object.values(this.xScaleOptions)[0];
    this.treatmentArmOptions = study.domains.dm.col(TREATMENT_ARM).categories;
    this.initialGuide = ui.info('Please select values for treatment/placebo arms in a property panel', '', false);
    updateDivInnerHTML(this.riskAssessmentDiv, this.initialGuide);
    this.root.append(this.riskAssessmentDiv);
  }

  createVolcanoPlotDiv() {
    if (this.placeboArm.length && this.activeArm.length && !this.incorrectTreatmentArms) {
      this.riskAssessmentDataframe = createAERiskAssessmentDataframe(study.domains.ae, study.domains.dm.clone(), this.placeboArm, this.activeArm);
      this.riskAssessmentDataframe.col(this.volcanoPlotXAxis).tags[DG.TAGS.COLOR_CODING_TYPE] = 'Linear';
      this.riskAssessmentDataframe.col(this.volcanoPlotXAxis).tags[DG.TAGS.COLOR_CODING_LINEAR] = `[${DG.Color.blue}, ${DG.Color.red}]`;
      this.updateVolcanoPlot();
    } else {
      this.riskAssessmentDataframe = null;
      updateDivInnerHTML(this.riskAssessmentDiv, this.initialGuide);
    }
  }

  updateVolcanoPlot() {
    this.volcanoPlot = this.riskAssessmentDataframe.plot.scatter({
      x: this.volcanoPlotXAxis,
      y: this.volcanoPlotYAxis,
      color: this.volcanoPlotXAxis,
      showViewerFormulaLines: true
    });
    this.volcanoPlot.meta.addFormulaLine({
      formula: `\${${this.volcanoPlotYAxis}} = ${-Math.log10(this.pValueLimit)}`,
      color: "#ff0000",
      width: 0.5
    });
    updateDivInnerHTML(this.riskAssessmentDiv, ui.splitV([this.volcanoPlot.root]))
  }

  override async propertyPanel() {

    let acc = this.createAccWithTitle(this.name);
    let dfDiv = ui.div();

    let updateDf = () => { updateDivInnerHTML(dfDiv, this.riskAssessmentDataframe ? this.riskAssessmentDataframe.plot.grid().root : ''); }

    const treatmentArmPane = (arm, armToCheck, name) => {
      acc.addPane(name, () => {
        //@ts-ignore
        const armChoices = ui.multiChoiceInput('', this[arm], this.treatmentArmOptions);
        armChoices.onChanged((v) => {
          if (this[armToCheck].filter(it => armChoices.value.includes(it)).length) {
            grok.shell.error(`Some products are selected for both treatment arms. One product can be selected only for one treatment arm.`);
            this.incorrectTreatmentArms = true;
          } else {
            this.incorrectTreatmentArms = false;;
          }
          this[arm] = armChoices.value;
          this.createVolcanoPlotDiv();
          updateDf();
        });
        return armChoices.root;
      });
    }

    treatmentArmPane('placeboArm', 'activeArm', `Placebo arm`);
    treatmentArmPane('activeArm', 'placeboArm', `Active arm`);

    acc.addPane(`AE risk dataframe`, () => {
      updateDf();
      return dfDiv;
    });

    acc.addPane('Parameters', () => {
      let xChoices = ui.choiceInput('X axis type', Object.keys(this.xScaleOptions)[0], Object.keys(this.xScaleOptions));
      xChoices.onChanged((v) => {
        this.volcanoPlotXAxis = this.xScaleOptions[xChoices.value];
        this.updateVolcanoPlot();
      });
      let pValue = ui.stringInput('Significance level', this.pValueLimit.toString());
      pValue.onChanged((v) => {
        this.pValueLimit = parseFloat(pValue.value);
        if (!isNaN(this.pValueLimit)) {
          this.updateVolcanoPlot();
        }
      });
      return ui.inputs([
        xChoices,
        pValue,
      ]);
    })
    return acc.root;
  }

}