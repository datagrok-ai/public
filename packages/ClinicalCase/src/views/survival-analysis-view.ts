import * as DG from "datagrok-api/dg";
import * as grok from 'datagrok-api/grok';
import * as ui from "datagrok-api/ui";
import { getColNames } from "../data-preparation/utils";

export class SurvivalAnalysisView extends DG.ViewBase {

survivalPlotDiv = ui.div();
covariatesPlotDiv = ui.div();
survivalColumns: string[];
confIntervals = [0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99];
confInterval = 0.7;
strata = '';



  constructor() {
    super();

    let survivalDataframe = grok.shell.table('survival');
    this.survivalColumns = getColNames(survivalDataframe);
    let survivalChoices = [''].concat(this.survivalColumns);


    let confIntChoices = ui.choiceInput('Confidence Intreval', this.confIntervals[0], this.confIntervals);
    confIntChoices.onChanged((v) => {
      this.confInterval = confIntChoices.value;
      grok.functions.call(
        "Clinicalcase:survivalPlot", {
          "inputStrata": this.strata,
          "confInt": this.confInterval.toString()
        }).then((result) => {
          this.updatePlotDiv(result, this.survivalPlotDiv);
      });     
    });

    let strataChoices = ui.choiceInput('Strata', survivalChoices[0], survivalChoices);
    strataChoices.onChanged((v) => {
      this.strata = strataChoices.value;
      grok.functions.call(
        "Clinicalcase:survivalPlot", {
          "inputStrata": this.strata,
          "confInt": this.confInterval.toString()
        }).then((result) => {
          this.updatePlotDiv(result, this.survivalPlotDiv);
      });     
    });

/*     let covariatesChoices = ui.multiChoiceInput('Value', ['A'], ['A', 'B', 'C']);
    covariatesChoices.onChanged((v) => {
      grok.functions.call(
        "Clinicalcase:covariates", {
          "coVariates": "trt"
        }).then((result) => {
          this.updatePlotDiv(result, this.covariatesPlotDiv);
      });     
    }); */

    grok.functions.call(
        "Clinicalcase:survivalPlot", {
          "inputStrata": this.strata,
          "confInt": this.confInterval.toString()

        }).then((survivalResult) => {

        grok.functions.call(
          "Clinicalcase:covariates", {
            "coVariates": "trt"
        }).then((covariatesResult) => {

          this.root.appendChild(
            ui.splitV([ 
              ui.tabControl({'Survival chart': 
                ui.divH([
                  ui.divV([
                    confIntChoices.root,
                    strataChoices.root
                  ]), 
                  this.survivalPlotDiv]), 
              'Co-variates': 
              this.covariatesPlotDiv}).root
            ], { style: { width: '100%', height: '100%' } })
          );
          this.updatePlotDiv(survivalResult, this.survivalPlotDiv);
          this.updatePlotDiv(covariatesResult, this.covariatesPlotDiv);

        });

      });

     

  } 

  private updatePlotDiv(img: string, div: HTMLDivElement) {
    div.innerHTML = '';
    div.append(ui.image(`data:image/png;base64,${img}`, 700, 700));
  }
}