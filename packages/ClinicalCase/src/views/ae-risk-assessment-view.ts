import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study, ClinRow } from "../clinical-study";
import { createAERiskAssessmentDataframe } from '../data-preparation/data-preparation';

export class AERiskAssessmentView extends DG.ViewBase {

  riskAssessmentDataframe: DG.DataFrame;

  constructor(name) {
    super(name);
    this.name = name;
    this.riskAssessmentDataframe = createAERiskAssessmentDataframe(study.domains.ae, study.domains.ex, 'PLACEBO', 'XANOMELINE');
    let grid = this.riskAssessmentDataframe.plot.grid();
   
    this.root.className = 'grok-view ui-box';
    this.root.appendChild(grid.root);

  } 
}