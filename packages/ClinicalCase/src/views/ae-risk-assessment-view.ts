import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study, ClinRow } from "../clinical-study";
import { createAERiskAssessmentDataframe } from '../data-preparation/data-preparation';

export class AERiskAssessmentView extends DG.ViewBase {

  riskAssessmentDataframe: DG.DataFrame;

  constructor() {
    super();

    this.riskAssessmentDataframe = createAERiskAssessmentDataframe(study.domains.ae, study.domains.ex);
    let grid = this.riskAssessmentDataframe.plot.grid();

    this.root.appendChild(
        ui.divH([ grid.root ], { style: { height: '100%', width: '100%' } }));

  } 
}