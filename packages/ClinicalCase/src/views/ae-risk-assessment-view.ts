import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study, ClinRow } from "../clinical-study";
import { createAERiskAssessmentDataframe } from '../data-preparation/data-preparation';
import { ILazyLoading } from '../lazy-loading/lazy-loading';
import { checkMissingDomains } from './utils';
import { requiredColumnsByView } from '../constants';

export class AERiskAssessmentView extends DG.ViewBase implements ILazyLoading  {

  riskAssessmentDataframe: DG.DataFrame;

  constructor(name) {
    super({});
    this.name = name;
  } 

  loaded: boolean;

  load(): void {
    checkMissingDomains(requiredColumnsByView[this.name], false, this);
 }

  createView(): void {
    this.riskAssessmentDataframe = createAERiskAssessmentDataframe(study.domains.ae, study.domains.ex, 'PLACEBO', 'XANOMELINE');
    let grid = this.riskAssessmentDataframe.plot.grid();
   
    this.root.className = 'grok-view ui-box';
    this.root.appendChild(grid.root);

  }
}