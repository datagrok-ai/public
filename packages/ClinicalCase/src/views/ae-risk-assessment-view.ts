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

  loaded = false;

  load(): void {
    checkMissingDomains(requiredColumnsByView[this.name], this);
 }

  createView(): void {
    this.riskAssessmentDataframe = createAERiskAssessmentDataframe(study.domains.ae, study.domains.ex, 'PLACEBO', 'XANOMELINE');
    let grid = this.riskAssessmentDataframe.plot.grid();
    let volcanoPlot = this.riskAssessmentDataframe.plot.scatter({
      x: 'RD',
      y: '-log 10(p-value)',
    });
   
    this.root.className = 'grok-view ui-box';
    this.root.append(ui.splitV([
      grid.root,
      volcanoPlot.root
    ]))

  }
}