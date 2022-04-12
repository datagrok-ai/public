import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {CHEM_SIMILARITY_METRICS} from '@datagrok-libraries/utils/src/similarity-metrics';
import {updateDivInnerHTML} from '../utils/ui-utils';

export class ChemSearchBaseViewer extends DG.JsViewer {
  distanceMetric: string;
  limit: number;
  fingerprint: string;
  metricsProperties = ['distanceMetric', 'fingerprint'];
  fingerprintChoices = ['Morgan', 'RDKit', 'Pattern'];
  sizesMap = {
    'small': {height: 60, width: 120},
    'normal': {height: 100, width: 200},
    'large': {height: 150, width: 300}};
  size: string;

  constructor() {
    super();
    this.fingerprint = this.string('fingerprint', this.fingerprintChoices[0], {choices: this.fingerprintChoices});
    this.limit = this.int('limit', 10);
    this.distanceMetric = this.string('distanceMetric', CHEM_SIMILARITY_METRICS[0], {choices: CHEM_SIMILARITY_METRICS});
    this.size = this.string('size', Object.keys(this.sizesMap)[0], {choices: Object.keys(this.sizesMap)});
  }

  updateMetricsLink(metricsDiv: HTMLDivElement, object: any, options: any) {
    const metricsButton = ui.button(`${this.distanceMetric}/${this.fingerprint}`, () => {
      if (!grok.shell.windows.showProperties)
        grok.shell.windows.showProperties = true;
      grok.shell.o = object;
    });
    //@ts-ignore
    Object.keys(options).forEach((it) => metricsButton.style[it] = options[it]);
    updateDivInnerHTML(metricsDiv, metricsButton);
  }
}
