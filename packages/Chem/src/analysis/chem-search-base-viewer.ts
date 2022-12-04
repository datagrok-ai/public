import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {CHEM_SIMILARITY_METRICS} from '@datagrok-libraries/utils/src/similarity-metrics';
import {updateDivInnerHTML} from '../utils/ui-utils';

export class ChemSearchBaseViewer extends DG.JsViewer {
  name: string = '';
  distanceMetric: string;
  limit: number;
  fingerprint: string;
  metricsProperties = ['distanceMetric', 'fingerprint'];
  fingerprintChoices = ['Morgan'];
  sizesMap = {
    'small': {height: 60, width: 120},
    'normal': {height: 100, width: 200},
    'large': {height: 150, width: 300}};
  size: string;
  moleculeColumn?: DG.Column|null;
  moleculeColumnName: string;
  initialized: boolean = false;
  metricsDiv = ui.div('', {style: {height: '10px', display: 'flex', justifyContent: 'right'}});

  constructor(name: string) {
    super();
    this.fingerprint = this.string('fingerprint', this.fingerprintChoices[0], {choices: this.fingerprintChoices});
    this.limit = this.int('limit', 10);
    this.distanceMetric = this.string('distanceMetric', CHEM_SIMILARITY_METRICS[0], {choices: CHEM_SIMILARITY_METRICS});
    this.size = this.string('size', Object.keys(this.sizesMap)[0], {choices: Object.keys(this.sizesMap)});
    this.moleculeColumnName = this.string('moleculeColumnName');
    this.name = name;
  }

  init(): void {
    this.initialized = true;
  }

  detach(): void {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  async onTableAttached(): Promise<void> {
    this.init();

    if (this.dataFrame) {
      this.subs.push(DG.debounce(this.dataFrame.onRowsRemoved, 50).subscribe(async (_: any) => await this.render()));
      const compute = this.name !== 'diversity';
      this.subs.push(DG.debounce(this.dataFrame.onCurrentRowChanged, 50)
        .subscribe(async (_: any) => await this.render(compute)));
      this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50)
        .subscribe(async (_: any) => await this.render(false)));
      this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe(async (_: any) => await this.render(false)));
      this.moleculeColumn = this.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE);
      this.moleculeColumnName = this.moleculeColumn?.name!;
      this.getProperty('limit')!.fromOptions({min: 1, max: this.dataFrame.rowCount});
    }
    await this.render();
  }

  onPropertyChanged(property: DG.Property): void {
    super.onPropertyChanged(property);
    if (!this.initialized)
      return;
    if (this.metricsProperties.includes(property.name))
      this.updateMetricsLink(this.metricsDiv, this, {fontSize: '10px', fontWeight: 'normal', paddingBottom: '15px'});
    if (property.name === 'moleculeColumnName') {
      const col = this.dataFrame.col(property.get(this))!;
      if (col.semType === DG.SEMTYPE.MOLECULE)
        this.moleculeColumn = col;
    }
    this.render();
  }

  updateMetricsLink(metricsDiv: HTMLDivElement, object: any, options: any): void {
    const metricsButton = ui.button(`${this.distanceMetric}/${this.fingerprint}`, () => {
      if (!grok.shell.windows.showProperties)
        grok.shell.windows.showProperties = true;
      grok.shell.o = object;
    });
    //@ts-ignore
    Object.keys(options).forEach((it) => metricsButton.style[it] = options[it]);
    updateDivInnerHTML(metricsDiv, metricsButton);
  }

  async render(computeData = true) {

  }

  beforeRender() {
    if (!this.initialized)
      return false;
    if (this.dataFrame && this.moleculeColumnName &&
          this.dataFrame.col(this.moleculeColumnName)!.semType !== DG.SEMTYPE.MOLECULE) {
      grok.shell.error(`${this.moleculeColumnName} is not Molecule type`);
      return false;
    }
    return true;
  }
}
