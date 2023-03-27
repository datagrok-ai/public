import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {CHEM_SIMILARITY_METRICS} from '@datagrok-libraries/ml/src/distance-metrics-methods';
import {updateDivInnerHTML} from '../utils/ui-utils';

const BACKGROUND = 'background';
const TEXT = 'text';
export class ChemSearchBaseViewer extends DG.JsViewer {
  name: string = '';
  distanceMetric: string;
  limit: number;
  fingerprint: string;
  metricsProperties = ['distanceMetric', 'fingerprint'];
  fingerprintChoices = ['Morgan'];
  sizesMap: {[key: string]: {[key: string]: number}} = {
    'small': {height: 60, width: 120},
    'normal': {height: 100, width: 200},
    'large': {height: 150, width: 300}};
  size: string;
  moleculeColumn?: DG.Column|null;
  moleculeColumnName: string;
  initialized: boolean = false;
  metricsDiv = ui.div('', {style: {height: '10px', display: 'flex', justifyContent: 'right'}});
  moleculeProperties: string[];
  applyColorTo: string;

  constructor(name: string, col?: DG.Column) {
    super();
    this.fingerprint = this.string('fingerprint', this.fingerprintChoices[0], {choices: this.fingerprintChoices});
    this.limit = this.int('limit', 10);
    this.distanceMetric = this.string('distanceMetric', CHEM_SIMILARITY_METRICS[0], {choices: CHEM_SIMILARITY_METRICS});
    this.size = this.string('size', Object.keys(this.sizesMap)[0], {choices: Object.keys(this.sizesMap)});
    this.moleculeColumnName = this.string('moleculeColumnName');
    this.name = name;
    this.moleculeProperties = this.columnList('moleculeProperties', []);
    this.applyColorTo = this.string('applyColorTo', BACKGROUND, {choices: [BACKGROUND, TEXT]});
    if (col) {
      this.moleculeColumn = col;
      this.moleculeColumnName = col.name!;
    }
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
      this.moleculeColumn ??= this.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE);
      this.moleculeColumnName ??= this.moleculeColumn?.name!;
      this.getProperty('limit')!.fromOptions({min: 1, max: this.dataFrame.rowCount});
      if (this.limit > this.dataFrame.rowCount)
        this.limit = this.dataFrame.rowCount;
    }
    await this.render(true);
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

  updateMetricsLink(metricsDiv: HTMLDivElement, object: any, options: {[key: string]: string}): void {
    const metricsButton = ui.button(`${this.distanceMetric}/${this.fingerprint}`, () => {
      if (!grok.shell.windows.showProperties)
        grok.shell.windows.showProperties = true;
      grok.shell.o = object;
    });
    Object.keys(options).forEach((it: any) => metricsButton.style[it] = options[it]);
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

  pickTextColorBasedOnBgColor(bgColor: string, lightColor: string, darkColor: string) {
    const color = (bgColor.charAt(0) === '#') ? bgColor.substring(1, 7) : bgColor;
    const r = parseInt(color.substring(0, 2), 16); // hexToR
    const g = parseInt(color.substring(2, 4), 16); // hexToG
    const b = parseInt(color.substring(4, 6), 16); // hexToB
    return (((r * 0.299) + (g * 0.587) + (b * 0.114)) > 186) ?
      darkColor : lightColor;
  }

  createMoleculePropertiesDiv(idx: number, similarity?: number): HTMLDivElement{
    const propsDict: {[key: string]: any} = {};
    const grid = grok.shell.tv.grid;
    if (similarity)
      propsDict['similarity'] = {val: similarity};
    for (const col of this.moleculeProperties) {  
        propsDict[col] = {val: this.moleculeColumn!.dataFrame.col(col)!.getString(idx)};
        if (this.moleculeColumn!.dataFrame.col(col)!.tags[DG.TAGS.COLOR_CODING_TYPE]) {
            propsDict[col].color = grid.cell(col, idx).color;
        }
    }
    const div = ui.divV([]);
    for (const key of Object.keys(propsDict)) {
      const label = ui.divText(`${key}`, 'similarity-prop-label');
      ui.tooltip.bind(label, key);
      const value = ui.divText(`${propsDict[key].val}`, 'similarity-prop-value');
      if (propsDict[key].color) {
        const bgColor = DG.Color.toHtml(propsDict[key].color);
        if (this.applyColorTo === BACKGROUND) {
          value.style.backgroundColor = bgColor,
          value.style.color = this.pickTextColorBasedOnBgColor(bgColor, '#FFFFFF', '#000000');;
        } else
          value.style.color = bgColor;
      }
      const item = ui.divH([
        label,
        value
      ], 'similarity-prop-item')
      div.append(item);
    }
    return div; 
  }

  getPropsColumnsNames(): string[] {
    let fingerprintTag = '';
    for (let t of this.moleculeColumn!.tags.keys()) {
      if (t.endsWith('.Column')) {
        fingerprintTag = t;
        break;
      }
    }
    return this.moleculeColumn!.dataFrame.columns.names()
      .filter((name) => name !== this.moleculeColumn!.getTag(fingerprintTag) && name !== this.moleculeColumn!.name);
  }
}
