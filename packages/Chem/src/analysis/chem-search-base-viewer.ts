import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {CHEM_SIMILARITY_METRICS} from '@datagrok-libraries/ml/src/distance-metrics-methods';
import '../../css/chem.css';
import {Fingerprint} from '../utils/chem-common';

const BACKGROUND = 'background';
const TEXT = 'text';
export const SIMILARITY = 'similarity';
export const DIVERSITY = 'diversity';
export const MAX_LIMIT = 50;
export class ChemSearchBaseViewer extends DG.JsViewer {
  isEditedFromSketcher: boolean = false;
  gridSelect: boolean = false;
  name: string = '';
  distanceMetric: string;
  limit: number;
  fingerprint: string;
  metricsProperties = ['distanceMetric', 'fingerprint'];
  fingerprintChoices = [Fingerprint.Morgan];
  sizesMap: {[key: string]: {[key: string]: number}} = {
    'small': {height: 60, width: 120},
    'normal': {height: 100, width: 200},
    'large': {height: 150, width: 300}};
  size: string;
  moleculeColumn?: DG.Column|null;
  moleculeColumnName: string;
  initialized: boolean = false;
  metricsDiv: HTMLElement | null;
  moleculeProperties: string[];
  applyColorTo: string;

  constructor(name: string, col?: DG.Column) {
    super();
    this.fingerprint = this.string('fingerprint', this.fingerprintChoices[0], {choices: this.fingerprintChoices});
    this.limit = this.int('limit', 12, {min: 1, max: MAX_LIMIT});
    this.distanceMetric = this.string('distanceMetric', CHEM_SIMILARITY_METRICS[0], {choices: CHEM_SIMILARITY_METRICS});
    this.size = this.string('size', Object.keys(this.sizesMap)[0], {choices: Object.keys(this.sizesMap)});
    this.moleculeColumnName = this.string('moleculeColumnName');
    this.name = name;
    this.moleculeProperties = this.columnList('moleculeProperties', [],
      {description: 'Adds selected fields from the grid to similarity search viewer'});
    this.applyColorTo = this.string('applyColorTo', BACKGROUND, {choices: [BACKGROUND, TEXT],
      description: 'Applies to data added via Molecule Properties control (color-code the column in the grid first)'});
    if (col) {
      this.moleculeColumn = col;
      this.moleculeColumnName = col.name!;
    }
    const header = this.name === DIVERSITY ? `Most diverse structures` : `Most similar structures`;
    this.metricsDiv = ui.divH([ui.divText(header)], 'similarity-header');
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
      const compute = this.name !== DIVERSITY;
      this.subs.push(DG.debounce(this.dataFrame.onCurrentRowChanged, 50)
        .subscribe(async (_: any) => {
          if (this.isEditedFromSketcher && !this.gridSelect)
            this.isEditedFromSketcher = false;
          await this.render(compute);
        }));
      this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50)
        .subscribe(async (_: any) => await this.render(false)));
      this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50)
        .subscribe(async (_: any) => await this.render(false)));
      this.subs.push(DG.debounce((this.dataFrame.onMetadataChanged), 50)
        .subscribe(async (_: any) => await this.render(false)));
      this.moleculeColumn = this.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE);
      this.moleculeColumnName = this.moleculeColumn?.name ?? '';
    }
    await this.render(true);
  }

  onPropertyChanged(property: DG.Property): void {
    super.onPropertyChanged(property);
    if (!this.initialized)
      return;
    if (this.metricsProperties.includes(property.name))
      this.updateMetricsLink(this, {});
    if (property.name === 'moleculeColumnName') {
      const col = this.dataFrame.col(property.get(this));
      if (col?.semType === DG.SEMTYPE.MOLECULE)
        this.moleculeColumn = col;
    }
    if (property.name === 'limit' && property.get(this) > MAX_LIMIT )
      this.limit = MAX_LIMIT;
    if (property.name === 'moleculeProperties') {
      this.render(false);
      return;
    }
    this.render();
  }

  updateMetricsLink(object: any, options: {[key: string]: string}): void {
    const metricsButton = ui.link(` ${this.distanceMetric}, ${this.fingerprint}`, () => {
      if (!grok.shell.windows.showProperties)
        grok.shell.windows.showProperties = true;
      grok.shell.o = object;
    }, 'Distance metric and fingerprint', '');
    Object.keys(options).forEach((it: any) => metricsButton.style[it] = options[it]);
    if (this.metricsDiv!.children.length > 1)
      this.metricsDiv!.removeChild(this.metricsDiv!.children[1]);
    this.metricsDiv!.appendChild(metricsButton);
  }

  async render(compute = true) {

  }

  beforeRender() {
    if (!this.initialized || !this.dataFrame)
      return false;
    if (this.dataFrame && this.moleculeColumnName &&
          this.dataFrame.col(this.moleculeColumnName)?.semType !== DG.SEMTYPE.MOLECULE) {
      grok.shell.error(`${this.moleculeColumnName} is not Molecule type or missing`);
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

  createMoleculePropertiesDiv(idx: number, refMolecule: boolean, similarity?: number): HTMLDivElement {
    const propsDict: {[key: string]: any} = {};
    const grid = grok.shell.tv.grid;
    if (similarity) {
      if (refMolecule)
        propsDict['Reference'] = {val: ''};
      else
        propsDict[SIMILARITY] = {val: similarity};
    }
    for (const col of this.moleculeProperties) {
      propsDict[col] = {val: this.moleculeColumn!.dataFrame.col(col)!.getString(idx)};
      const colorCoding = this.moleculeColumn!.dataFrame.col(col)!.tags[DG.TAGS.COLOR_CODING_TYPE];
      if (colorCoding && colorCoding !== DG.COLOR_CODING_TYPE.OFF)
        propsDict[col].color = grid.cell(col, idx).color;
    }
    //const item = ui.divH([], 'similarity-prop-item');
    const div = ui.divV([], {style: {marginTop: '5px'}});
    for (const key of Object.keys(propsDict)) {
      const labelName = key === SIMILARITY ? '' : key;
      const label = ui.divText(`${labelName}`, 'similarity-prop-label');
      const value = ui.divText(`${propsDict[key].val}`, 'similarity-prop-value');
      ui.tooltip.bind(value, key);
      if (propsDict[key].color) {
        const bgColor = DG.Color.toHtml(propsDict[key].color);
        if (this.applyColorTo === BACKGROUND) {
          value.style.backgroundColor = bgColor,
          value.style.color = this.pickTextColorBasedOnBgColor(bgColor, '#FFFFFF', '#000000'); ;
        } else
          value.style.color = bgColor;
      }
      const item = ui.divH([
        label,
        value,
      ], 'similarity-prop-item');
      div.append(item);
    }
    return div;
  }

  getPropsColumnsNames(): string[] {
    let fingerprintTag = '';
    for (const t of this.moleculeColumn!.tags.keys()) {
      if (t.endsWith('.Column')) {
        fingerprintTag = t;
        break;
      }
    }
    return this.moleculeColumn!.dataFrame.columns.names()
      .filter((name) => name !== this.moleculeColumn!.getTag(fingerprintTag) && name !== this.moleculeColumn!.name);
  }
}
