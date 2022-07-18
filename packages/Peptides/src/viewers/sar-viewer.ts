import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import * as C from '../utils/constants';
import {PeptidesModel} from '../model';

let IS_PROPERTY_CHANGING = false;

export class SARViewerBase extends DG.JsViewer {
  tempName!: string;
  viewerGrid!: DG.Grid;
  sourceGrid!: DG.Grid;
  model!: PeptidesModel;
  scaling: string;
  bidirectionalAnalysis: boolean;
  showSubstitution: boolean;
  maxSubstitutions: number;
  minActivityDelta: number;
  _titleHost = ui.divText('SAR Viewer', {id: 'pep-viewer-title'});
  initialized = false;
  isPropertyChanging: boolean = false;

  constructor() {
    super();

    this.scaling = this.string('scaling', 'none', {choices: ['none', 'lg', '-lg']});
    this.bidirectionalAnalysis = this.bool('bidirectionalAnalysis', false);
    this.showSubstitution = this.bool('showSubstitution', true);
    this.maxSubstitutions = this.int('maxSubstitutions', 2);
    this.minActivityDelta = this.float('minActivityDelta', 1);
  }

  async onTableAttached(): Promise<void> {
    super.onTableAttached();
    this.sourceGrid = this.view?.grid ?? (grok.shell.v as DG.TableView).grid;
    this.model = await PeptidesModel.getInstance(this.dataFrame);
    // this.model.init(this.dataFrame);
    this.helpUrl = '/help/domains/bio/peptides.md';
    // await this.requestDataUpdate();

    this.initProperties();
  }

  initProperties() {
    const props = this.model.usedProperties;
    IS_PROPERTY_CHANGING = true;
    for (const [propName, propVal] of Object.entries(props))
      this.props.set(propName, propVal as any as object);
    IS_PROPERTY_CHANGING = false;
  }

  detach(): void {this.subs.forEach((sub) => sub.unsubscribe());}

  render(refreshOnly = false): void {
    if (!this.initialized)
      return;
    if (!refreshOnly) {
      $(this.root).empty();
      const viewerRoot = this.viewerGrid.root;
      viewerRoot.style.width = 'auto';
      this.root.appendChild(ui.divV([this._titleHost, viewerRoot]));
    }
    this.viewerGrid?.invalidate();
  }

  // async requestDataUpdate(): Promise<void> {
  //   await this.model.updateData(this.scaling, this.sourceGrid, this.bidirectionalAnalysis,
  //     this.minActivityDelta, this.maxSubstitutions, this.showSubstitution);
  // }

  onPropertyChanged(property: DG.Property): void {
    super.onPropertyChanged(property);
    this.dataFrame.tags[property.name] = `${property.get(this)}`;
    if (!this.initialized || IS_PROPERTY_CHANGING)
      return;

    const propName = property.name;

    if (propName === 'scaling' && typeof this.dataFrame !== 'undefined') {
      const activityCol = this.dataFrame.columns.bySemType(C.SEM_TYPES.ACTIVITY)!;
      const minActivity = activityCol.stats.min;
      if (minActivity && minActivity <= 0 && this.scaling !== 'none') {
        grok.shell.warning(`Could not apply ${this.scaling}: ` +
          `activity column ${activityCol.name} contains zero or negative values, falling back to 'none'.`);
        property.set(this, 'none');
        return;
      }
    }

    if (!this.showSubstitution && ['maxSubstitutions', 'activityLimit'].includes(propName))
      return;

    // await this.requestDataUpdate();
    this.model.updateDefault();
    this.render(true);
  }
}

/**
 * Structure-activity relationship viewer.
 */
export class SARViewer extends SARViewerBase {
  _titleHost = ui.divText('Monomer-Positions', {id: 'pep-viewer-title'});
  _name = 'Structure-Activity Relationship';

  constructor() {super();}

  get name(): string {return this._name;}

  async onTableAttached(): Promise<void> {
    await super.onTableAttached();
    this.model.sarViewer ??= this;
    // this.dataFrame.temp['sarViewer'] = this;

    this.subs.push(this.model.onSARGridChanged.subscribe((data) => {
      this.viewerGrid = data;
      this.render();
    }));

    this.model.updateDefault(true);
    this.viewerGrid = this.model._sarGrid;
    this.initialized = true;
    this.render();
  }

  isInitialized(): DG.Grid {return this.model?._sarGrid;}

  //1. debouncing in rxjs; 2. flags?
  onPropertyChanged(property: DG.Property): void {
    if (!this.isInitialized() || IS_PROPERTY_CHANGING)
      return;

    super.onPropertyChanged(property);
    IS_PROPERTY_CHANGING = true;
    this.model.syncProperties(true);
    IS_PROPERTY_CHANGING = false;
  }
}

/** Vertical structure activity relationship viewer. */
export class SARViewerVertical extends SARViewerBase {
  _name = 'Sequence-Activity relationship';
  _titleHost = ui.divText('Most Potent Residues', {id: 'pep-viewer-title'});

  constructor() {
    super();
  }

  get name(): string {return this._name;}

  async onTableAttached(): Promise<void> {
    await super.onTableAttached();
    this.model.sarViewerVertical ??= this;

    this.subs.push(this.model.onSARVGridChanged.subscribe((data) => {
      this.viewerGrid = data;
      this.render();
    }));
    
    this.model.updateDefault(true);
    this.viewerGrid = this.model._sarVGrid;

    this.initialized = true;
    this.render();
  }

  isInitialized(): DG.Grid {return this.model?._sarVGrid;}

  onPropertyChanged(property: DG.Property): void {
    if (!this.isInitialized() || IS_PROPERTY_CHANGING)
      return;

    super.onPropertyChanged(property);
    IS_PROPERTY_CHANGING = true;
    this.model.syncProperties(false);
    IS_PROPERTY_CHANGING = false;
  }
}
