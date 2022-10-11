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
  maxSubstitutions: number;
  minActivityDelta: number;
  _titleHost = ui.divText('SAR Viewer', {id: 'pep-viewer-title'});
  initialized = false;
  isPropertyChanging: boolean = false;
  _isVertical = false;
  isModeChanging = false;

  constructor() {
    super();

    this.scaling = this.string('scaling', 'none', {choices: ['none', 'lg', '-lg']});
    this.bidirectionalAnalysis = this.bool('bidirectionalAnalysis', false);
    this.maxSubstitutions = this.int('maxSubstitutions', 1);
    this.minActivityDelta = this.float('minActivityDelta', 0);
  }

  get name(): string {return '';}

  async onTableAttached(): Promise<void> {
    super.onTableAttached();
    this.sourceGrid = this.view?.grid ?? (grok.shell.v as DG.TableView).grid;
    this.model = await PeptidesModel.getInstance(this.dataFrame);
    this.helpUrl = '/help/domains/bio/peptides.md';

    this.initProperties();
  }

  initProperties(): void {
    const props = this.model.usedProperties;
    IS_PROPERTY_CHANGING = true;
    for (const [propName, propVal] of Object.entries(props))
      this.props.set(propName, propVal as any as object);
    IS_PROPERTY_CHANGING = false;
  }

  detach(): void {this.subs.forEach((sub) => sub.unsubscribe());}

  get state(): string {
    return this.dataFrame.getTag(C.TAGS.SAR_MODE) ?? '10';
  }
  set state(s: string) {
    this.dataFrame.setTag(C.TAGS.SAR_MODE, s);
  }

  render(refreshOnly = false): void {
    if (!this.initialized)
      return;
    if (!refreshOnly) {
      $(this.root).empty();
      let switchHost = ui.div();
      if (this.name == 'MC') {
        const mutationCliffsMode = ui.boolInput('', this.state[0] === '1', () => {
          if (this.isModeChanging)
            return;
          this.isModeChanging = true;
          invariantMapMode.value = !invariantMapMode.value;
          this.isModeChanging = false;
          this._titleHost.innerText = 'Mutation Cliffs';
          this.model.isInvariantMap = false;
          this.viewerGrid.invalidate();
        });
        mutationCliffsMode.addPostfix('Mutation Cliffs');
        const invariantMapMode = ui.boolInput('', this.state[1] === '1', () => {
          if (this.isModeChanging)
            return;
          this.isModeChanging = true;
          mutationCliffsMode.value = !mutationCliffsMode.value;
          this.isModeChanging = false;
          this._titleHost.innerText = 'Invariant Map';
          this.model.isInvariantMap = true;
          this.viewerGrid.invalidate();
        });
        invariantMapMode.addPostfix('Invariant Map');
        const setDefaultProperties = (input: DG.InputBase): void => {
          $(input.root).find('.ui-input-editor').css('margin', '0px').attr('type', 'radio');
          $(input.root).find('.ui-input-description').css('padding', '0px').css('padding-left', '5px');
        };
        setDefaultProperties(mutationCliffsMode);
        setDefaultProperties(invariantMapMode);
        $(mutationCliffsMode.root).css('padding-right', '10px').css('padding-left', '5px');

        switchHost = ui.divH([mutationCliffsMode.root, invariantMapMode.root]);
        switchHost.style.position = 'absolute';
      }
      const viewerRoot = this.viewerGrid.root;
      viewerRoot.style.width = 'auto';
      this.root.appendChild(ui.divV([ui.divH([switchHost, this._titleHost]), viewerRoot]));
    }
    this.viewerGrid?.invalidate();
  }

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

    this.model.updateDefault();
    this.render(true);
  }
}

/**
 * Structure-activity relationship viewer.
 */
export class MutationCliffsViewer extends SARViewerBase {
  _titleHost = ui.divText('Mutation Cliffs', {id: 'pep-viewer-title'});
  _name = 'MC';
  _isVertical = false;

  constructor() {super();}

  get name(): string {return this._name;}

  async onTableAttached(): Promise<void> {
    await super.onTableAttached();
    this.model.mutationCliffsViewer ??= this;

    this.subs.push(this.model.onMutationCliffsGridChanged.subscribe((data) => {
      this.viewerGrid = data;
      this.render();
    }));

    this.model.updateDefault();
    this.viewerGrid = this.model.mutationCliffsGrid;
    this.initialized = true;
    this.render();
  }

  isInitialized(): DG.Grid {return this.model?.mutationCliffsGrid;}

  //1. debouncing in rxjs; 2. flags?
  onPropertyChanged(property: DG.Property): void {
    if (!this.isInitialized() || IS_PROPERTY_CHANGING)
      return;

    if (property.name == 'invariantMap')
      this._titleHost = ui.divText(property.get(this) ? 'Invariant Map' : 'Mutation Cliffs', {id: 'pep-viewer-title'});

    super.onPropertyChanged(property);
    IS_PROPERTY_CHANGING = true;
    this.model.syncProperties(true);
    IS_PROPERTY_CHANGING = false;
  }
}

/** Vertical structure activity relationship viewer. */
export class MostPotentResiduesViewer extends SARViewerBase {
  _name = 'MPR';
  _titleHost = ui.divText('Most Potent Residues', {id: 'pep-viewer-title'});
  _isVertical = true;

  constructor() {
    super();
  }

  get name(): string {return this._name;}

  async onTableAttached(): Promise<void> {
    await super.onTableAttached();
    this.model.mostPotentResiduesViewer ??= this;

    this.subs.push(this.model.onMostPotentResiduesGridChanged.subscribe((data) => {
      this.viewerGrid = data;
      this.render();
    }));

    this.model.updateDefault();
    this.viewerGrid = this.model.mostPotentResiduesGrid;

    this.initialized = true;
    this.render();
  }

  isInitialized(): DG.Grid {return this.model?.mostPotentResiduesGrid;}

  onPropertyChanged(property: DG.Property): void {
    if (!this.isInitialized() || IS_PROPERTY_CHANGING)
      return;

    super.onPropertyChanged(property);
    IS_PROPERTY_CHANGING = true;
    this.model.syncProperties(false);
    IS_PROPERTY_CHANGING = false;
  }
}
