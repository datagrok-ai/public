import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import * as C from '../utils/constants';
import {PeptidesModel} from '../model';

export class SARViewerBase extends DG.JsViewer {
  tempName!: string;
  viewerGrid!: DG.Grid;
  sourceGrid!: DG.Grid;
  model!: PeptidesModel;
  initialized = false;
  isPropertyChanging: boolean = false;
  _isVertical = false;
  isModeChanging = false;

  constructor() {
    super();
  }

  get name(): string {return '';}

  async onTableAttached(): Promise<void> {
    super.onTableAttached();
    this.sourceGrid = this.view?.grid ?? (grok.shell.v as DG.TableView).grid;
    this.model = await PeptidesModel.getInstance(this.dataFrame);
    this.helpUrl = '/help/domains/bio/peptides.md';
  }

  detach(): void {this.subs.forEach((sub) => sub.unsubscribe());}

  get isMutationCliffsMode(): string {
    return this.dataFrame.getTag(C.TAGS.SAR_MODE) ?? '1';
  }
  set isMutationCliffsMode(s: string) {
    this.dataFrame.setTag(C.TAGS.SAR_MODE, s);
  }

  render(refreshOnly = false): void {
    if (!this.initialized)
      return;
    if (!refreshOnly) {
      $(this.root).empty();
      let switchHost = ui.div();
      if (this.name == 'MC') {
        const mutationCliffsMode = ui.boolInput('', this.isMutationCliffsMode === '1', () => {
          if (this.isModeChanging)
            return;
          this.isModeChanging = true;
          invariantMapMode.value = !invariantMapMode.value;
          this.isMutationCliffsMode = '1';
          this.isModeChanging = false;
          this.model.isInvariantMap = false;
          this.viewerGrid.invalidate();
        });
        mutationCliffsMode.addPostfix('Mutation Cliffs');
        const invariantMapMode = ui.boolInput('', this.isMutationCliffsMode === '0', () => {
          if (this.isModeChanging)
            return;
          this.isModeChanging = true;
          mutationCliffsMode.value = !mutationCliffsMode.value;
          this.isMutationCliffsMode = '0';
          this.isModeChanging = false;
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

        switchHost = ui.divH([mutationCliffsMode.root, invariantMapMode.root], {id: 'pep-viewer-title'});
        $(switchHost).css('width', 'auto').css('align-self', 'center');
      }
      const viewerRoot = this.viewerGrid.root;
      viewerRoot.style.width = 'auto';
      this.root.appendChild(ui.divV([switchHost, viewerRoot]));
    }
    this.viewerGrid?.invalidate();
  }

  onPropertyChanged(property: DG.Property): void {
    super.onPropertyChanged(property);

    if (!this.initialized)
      return;

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
    if (!this.isInitialized())
      return;

    super.onPropertyChanged(property);
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
    if (!this.isInitialized())
      return;

    super.onPropertyChanged(property);
  }
}
