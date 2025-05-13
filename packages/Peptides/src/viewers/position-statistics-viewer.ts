/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import wu from 'wu';
import $ from 'cash-dom';
import {PeptideUtils} from '../peptideUtils';

export const POSITION_HIDDEN_NAME = '~sequence_position_monomers';

export class SequencePositionStatsViewer extends DG.JsViewer {
  public position: number;
  public sequenceColumnName: string;
  private _positionColumn?: DG.Column;
  private _boxPlotViewer?: DG.Viewer;
  public valueColumnName: string;

  constructor() {
    super();
    this.position = this.int('position', 0, {nullable: false});
    this.sequenceColumnName = this.column('sequence', {semType: DG.SEMTYPE.MACROMOLECULE, nullable: false});
    this.valueColumnName = this.column('value', {columnTypeFilter: 'numerical', nullable: false});
    grok.events.onContextMenu.subscribe((e) => {
      if (e.causedBy && e.causedBy.target && this._boxPlotViewer?.root.contains(e.causedBy.target)) {
        e.causedBy.preventDefault();
        e.causedBy.stopPropagation();
        e.causedBy.stopImmediatePropagation();
      }
    });
  }

  onTableAttached(): void {
    super.onTableAttached();
    if (this.dataFrame.columns.bySemType(DG.SEMTYPE.MACROMOLECULE) == null) {
      grok.shell.error('No sequence column found');
      throw new Error('No sequence column found');
    }
    const positionColumnName = this.dataFrame.columns.getUnusedName(POSITION_HIDDEN_NAME);

    this._positionColumn = this.dataFrame.columns.addNewString(
      positionColumnName);
    this._positionColumn.semType = 'Monomer';
    this._positionColumn.setTag('cell.renderer', 'Monomer');

    this.getProperty('sequenceColumnName')!.set(this, this.dataFrame.columns.bySemType(DG.SEMTYPE.MACROMOLECULE)!.name);
    this.getProperty('valueColumnName')!.set(this, wu(this.dataFrame.columns.numerical).next().value.name);
    this.getProperty('position')!.set(this, 0);
  }

  render(): void {
    if (this.dataFrame == null || !this.sequenceColumnName || (this.position ?? -1) < 0 || !this._positionColumn || !this.valueColumnName)
      return;

    $(this.root).empty();
    this._boxPlotViewer?.detach();
    const seqHelper = PeptideUtils.getSeqHelper();
    const sequenceColumn = this.dataFrame.col(this.sequenceColumnName)!;
    const seqHandler = seqHelper.getSeqHandler(sequenceColumn);

    this._positionColumn.init((i) => {
      try {
        return seqHandler.getSplitted(i).getCanonical(this.position);
      } catch (_) {
        return seqHandler.defaultGapOriginal;
      }
    });

    this._boxPlotViewer = this.dataFrame.plot.box({categoryColumnNames: [this._positionColumn.name], plotStyle: 'violin',
      valueColumnName: this.valueColumnName, colorColumnName: this._positionColumn.name, showColorSelector: false, showSizeSelector: false, showCategorySelector: false,
    });

    this._boxPlotViewer.root.style.width = '100%';
    this._boxPlotViewer.root.style.height = '100%';
    this.root.appendChild(this._boxPlotViewer.root);
  }

  onPropertyChanged(property: DG.Property | null): void {
    if (this.dataFrame == null || !this.sequenceColumnName)
      return;

    this.render();
  }
}
