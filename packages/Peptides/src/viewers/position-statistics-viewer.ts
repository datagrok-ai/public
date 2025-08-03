/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import wu from 'wu';
import $ from 'cash-dom';
import {PeptideUtils} from '../peptideUtils';
import {TAGS as bioTAGS, MONOMER_MOTIF_SPLITTER} from '@datagrok-libraries/bio/src/utils/macromolecule';

export const POSITION_HIDDEN_NAME = '~sequence_position_monomers';

export class SequencePositionStatsViewer extends DG.JsViewer {
  public position: number;
  public sequenceColumnName: string;
  private _positionColumn?: DG.Column;
  private _boxPlotViewer?: DG.Viewer;
  public valueColumnName: string;
  public leftMotifLength: number = 0;
  public rightMotifLength: number = 0;
  public showPositionInfo: boolean = true;
  constructor() {
    super();
    this.position = this.int('position', 0, {nullable: false});
    this.sequenceColumnName = this.column('sequence', {semType: DG.SEMTYPE.MACROMOLECULE, nullable: false});
    this.valueColumnName = this.column('value', {columnTypeFilter: 'numerical', nullable: false});
    this.leftMotifLength = this.int('leftMotifLength', 0, {nullable: false, min: 0, max: 10});
    this.rightMotifLength = this.int('rightMotifLength', 0, {nullable: false, min: 0, max: 10});
    this.showPositionInfo = this.bool('showPositionInfo', true, {nullable: false, defaultValue: true, description: 'Show position and overhangs info in the viewer header'});
    grok.events.onContextMenu.subscribe((e) => {
      if (e.causedBy && e.causedBy.target && this._boxPlotViewer?.root.contains(e.causedBy.target)) {
        e.causedBy.preventDefault();
        e.causedBy.stopPropagation();
        e.causedBy.stopImmediatePropagation();
      }
    });
  }

  getPositionFromColumn(): number {
    const col = this.dataFrame.col(this.sequenceColumnName);
    if (col == null)
      return 0;
    const positionString = col.getTag(bioTAGS.selectedPosition) ?? '1';
    const position = parseInt(positionString);
    if (Number.isNaN(position))
      return 0;
    return Math.max(0, position - 1);
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
    this.getProperty('position')!.set(this, this.getPositionFromColumn());

    this.subs.push(DG.debounce(this.dataFrame.onMetadataChanged, 200).subscribe((_) => {
      const curPosition = this.getPositionFromColumn();
      if (this.position !== curPosition)
        this.getProperty('position')!.set(this, curPosition);
    }));
  }

  render(): void {
    if (this.dataFrame == null || !this.sequenceColumnName || (this.position ?? -1) < 0 || !this._positionColumn || !this.valueColumnName)
      return;

    $(this.root).empty();
    this._boxPlotViewer?.detach();
    const seqHelper = PeptideUtils.getSeqHelper();
    const sequenceColumn = this.dataFrame.col(this.sequenceColumnName)!;
    const seqHandler = seqHelper.getSeqHandler(sequenceColumn);
    const leftOverhang = Math.min(Math.max(this.leftMotifLength ?? 0, 0), 10);
    const rightOverhang = Math.min(Math.max(this.rightMotifLength ?? 0, 0), 10);
    const start = Math.max(0, this.position - leftOverhang);
    const end = rightOverhang + this.position;
    const canonicals = Array.from({length: end - start + 1}).fill('')
      .map((_, i) => seqHandler.getMonomersAtPosition(start + i, true));
    this._positionColumn.init((i) => canonicals.map((c) => c[i]).join(MONOMER_MOTIF_SPLITTER));

    this._boxPlotViewer = this.dataFrame.plot.box({categoryColumnNames: [this._positionColumn.name], plotStyle: 'violin',
      valueColumnName: this.valueColumnName, colorColumnName: this._positionColumn.name, showColorSelector: false, showSizeSelector: false, showCategorySelector: false,
      legendVisibility: DG.VisibilityMode.Never, markerColorColumnName: this._positionColumn.name, title: 'Sequence Position Statistics',
      autoLayout: false, labelOrientation: 'Vert',
    });

    const leftOverhangInput = ui.input.int('Left Overhang', {value: leftOverhang, min: 0, max: 10, step: 1, showSlider: false, showPlusMinus: true,
      onValueChanged: (v) => {
        this.getProperty('leftMotifLength')!.set(this, leftOverhangInput.value);
      }, tooltipText: 'Left overhang motif length from the selected position',
    });
    const rightOverhangInput = ui.input.int('Right Overhang', {value: rightOverhang, min: 0, max: 10, step: 1, showSlider: false, showPlusMinus: true,
      onValueChanged: (v) => {
        this.getProperty('rightMotifLength')!.set(this, rightOverhangInput.value);
      }, tooltipText: 'Right overhang motif length from the selected position',
    });

    const descriptionDiv = ui.divH([
      leftOverhangInput.root, ui.h2(`${this.sequenceColumnName}: Position ${this.position + 1}`),
      rightOverhangInput.root,
    ], {style: {alignItems: 'center', justifyContent: 'space-around', width: '100%'}});
    if (this.showPositionInfo) {
      this.root.appendChild(descriptionDiv);
      leftOverhangInput.input.style.width = '20px';
      rightOverhangInput.input.style.width = '20px';
    }
    // setTimeout(() => {
    //   this._boxPlotViewer!.props.title = 'Sequence Position Statistics';
    //   this._boxPlotViewer!.props.description = `${this.sequenceColumnName}: Position ${this.position + 1}`;
    // }, 200);

    this._boxPlotViewer.props.statistics = ['min', 'max', 'avg', 'med', 'count'];

    this._boxPlotViewer.root.style.width = '100%';
    this._boxPlotViewer.root.style.height = '100%';
    this.root.appendChild(this._boxPlotViewer.root);
    this._boxPlotViewer.sub(this._boxPlotViewer.onPropertyValueChanged.subscribe((_e) => {
      if (this._boxPlotViewer?.props?.valueColumnName && this._boxPlotViewer?.props?.valueColumnName !== this.valueColumnName) {
        const value = this._boxPlotViewer.props.valueColumnName;
        setTimeout(() => this.getProperty('valueColumnName')!.set(this, value), 10);
      }
    }));
  }

  onPropertyChanged(property: DG.Property | null): void {
    if (this.dataFrame == null || !this.sequenceColumnName)
      return;

    this.render();
  }
}
