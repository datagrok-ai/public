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
  public positions: string;
  public sequenceColumnName: string;
  private _positionColumn?: DG.Column;
  private _boxPlotViewer?: DG.Viewer;
  public valueColumnName: string;
  public showPositionInfo: boolean = true;
  constructor() {
    super();
    this.positions = this.string('positions', '1', {description: 'Comma-separated sequence positions (1-based) to analyze'});
    this.sequenceColumnName = this.column('sequence', {semType: DG.SEMTYPE.MACROMOLECULE, nullable: false});
    this.valueColumnName = this.column('value', {columnTypeFilter: 'numerical', nullable: false});
    this.showPositionInfo = this.bool('showPositionInfo', true, {nullable: false, defaultValue: true, description: 'Show position selector in the viewer header'});
    this.subs.push(grok.events.onContextMenu.subscribe((e) => {
      if (e.causedBy && e.causedBy.target && this._boxPlotViewer?.root.contains(e.causedBy.target)) {
        e.causedBy.preventDefault();
        e.causedBy.stopPropagation();
        e.causedBy.stopImmediatePropagation();
      }
    }));
  }

  getPositionFromColumn(): number {
    // also starting from 1!
    const col = this.dataFrame.col(this.sequenceColumnName);
    if (col == null)
      return 1;
    const positionString = col.getTag(bioTAGS.selectedPosition) ?? '1';
    const position = parseInt(positionString);
    if (Number.isNaN(position))
      return 1;
    return Math.max(1, position);
  }

  parsePositions(): number[] {
    if (!this.positions) return [1];
    const result = this.positions.split(',')
      .map((s) => parseInt(s.trim()))
      .filter((n) => !Number.isNaN(n) && n >= 1);
    return [...new Set(result)].sort((a, b) => a - b);
  }

  private _setPositions(positions: number[]): void {
    const unique = [...new Set(positions)].sort((a, b) => a - b);
    this.getProperty('positions')!.set(this, unique.join(', '));
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
    this.getProperty('positions')!.set(this, String(this.getPositionFromColumn()));

    this.subs.push(DG.debounce(this.dataFrame.onMetadataChanged, 200).subscribe((_) => {
      const curPosition = this.getPositionFromColumn();
      const currentPositions = this.parsePositions();
      if (!currentPositions.includes(curPosition))
        this._setPositions([curPosition]);
    }));
  }

  render(): void {
    const positions = this.parsePositions();
    if (this.dataFrame == null || !this.sequenceColumnName || !positions.length || !this._positionColumn || !this.valueColumnName)
      return;

    $(this.root).empty();
    this._boxPlotViewer?.detach();
    const seqHelper = PeptideUtils.getSeqHelper();
    const sequenceColumn = this.dataFrame.col(this.sequenceColumnName)!;
    const seqHandler = seqHelper.getSeqHandler(sequenceColumn);
    const maxPos = seqHandler.maxLength;

    const originals = positions.map((p) => seqHandler.getMonomersAtPosition(p - 1, false));
    this._positionColumn.init((i) => originals.map((c) => c[i]).join(MONOMER_MOTIF_SPLITTER));

    this._boxPlotViewer = this.dataFrame.plot.box({categoryColumnNames: [this._positionColumn.name], plotStyle: 'violin',
      valueColumnName: this.valueColumnName, colorColumnName: this._positionColumn.name, showColorSelector: false, showSizeSelector: false, showCategorySelector: false,
      legendVisibility: DG.VisibilityMode.Never, markerColorColumnName: this._positionColumn.name, title: 'Sequence Position Statistics',
      autoLayout: false, labelOrientation: 'Vert',
    });

    if (this.showPositionInfo) {
      const selectorDiv = this._renderPositionSelector(positions, maxPos);
      this.root.appendChild(selectorDiv);
    }

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

  private _renderPositionSelector(positions: number[], maxPos: number): HTMLElement {
    const container = ui.divH([], {style: {alignItems: 'center', justifyContent: 'center', gap: '4px', flexWrap: 'wrap', padding: '4px 8px', width: '100%'}});

    const label = ui.label(`${this.sequenceColumnName}:`);
    label.style.fontWeight = 'bold';
    label.style.marginRight = '4px';
    container.appendChild(label);

    for (const pos of positions)
      container.appendChild(this._createChip(pos, positions, maxPos));

    container.appendChild(this._createAddButton(container, positions, maxPos));
    return container;
  }

  private _createChip(pos: number, allPositions: number[], maxPos: number): HTMLElement {
    const isSingle = allPositions.length === 1;
    const posLabel = document.createElement('span');
    posLabel.textContent = String(pos);
    posLabel.style.marginRight = '3px';

    const actionIcon = isSingle ?
      ui.iconFA('pencil', () => this._showInlineEdit(chip, pos, maxPos), 'Edit position') :
      ui.iconFA('times', () => this._setPositions(allPositions.filter((p) => p !== pos)), 'Remove position');
    actionIcon.style.cursor = 'pointer';
    actionIcon.style.fontSize = '10px';
    actionIcon.style.opacity = '0.6';

    const chip = ui.divH([posLabel, actionIcon], {style: {
      display: 'inline-flex', alignItems: 'center',
      background: '#f0f0f0', borderRadius: '12px', padding: '2px 8px',
      fontSize: '12px', cursor: 'default', border: '1px solid #d0d0d0',
    }});
    return chip;
  }

  private _showInlineEdit(chip: HTMLElement, currentPos: number, maxPos: number): void {
    const input = document.createElement('input');
    input.type = 'number';
    input.min = '1';
    input.max = String(maxPos);
    input.value = String(currentPos);
    input.style.cssText = 'width:45px;height:22px;font-size:12px;text-align:center;border-radius:12px;border:1px solid #aaa;outline:none;padding:0 4px;';

    let handled = false;
    const confirm = () => {
      if (handled) return;
      handled = true;
      const val = parseInt(input.value);
      const positions = this.parsePositions();
      if (!Number.isNaN(val) && val >= 1 && val <= maxPos && val !== currentPos)
        this._setPositions(positions.map((p) => p === currentPos ? val : p));
      else
        input.replaceWith(chip);
    };

    input.addEventListener('keydown', (e) => {
      if (e.key === 'Enter') confirm();
      if (e.key === 'Escape') {handled = true; chip.parentElement?.replaceChild(chip, input);}
    });
    input.addEventListener('blur', confirm);

    chip.parentElement?.replaceChild(input, chip);
    input.focus();
    input.select();
  }

  private _createAddButton(container: HTMLElement, positions: number[], maxPos: number): HTMLElement {
    const addBtn = ui.iconFA('plus-circle', () => {
      if (container.querySelector('.position-add-input'))
        return;
      const input = document.createElement('input');
      input.type = 'number';
      input.min = '1';
      input.max = String(maxPos);
      input.className = 'position-add-input';
      input.style.cssText = 'width:45px;height:22px;font-size:12px;text-align:center;border-radius:12px;border:1px solid #aaa;outline:none;padding:0 4px;';

      let handled = false;
      const confirm = () => {
        if (handled) return;
        handled = true;
        const val = parseInt(input.value);
        if (!Number.isNaN(val) && val >= 1 && val <= maxPos && !positions.includes(val))
          this._setPositions([...positions, val]);
        else
          input.remove();
      };

      input.addEventListener('keydown', (e) => {
        if (e.key === 'Enter') confirm();
        if (e.key === 'Escape') {handled = true; input.remove();}
      });
      input.addEventListener('blur', confirm);

      container.insertBefore(input, addBtn);
      input.focus();
    }, 'Add position');
    addBtn.style.cursor = 'pointer';
    addBtn.style.fontSize = '14px';
    addBtn.style.color = '#2083d5';
    return addBtn;
  }

  onPropertyChanged(property: DG.Property | null): void {
    if (this.dataFrame == null || !this.sequenceColumnName)
      return;

    this.render();
  }
}
