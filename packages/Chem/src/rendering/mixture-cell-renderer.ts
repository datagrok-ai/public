import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {RDKitCellRenderer} from './rdkit-cell-renderer';
import {tableFromMap} from 'datagrok-api/ui';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {addStructureFields, Mixfile, MixfileComponent, STRUCTURE_FIELDS} from '../utils/mixfile';
import {trimText} from '@datagrok-libraries/gridext/src/utils/TextUtils';

export interface Component {
  x: number;
  y: number;
  w: number;
  h: number;
  comp: MixfileComponent;
}

export class MixtureCellRenderer extends RDKitCellRenderer {
  private _componentsCache: DG.LruCache<string, Component[]> = new DG.LruCache<string, Component[]>(100);
  private hoveredRowIdx: number | undefined = undefined;
  private hoveredCompIdx: number | undefined = undefined;

  get name() {return 'Mixture cell renderer';}
  get cellType() {return 'ChemicalMixture';}

  constructor(rdKitModule: RDModule) {
    super(rdKitModule);
  }

  private _cellKey(gridCell: DG.GridCell): string {
    // Use grid root id (if available), column name, and row index for uniqueness
    const gridId = (gridCell.grid.root && gridCell.grid.root.id) ? gridCell.grid.root.id : gridCell.grid.toString();
    const colName = gridCell.cell.column?.name ?? '';
    const rowIdx = gridCell.tableRowIndex;
    return `${gridId}:${colName}:${rowIdx}`;
  }

  // Helper to draw brackets around a region
  private drawBrackets(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number) {
    x = x + 5;
    y = y + 5;
    w = w - 10;
    h = h - 10;
    const bracketLen = Math.min(10, h / 6);
    g.save();
    g.setTransform(1, 0, 0, 1, 0, 0);
    g.strokeStyle = '#222';
    g.lineWidth = 0.5;
    // Left bracket
    g.beginPath();
    g.moveTo(x, y);
    g.lineTo(x, y + h);
    g.moveTo(x, y);
    g.lineTo(x + bracketLen, y);
    g.moveTo(x, y + h);
    g.lineTo(x + bracketLen, y + h);
    g.stroke();
    // Right bracket
    g.beginPath();
    g.moveTo(x + w, y);
    g.lineTo(x + w, y + h);
    g.moveTo(x + w, y);
    g.lineTo(x + w - bracketLen, y);
    g.moveTo(x + w, y + h);
    g.lineTo(x + w - bracketLen, y + h);
    g.stroke();
    g.restore();
  }

  private createAnnotation(
    comp: MixfileComponent,
    g: CanvasRenderingContext2D,
    curX: number,
    compW: number,
    compH: number,
    minY: number,
    annotationFrac: number,
    molH?: number,
  ): void {
    let annotation = '';
    if (comp.relation)
      annotation = `${comp.relation}`;
    if (comp.quantity && comp.units)
      annotation += (annotation ? ' ' : '') + `${comp.quantity} ${comp.units}`;
    if (comp.ratio && comp.ratio.length > 1)
      annotation += (annotation ? ', ' : '') + `ratio ${comp.ratio[0]}/${comp.ratio[1]}`;
    if (comp.error)
      annotation += (annotation ? ', ' : '') + `SE: ${comp.error}`;
    if (annotation) {
      g.save();
      g.setTransform(1, 0, 0, 1, 0, 0);
      g.font = '9px sans-serif';
      g.fillStyle = '#222';
      g.textAlign = 'center';
      // Trim annotation to fit width
      const trimmed = trimText(annotation, g, compW - 4);
      g.fillText(trimmed, curX + compW / 2, minY + (molH ?? compH) + compH * annotationFrac * 0.7);
      g.restore();
    }
  }

  // Helper to count all leaf components in a mixture tree
  private countLeafComponents(mixture: MixfileComponent): number {
    if (!mixture.contents || !Array.isArray(mixture.contents) || mixture.contents.length === 0)
      return 1;
    return mixture.contents.map((c) => this.countLeafComponents(c)).reduce((a, b) => a + b, 0);
  }

  // Helper to draw a single component (molfile > smiles > name as text)
  private drawComponent(
    g: CanvasRenderingContext2D,
    x: number,
    y: number,
    w: number,
    h: number,
    comp: MixfileComponent,
    cellStyle: DG.GridCellStyle,
    annotationFrac: number,
    highlight?: boolean,
  ): void {
    const molH = h * (1 - annotationFrac);
    if (comp.molfile || comp.smiles) {
      this._drawMolecule(x, y, w, molH, g.canvas, comp.molfile ?? comp.smiles!,
        [], false, false, cellStyle, false, {}, undefined);
    } else {
      // Draw the name as text (centered)
      g.save();
      g.setTransform(1, 0, 0, 1, 0, 0);
      g.font = '12px sans-serif';
      g.fillStyle = '#444';
      g.textAlign = 'center';
      g.textBaseline = 'middle';
      const text = comp.name || 'Component';
      // Trim name to fit width
      const trimmed = trimText(text, g, w - 4);
      g.fillText(trimmed, x + w / 2, y + molH / 2);
      g.restore();
    }
    if (highlight === true) {
      g.save();
      g.setTransform(1, 0, 0, 1, 0, 0);
      g.strokeStyle = '#1976d2';
      g.lineWidth = 0.5;
      g.strokeRect(x + 3, y + 3, w - 10, molH - 10);
      g.restore();
    }
  }

  // Recursive mixture/component drawing, with equal width for all leaf components
  private drawMixture(
    g: CanvasRenderingContext2D,
    x: number,
    y: number,
    w: number,
    h: number,
    mixture: MixfileComponent,
    cellStyle: DG.GridCellStyle,
    boxes: Component[],
    leafWidth?: number,
    hoveredIdx?: number,
  ): void {
    const margin = 4;
    const annotationFrac = 0.1; // 10% for annotation, 90% for molecule
    const molFrac = 1 - annotationFrac;
    if (mixture.contents && Array.isArray(mixture.contents) && mixture.contents.length > 0) {
      // Count all leaf components in this subtree
      const totalLeaves = this.countLeafComponents(mixture);
      const compH = h - 2 * margin;
      let curX = x + margin;
      const remainingW = w - 2 * margin;
      for (const comp of mixture.contents) {
        const compLeaves = this.countLeafComponents(comp);
        const compW = (compLeaves / totalLeaves) * remainingW;
        if (comp && comp.contents && Array.isArray(comp.contents) && comp.contents.length > 0) {
          this.drawMixture(g, curX, y + margin, compW, compH, comp, cellStyle, boxes, leafWidth, hoveredIdx);
          this.drawBrackets(g, curX, y + margin, compW, compH);
          this.createAnnotation(comp, g, curX, compW, compH, y + margin, annotationFrac);
          boxes.push({x: curX, y: y + margin, w: compW, h: compH, comp});
        } else {
          // Single component
          const idx = boxes.length;
          const isHighlighted = hoveredIdx !== undefined && hoveredIdx === idx;
          this.drawComponent(g, curX, y + margin, compW, compH, comp, cellStyle, annotationFrac, isHighlighted);
          this.createAnnotation(comp, g, curX, compW, compH, y + margin, annotationFrac, compH * molFrac);
          boxes.push({x: curX, y: y + margin, w: compW, h: compH * molFrac, comp});
        }
        curX += compW;
      }
    } else {
      // Single component
      const idx = boxes.length;
      const isHighlighted = hoveredIdx !== undefined && hoveredIdx === idx;
      this.drawComponent(g, x, y, w, h, mixture, cellStyle, annotationFrac, isHighlighted);
      this.createAnnotation(mixture, g, x, w, h, y, annotationFrac, h * molFrac);
      boxes.push({x: x, y: y, w: w, h: h * molFrac, comp: mixture});
    }
  }

  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle): void {
    const value = gridCell.cell.value;
    if (!value) return;

    const r = window.devicePixelRatio;
    x = r * x;
    y = r * y;
    w = r * w;
    h = r * h;
    const mixture: Mixfile = typeof value === 'string' ? JSON.parse(value) : value;
    if (!mixture.contents || !Array.isArray(mixture.contents)) return;
    const boxes: Component[] = [];
    this.drawMixture(g, x, y, w, h, mixture, cellStyle, boxes, undefined,
      this.hoveredRowIdx === gridCell.tableRowIndex ? this.hoveredCompIdx : undefined);
    // Store in LRU cache
    this._componentsCache.set(this._cellKey(gridCell), boxes);
  }

  onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    this.hitTest((comp: MixfileComponent, x: number, y: number, idx: number) => {
      this.createTooltip(comp, x, y);
      const rowIdx = gridCell.tableRowIndex;
      let changed = false;
      if (rowIdx != null && rowIdx !== -1) {
        if (this.hoveredRowIdx !== rowIdx || this.hoveredCompIdx !== idx) {
          this.hoveredRowIdx = rowIdx;
          this.hoveredCompIdx = idx;
          changed = true;
        }
      } else {
        if (this.hoveredRowIdx !== undefined || this.hoveredCompIdx !== undefined) {
          this.hoveredRowIdx = undefined;
          this.hoveredCompIdx = undefined;
          changed = true;
        }
      }
      //invalidate only in case hovered component has changed
      if (changed)
        gridCell.grid.invalidate();
    }, e, gridCell, true);
  }

  onMouseLeave(gridCell: DG.GridCell, _e: MouseEvent): void {
    //invalidating only in case we were hovering over a component before leaving cell
    if (this.hoveredRowIdx !== undefined || this.hoveredCompIdx !== undefined) {
      this.hoveredRowIdx = undefined;
      this.hoveredCompIdx = undefined;
      gridCell.grid.invalidate();
    }
  }

  onClick(gridCell: DG.GridCell, e: MouseEvent): void {
    this.hitTest((comp: MixfileComponent) => {
      const struct = comp.molfile ?? comp.smiles;
      if (struct)
        grok.shell.o = DG.SemanticValue.fromValueType(struct, DG.SEMTYPE.MOLECULE);
    }, e, gridCell);
  }

  hitTest(action: Function, e: MouseEvent, gridCell: DG.GridCell, mouseMove?: boolean) {
    const boxes = this._componentsCache.get(this._cellKey(gridCell));
    if (!boxes) return;
    const rect = (e.target as HTMLElement).getBoundingClientRect();
    // Scale mouse coordinates to match the scaled box coordinates
    const devicePixelRatio = window.devicePixelRatio || 1;
    const x = (e.clientX - rect.left) * devicePixelRatio;
    const y = (e.clientY - rect.top) * devicePixelRatio;
    let overComponent = false;
    let isStructure = false;
    for (let i = 0; i < boxes.length; ++i) {
      const box = boxes[i];
      // Calculate center and 1/3 region
      const cx = box.x + box.w / 2;
      const cy = box.y + box.h / 2;
      const dx = Math.abs(x - cx);
      const dy = Math.abs(y - cy);
      const wx = box.w / 3;
      const hy = box.h / 3;
      if (dx <= wx && dy <= hy) {
        overComponent = true;
        isStructure = !!(box.comp.molfile || box.comp.smiles);
        action(box.comp, e.clientX, e.clientY, i);
        break;
      }
    }
    if (!mouseMove)
      return;
    if (!overComponent) {
      if (this.hoveredRowIdx !== undefined || this.hoveredCompIdx !== undefined) {
        this.hoveredRowIdx = undefined;
        this.hoveredCompIdx = undefined;
        ui.tooltip.hide();
        gridCell.grid.invalidate();
      }
    }
    // Set pointer cursor if over a box, otherwise default
    const gridRoot = gridCell.grid.root as HTMLElement;
    if (gridRoot)
      gridRoot.style.cursor = overComponent && isStructure ? 'pointer' : '';
  }

  createTooltip(comp: MixfileComponent, x: number, y: number) {
    const propsToView: {[key: string]: any} = {};
    Object.keys(comp).forEach((key) => {
      if (!STRUCTURE_FIELDS.includes(key) && key !== 'contents')
        propsToView[key] = (comp as any)[key];
    });
    addStructureFields(propsToView, comp);
    if (Object.keys(propsToView).length)
      ui.tooltip.show(tableFromMap(propsToView), x, y);
  }
}
