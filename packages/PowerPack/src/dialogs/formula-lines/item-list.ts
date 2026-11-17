import * as DG from 'datagrok-api/dg';

import {EditorItem, isAnnotationRegionType} from './types';

/** Where a viewer or a table keeps its formula lines and annotation regions. */
interface ItemsStorage {
  formulaLines: DG.FormulaLinesHelper,
  annotationRegions: DG.AnnotationRegionsHelper,
}

/**
 * Formula lines and annotation regions of one source (a viewer or a table).
 * A row index addresses them as one list: the formula lines first, then the annotation regions.
 */
export class ItemList {
  public readonly formulaLines: DG.FormulaLine[];
  public readonly annotationRegions: DG.AnnotationRegion[];

  constructor(private storage: ItemsStorage) {
    this.formulaLines = storage.formulaLines.items;
    this.annotationRegions = storage.annotationRegions.items;
  }

  public get length(): number {return this.formulaLines.length + this.annotationRegions.length;}

  public isFormulaLine(rowIdx: number): boolean {return rowIdx < this.formulaLines.length;}

  /** Row of the item at [itemIdx] among the formula lines, or among the annotation regions. */
  public rowIdx(itemIdx: number, isFormulaLine: boolean): number {
    return isFormulaLine ? itemIdx : itemIdx + this.formulaLines.length;
  }

  public at(rowIdx: number): EditorItem | null {
    if (rowIdx < 0)
      return null;
    return this.isFormulaLine(rowIdx) ?
      this.formulaLines[rowIdx] :
      this.annotationRegions[rowIdx - this.formulaLines.length] ?? null;
  }

  /** Puts the item first among the items of its kind and returns its row. */
  public add(item: EditorItem): number {
    if (isAnnotationRegionType(item.type)) {
      this.annotationRegions.unshift(item);
      return this.formulaLines.length;
    }
    this.formulaLines.unshift(item);
    return 0;
  }

  public removeAt(rowIdx: number): void {
    if (this.isFormulaLine(rowIdx))
      this.formulaLines.splice(rowIdx, 1);
    else
      this.annotationRegions.splice(rowIdx - this.formulaLines.length, 1);
  }

  /** Writes the items back to the source. */
  public save(): void {
    this.storage.formulaLines.clear();
    this.storage.formulaLines.addAll(this.formulaLines);

    this.storage.annotationRegions.clear();
    this.storage.annotationRegions.addAll(this.annotationRegions);
  }
}
