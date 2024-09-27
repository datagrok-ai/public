import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {Unsubscribable} from 'rxjs';
import {PolyToolPlaceholders, PolyToolPlaceholdersBreadth} from './types';
import {parseMonomerSymbolList} from './pt-placeholders-input';

export class PolyToolPlaceholdersBreadthInput extends DG.JsInputBase<DG.DataFrame> {
  get inputType(): string { return 'Breadth'; }

  get dataType(): string { return DG.TYPE.DATA_FRAME; }

  getInput(): HTMLElement { return this.gridHost; }

  getValue(): DG.DataFrame { return this.grid.dataFrame; }

  setValue(value: DG.DataFrame): void { this.grid.dataFrame = value; }

  getStringValue(): string { return this.grid.dataFrame.toCsv(); }

  setStringValue(str: string): void { this.grid.dataFrame = DG.DataFrame.fromCsv(str); }

  get placeholdersBreadthValue(): PolyToolPlaceholdersBreadth {
    return dfToPlaceholdersBreadth(this.grid.dataFrame);
  }

  private readonly gridHost: HTMLDivElement;
  public readonly grid: DG.Grid;

  private subs: Unsubscribable[] = [];

  protected constructor(name: string | undefined, grid: DG.Grid, heightRowCount?: number) {
    super();

    if (name) this.captionLabel.innerText = name;

    this.gridHost = ui.div([], {
      classes: 'ui-input-editor',
      style: {width: '100%', height: '100%', marginTop: '-8px', marginBottom: '8px', paddingBottom: '4px'},
    });

    this.grid = grid;
    this.gridHost.append(this.grid.root);

    if (heightRowCount != null) {
      this.updateGridHeight(heightRowCount + 0.7);
    } else {
      this.updateGridHeight(this.grid.dataFrame.rowCount + 0.6);
      this.subs.push(this.grid.dataFrame.onRowsAdded
        .subscribe(() => { this.updateGridHeight(this.grid.dataFrame.rowCount + 0.6); }));
    }
    this.grid.root.style.width = `100%`;

    this.subs.push(this.grid.dataFrame.onDataChanged.subscribe(() => {
      this.fireChanged();
    }));

    this.subs.push(ui.onSizeChanged(this.grid.root).subscribe(() => {
      this.grid.columns.byIndex(3)!.width = this.grid.root.clientWidth - this.grid.horzScroll.root.offsetWidth -
        this.grid.columns.byIndex(0)!.width - this.grid.columns.byIndex(1)!.width - this.grid.columns.byIndex(2)!.width - 10;
    }));

    this.root.classList.add('ui-input-polytool-pos-grid');
    this.root.append(this.gridHost);
  }

  detach(): void {
    for (const sub of this.subs) sub.unsubscribe();
  }

  public static async create(
    name?: string, options?: {}, heightRowCount?: number
  ): Promise<PolyToolPlaceholdersBreadthInput> {
    const df: DG.DataFrame = DG.DataFrame.fromColumns([
      DG.Column.fromType(DG.COLUMN_TYPE.INT, 'Start', 0),
      DG.Column.fromType(DG.COLUMN_TYPE.INT, 'End', 0),
      DG.Column.fromType(DG.COLUMN_TYPE.STRING, 'Monomers', 0),
    ])!;
    const grid = (await df.plot.fromType(DG.VIEWER.GRID, options)) as DG.Grid;
    grid.sort(['Start', 'End']);
    return new PolyToolPlaceholdersBreadthInput(name, grid, heightRowCount);
  }

  // -- Update view --

  private updateGridHeight(visibleRowCount: number): void {
    const gridHeight = this.grid.colHeaderHeight + visibleRowCount * this.grid.props.rowHeight + 6 + 2;
    this.grid.root.style.height = `${gridHeight}px`;
  }

  // -- Handle events --

  private gridRootOnSizeChanged(): void {
    this.grid.columns.byIndex(3)!.width = this.grid.root.clientWidth - this.grid.horzScroll.root.offsetWidth -
      this.grid.columns.byIndex(0)!.width - this.grid.columns.byIndex(1)!.width - this.grid.columns.byIndex(2)!.width - 10;
  }
}

export function dfToPlaceholdersBreadth(df: DG.DataFrame): PolyToolPlaceholdersBreadth {
  const res: PolyToolPlaceholdersBreadth = [];
  for (let rowI = 0; rowI < df.rowCount; rowI++) {
    const startPos = parseInt(df.get('Start', rowI)) - 1;
    const endPos = parseInt(df.get('End', rowI)) - 1;
    if (!isNaN(startPos) && !isNaN(endPos)) {
      const monomerSymbolList = parseMonomerSymbolList(df.get('Monomers', rowI));
      res.push({start: startPos, end: endPos, monomers: monomerSymbolList});
    }
  }
  return res;
}
