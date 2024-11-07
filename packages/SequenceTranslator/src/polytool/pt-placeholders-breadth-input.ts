import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {fromEvent, Unsubscribable} from 'rxjs';
import {PolyToolBreadthPlaceholder} from './types';
import {parseMonomerSymbolList} from './pt-placeholders-input';
import {GridCellRenderArgs} from 'datagrok-api/dg';

export class PolyToolPlaceholdersBreadthInput extends DG.JsInputBase<DG.DataFrame> {
  get inputType(): string { return 'Breadth'; }

  get dataType(): string { return DG.TYPE.DATA_FRAME; }

  getInput(): HTMLElement { return this.gridHost; }

  getValue(): DG.DataFrame { return this.grid.dataFrame; }

  setValue(value: DG.DataFrame): void {
    this.setDataFrame(value);
  }

  getStringValue(): string { return this.grid.dataFrame.toCsv(); }

  setStringValue(str: string): void {
    this.setDataFrame(DG.DataFrame.fromCsv(str));
  }

  get placeholdersBreadthValue(): PolyToolBreadthPlaceholder[] {
    return dfToPlaceholdersBreadth(this.grid.dataFrame);
  }

  private readonly gridHost: HTMLDivElement;
  public grid: DG.Grid;

  private subs: Unsubscribable[] = [];

  protected constructor(
    private readonly name: string | undefined,
    heightRowCount?: number, options?: {}
  ) {
    super();

    this.caption = name ?? '';
    this.root.classList.add('ui-input-polytool-pos-grid');
    this.root.append(this.gridHost = ui.div([], {
      classes: 'ui-input-editor',
      style: {width: '100%', height: '100%', marginTop: '-8px', marginBottom: '8px', paddingBottom: '4px'},
    }));
    this.render(heightRowCount, options).then(() => {});
  }

  detach(): void {
    for (const sub of this.subs) sub.unsubscribe();
    for (const sub of this.dataFrameSubs) sub.unsubscribe();
  }

  async render(heightRowCount?: number, options?: {}): Promise<void> {
    let removeCol: DG.Column<boolean>;
    const df: DG.DataFrame = DG.DataFrame.fromColumns([
      removeCol = DG.Column.fromType(DG.COLUMN_TYPE.STRING, 'Remove', 0),
      DG.Column.fromType(DG.COLUMN_TYPE.INT, 'Start', 0),
      DG.Column.fromType(DG.COLUMN_TYPE.INT, 'End', 0),
      DG.Column.fromType(DG.COLUMN_TYPE.STRING, 'Monomers', 0),
    ])!;
    removeCol.setTag(DG.TAGS.FRIENDLY_NAME, '');
    this.grid = (await df.plot.fromType(DG.VIEWER.GRID, options)) as DG.Grid;
    this.grid.sort(['Start', 'End']);
    this.grid.columns.byIndex(1)!.width = this.grid.props.rowHeight + 2;

    this.grid.onCellRender.subscribe((eventArgs: GridCellRenderArgs) => {
      const gridCell = eventArgs.cell;
      if (gridCell.tableColumn?.name == removeCol.name && (gridCell?.tableRowIndex ?? -1) >= 0) {
        gridCell.element = ui.div(ui.icons.delete(() => {
          this.grid.dataFrame.rows.removeAt(gridCell.tableRowIndex!);
        }, 'Delete'), {
          style: {
            height: `${gridCell.grid.props.rowHeight}px`, color: 'var(--grey-6)',
            display: 'flex', flexDirection: 'row', justifyContent: 'center'
          }
        });
      }
    });

    this.updateGridHeight(heightRowCount ?? this.grid.dataFrame.rowCount + 0.7);
    this.subs.push(ui.onSizeChanged(this.grid.root)
      .subscribe(this.gridRootOnSizeChanged.bind(this)));
    this.subs.push(fromEvent<KeyboardEvent>(this.grid.root, 'keydown')
      .subscribe((e: KeyboardEvent) => {
        if (e.key === 'Enter') e.stopPropagation();
      }));
    this.setDataFrame(df);

    this.grid.root.style.width = `100%`;
    this.gridHost.append(this.grid.root);
  }

  private dataFrameSubs: Unsubscribable[] = [];

  private setDataFrame(dataFrame: DG.DataFrame): void {
    for (const sub of this.dataFrameSubs) sub.unsubscribe();
    this.grid.dataFrame = dataFrame;

    this.dataFrameSubs.push(this.grid.dataFrame.onRowsRemoved
      .subscribe(() => {
        this.updateGridHeight(this.grid.dataFrame.rowCount + 0.7);
        this.fireChanged();
      }));
    this.dataFrameSubs.push(this.grid.dataFrame.onRowsAdded
      .subscribe(() => { this.updateGridHeight(this.grid.dataFrame.rowCount + 0.7); }));
    this.dataFrameSubs.push(this.grid.dataFrame.onDataChanged.subscribe(() => {
      this.fireChanged();
    }));

    this.updateGridHeight(this.grid.dataFrame.rowCount + 0.7);
    this.fireChanged();
  }

  public static async create(
    name?: string, options?: {}, heightRowCount?: number
  ): Promise<PolyToolPlaceholdersBreadthInput> {
    return new PolyToolPlaceholdersBreadthInput(name, heightRowCount, options);
  }

  // -- Update view --

  private updateGridHeight(visibleRowCount: number): void {
    const gridHeight = this.grid.colHeaderHeight + visibleRowCount * this.grid.props.rowHeight + 6 + 2;
    this.grid.root.style.height = `${gridHeight}px`;
  }

  // -- Handle events --

  private gridRootOnSizeChanged(): void {
    this.grid.columns.byIndex(4)!.width = this.grid.root.clientWidth - this.grid.horzScroll.root.offsetWidth -
      this.grid.columns.byIndex(0)!.width - this.grid.columns.byIndex(1)!.width -
      this.grid.columns.byIndex(2)!.width - this.grid.columns.byIndex(3)!.width - 10;
  }
}

export function dfToPlaceholdersBreadth(df: DG.DataFrame): PolyToolBreadthPlaceholder[] {
  const res: PolyToolBreadthPlaceholder[] = [];
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
