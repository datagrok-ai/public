import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {fromEvent, Unsubscribable} from 'rxjs';

import {PolyToolPlaceholder} from './types';
import {GridCellRenderArgs} from 'datagrok-api/dg';

export class PolyToolPlaceholdersInput extends DG.JsInputBase<DG.DataFrame> {
  get inputType(): string { return 'Positions'; }

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

  get placeholdersValue(): PolyToolPlaceholder[] {
    return dfToPlaceholders(this.grid.dataFrame);
  }

  private readonly gridHost: HTMLDivElement;
  private grid!: DG.Grid;

  private subs: Unsubscribable[] = [];

  protected constructor(
    private readonly name: string | undefined,
    heightRowCount?: number, options?: {},
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
    // this.captionLabel.appendChild(ui.divH([ui.divText(this.name ?? this.caption), ui.icons.add(() => {}, 'Add')]));
    let removeCol: DG.Column<boolean>;
    const df: DG.DataFrame = DG.DataFrame.fromColumns([
      removeCol = DG.Column.fromType(DG.COLUMN_TYPE.STRING, 'Remove', 0),
      DG.Column.fromType(DG.COLUMN_TYPE.INT, 'Position', 0),
      DG.Column.fromType(DG.COLUMN_TYPE.STRING, 'Monomers', 0),])!;
    removeCol.setTag(DG.TAGS.FRIENDLY_NAME, '');
    this.grid = (await df.plot.fromType(DG.VIEWER.GRID, options)) as DG.Grid;
    this.grid.sort(['Position']);
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

  /** @param {number} posOutIdx continuous (outer) position index, 0-based */
  addPosition(posOutIdx: number): void {
    const phDf = this.grid.dataFrame;
    const posList = phDf.columns.byName('Position').toList();
    let rowIdx = posList.indexOf(posOutIdx + 1);
    if (rowIdx === -1) {
      rowIdx = posList.findIndex((v) => isNaN(v));
      if (rowIdx === -1)
        rowIdx = phDf.rows.addNew(['', posOutIdx + 1, '']).idx;
      const tgtCell = this.grid.cell('Monomers', rowIdx);
    }
    phDf.currentCell = phDf.cell(rowIdx, 'Monomers');
    //const gridRowIdx = inputs.placeholders.grid.tableRowToGrid(rowIdx);
    //const monomersGCell = inputs.placeholders.grid.cell('Monomers', gridRowIdx);
    const k = 42;
  }

  public static async create(
    name?: string, options?: {}, heightRowCount?: number
  ): Promise<PolyToolPlaceholdersInput> {
    return new PolyToolPlaceholdersInput(name, heightRowCount, options);
  }

  // -- Update view --

  private updateGridHeight(visibleRowCount: number): void {
    const gridHeight = this.grid.colHeaderHeight + visibleRowCount * this.grid.props.rowHeight + 6 + 2;
    this.grid.root.style.height = `${gridHeight}px`;
  }

  // -- Handle events --

  private gridRootOnSizeChanged(): void {
    this.grid.columns.byIndex(3)!.width = this.grid.root.clientWidth - this.grid.horzScroll.root.offsetWidth -
      this.grid.columns.byIndex(0)!.width - this.grid.columns.byIndex(1)!.width -
      this.grid.columns.byIndex(2)!.width - 10;
  }
}

export function getPlaceholdersFromText(src: string): PolyToolPlaceholder[] {
  const res: PolyToolPlaceholder[] = [];
  for (const line of src.split('\n')) {
    const lineM = /^\s*(?<pos>\d+)\s*:\s*(?<monomers>.+)$/.exec(line);
    if (lineM) {
      const pos: number = parseInt(lineM.groups!['pos']) - 1;
      const monomerList: string[] = lineM.groups!['monomers'].split(',').map((m) => m.trim());
      res.push({position: pos, monomers: monomerList});
    }
  }
  return res;
}

export function dfToPlaceholders(df: DG.DataFrame): PolyToolPlaceholder[] {
  const res: PolyToolPlaceholder[] = [];
  for (let rowI = 0; rowI < df.rowCount; rowI++) {
    if (df.getCol('Position').isNone(rowI))
      continue;
    const pos = parseInt(df.get('Position', rowI)) - 1;
    if (!isNaN(pos)) {
      const monomerSymbolList = parseMonomerSymbolList(df.get('Monomers', rowI));
      res.push({position: pos, monomers: monomerSymbolList});
    }
  }
  return res;
}

export function parseMonomerSymbolList(src: string): string[] {
  // L, L-hArg(Et,Et), "hArg(Et,Et)"
  return src.split(/,(?![^(]*\))/)
    .map((s) => {
      s = s.trim();
      if (s.slice(0, 1) === `"` && s.slice(-1) === `"`) s = s.slice(1, -1);
      if (s.slice(0, 1) === `'` && s.slice(-1) === `'`) s = s.slice(1, -1);
      return s.trim();
    })
    .filter((s) => !!s);
}
