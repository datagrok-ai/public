import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {Unsubscribable} from 'rxjs';

import {PolyToolPlaceholders} from './types';

export class PolyToolPlaceholdersInput extends DG.JsInputBase<DG.DataFrame> {
  get inputType(): string { return 'Positions'; }

  get dataType(): string { return DG.TYPE.DATA_FRAME; }

  getInput(): HTMLElement { return this.gridHost; }

  getValue(): DG.DataFrame { return this.grid.dataFrame; }

  setValue(value: DG.DataFrame): void { this.grid.dataFrame = value; }

  getStringValue(): string { return this.grid.dataFrame.toCsv(); }

  setStringValue(str: string): void { this.grid.dataFrame = DG.DataFrame.fromCsv(str); }

  get placeholdersValue(): PolyToolPlaceholders {
    return dfToPlaceholders(this.grid.dataFrame);
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
      this.grid.columns.byIndex(2)!.width = this.grid.root.clientWidth - this.grid.horzScroll.root.offsetWidth -
        this.grid.columns.byIndex(0)!.width - this.grid.columns.byIndex(1)!.width - 10;
    }));

    this.root.classList.add('ui-input-polytool-pos-grid');
    this.root.append(this.gridHost);
  }

  detach(): void {
    for (const sub of this.subs) sub.unsubscribe();
  }

  public static async create(
    name?: string, options?: {}, heightRowCount?: number
  ): Promise<PolyToolPlaceholdersInput> {
    const df: DG.DataFrame = DG.DataFrame.fromColumns([
      DG.Column.fromType(DG.COLUMN_TYPE.INT, 'Position', 0),
      DG.Column.fromType(DG.COLUMN_TYPE.STRING, 'Monomers', 0),])!;
    const grid = (await df.plot.fromType(DG.VIEWER.GRID, options)) as DG.Grid;
    grid.sort(['Position']);
    return new PolyToolPlaceholdersInput(name, grid, heightRowCount);
  }

  // -- Update view --

  private updateGridHeight(visibleRowCount: number): void {
    const gridHeight = this.grid.colHeaderHeight + visibleRowCount * this.grid.props.rowHeight + 6 + 2;
    this.grid.root.style.height = `${gridHeight}px`;
  }

  // -- Handle events --

  private gridRootOnSizeChanged(): void {
    this.grid.columns.byIndex(2)!.width = this.grid.root.clientWidth - this.grid.horzScroll.root.offsetWidth -
      this.grid.columns.byIndex(0)!.width - this.grid.columns.byIndex(1)!.width - 10;
  }
}

export function getPlaceholdersFromText(src: string): PolyToolPlaceholders {
  const res: PolyToolPlaceholders = [];
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

export function dfToPlaceholders(df: DG.DataFrame): PolyToolPlaceholders {
  const res: PolyToolPlaceholders = [];
  for (let rowI = 0; rowI < df.rowCount; rowI++) {
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
