import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {PeptiHitTemplate} from './types';
import {PeptiHitHelmColName, TileCategoriesColName, ViDColName} from './consts';
import {getNewVid} from './utils/calculate-single-cell';
import '../../css/hit-triage.css';
import {_package} from '../package';
import {Subscription} from 'rxjs';
import {filter} from 'rxjs/operators';
import {HitDesignApp} from './hit-design-app';
import {PeptiHitInfoView} from './pepti-hits-views/info-view';

export class PeptiHitApp extends HitDesignApp<PeptiHitTemplate> {
  _helmColName: string = PeptiHitHelmColName;
  constructor(c: DG.FuncCall) {
    super(c, 'PeptiHit', (app) => new PeptiHitInfoView(app as PeptiHitApp));
  }

  public async setTemplate(template: PeptiHitTemplate, campaignId?: string): Promise<void> {
    await super.setTemplate(template, campaignId);
    this._helmColName = this.dataFrame!.columns.bySemType(DG.SEMTYPE.MACROMOLECULE)?.name ?? PeptiHitHelmColName;
  }

  get helmColName(): string {
    return this._helmColName ??
        this.dataFrame!.columns.bySemType(DG.SEMTYPE.MACROMOLECULE)?.name;
  }

  async performSingleCellCalculations(cellIndex: number, cellValue?: string) {
    try {
      if (!cellValue)
        throw new Error('No cell value provided');
      const col = DG.Column.fromStrings(this.helmColName, [cellValue]);
      const table = DG.DataFrame.fromColumns([col]);
      col.semType = DG.SEMTYPE.MACROMOLECULE;
      col.setTag('units', 'helm');
      col.setTag('.alphabetIsMultichar', 'true');
      await grok.functions.call('Bio:toAtomicLevel', {table: table, seqCol: col, nonlinear: true});
      const molCol = table.columns.names().find((c) => c.toLowerCase() !== this.helmColName.toLowerCase());
      const mol = molCol ? table.col(molCol)?.get(0) : null;
      if (!mol)
        throw new Error('Failed to convert sequence to atomic level');
      this.dataFrame!.set(this.molColName, cellIndex, mol);
      await super.performSingleCellCalculations(cellIndex, mol);
    } catch (e) {
      console.error(e);
    }
  }

  protected getDesignView(): DG.TableView {
    const subs: Subscription[] = [];
    const isNew = this.dataFrame!.col(this.helmColName)?.toList().every((m) => !m && m === '');
    const helmCol = this.dataFrame!.col(this.helmColName);
    if (!helmCol)
      throw new Error('No helm column found');
    helmCol.semType = DG.SEMTYPE.MACROMOLECULE;
    helmCol.setTag('units', 'helm');
    helmCol.setTag('.alphabetIsMultichar', 'true');
    helmCol.setTag('cell.renderer', 'helm');

    const view = grok.shell.addTableView(this.dataFrame!);
    this._designViewName = this.campaign?.name ?? this._designViewName;
    view.name = this._designViewName;
    view._onAdded();
    const layoutViewState = this._campaign?.layout ?? this.template?.layoutViewState;
    if (layoutViewState) {
      try {
        const layout = DG.ViewLayout.fromViewState(layoutViewState);
        view.loadLayout(layout);
      } catch (e) {
        grok.shell.error('Failed to apply layout. Falling back to default layout.');
        console.error(e);
      }
    }

    if (isNew)
      grok.functions.call('Helm:editMoleculeCell', {cell: view.grid.cell(this.helmColName, 0)});

    subs.push(this.dataFrame!.onRowsAdded.pipe(filter(() => !this.isJoining))
      .subscribe(() => { // TODO, insertion of rows in the middle
        try {
          if (this.template!.stages?.length > 0) {
            for (let i = 0; i < this.dataFrame!.rowCount; i++) {
              const colVal = this.dataFrame!.col(TileCategoriesColName)!.get(i);
              if (!colVal || colVal === '' || this.dataFrame!.col(TileCategoriesColName)?.isNone(i))
                this.dataFrame!.set(TileCategoriesColName, i, this.template!.stages[0]);
            }
          }
          let lastAddedCell: DG.GridCell | null = null;
          for (let i = 0; i < this.dataFrame!.rowCount; i++) {
            const cell = view.grid.cell(this.helmColName, i);
            if (!cell)
              continue;
            if (cell.cell.value === '' || cell.cell.value === null)
              lastAddedCell = cell;
          }
          if (lastAddedCell)
            grok.functions.call('Helm:editMoleculeCell', {cell: lastAddedCell});
        } catch (e) {
          console.error(e);
        }
      }));

    subs.push(grok.events.onContextMenu.subscribe((args) => {
      try {
        const viewer: DG.Viewer = args?.args?.context;
        if (!viewer)
          return;
        if (viewer?.type !== DG.VIEWER.GRID)
          return;
        if (!viewer.tableView || viewer.tableView.id !== view.id)
          return;
        if (args?.args?.item?.tableColumn?.name !== this.helmColName || !args?.args?.item?.isTableCell)
          return;
        const menu: DG.Menu = args?.args?.menu;
        if (!menu)
          return;
        menu.item('Add new row', () => {
              this.dataFrame!.rows.addNew(null, true);
        });
        menu.item('Duplicate Sequence', () => {
          try {
            const row = this.dataFrame!.rows.addNew(null, true);
            const cell = row.get(this.helmColName);
            if (cell != null && row.idx > -1)
              this.dataFrame!.cell(row.idx, this.helmColName).value = args?.args?.item?.cell?.value ?? '';
          } catch (e) {
            console.error(e);
          }
        });

        const cellIndex = args?.args?.item?.tableRowIndex;
        const cellValue = args?.args?.item?.cell?.value;
        if (cellValue && (cellIndex ?? -1) > -1) {
          menu.item('Re-Run Calculations', async () => {
            try {
              await this.performSingleCellCalculations(cellIndex, cellValue);
            } catch (e) {
              console.error(e);
            }
          });
        }
      } catch (e: any) {
        grok.log.error(e);
      }
    }));

    if (!view?.grid) {
      grok.shell.error('Applied layout created view without grid. Resetting layout.');
      view.resetLayout();
    }

    view?.grid && subs.push(view.grid.onCellValueEdited.subscribe(async (gc) => {
      try {
        if (gc.tableColumn?.name === TileCategoriesColName) {
          await this.saveCampaign(undefined, false);
          return;
        }
        if (gc.tableColumn?.name !== this.helmColName)
          return;
        const newValue = gc.cell.value;
        const newValueIdx = gc.tableRowIndex!;
        let newVid = this.dataFrame!.col(ViDColName)?.get(newValueIdx);
        let foundMatch = false;
        // try to find existing sequence
        if (newValue) {
          try {
            const canonicals = gc.tableColumn.toList();
            const canonicalNewValue = newValue;
            if (canonicals?.length === this.dataFrame!.rowCount) {
              for (let i = 0; i < canonicals.length; i++) {
                if (canonicals[i] === canonicalNewValue &&
                        i !== newValueIdx && this.dataFrame!.col(ViDColName)?.get(i)) {
                  newVid = this.dataFrame!.col(ViDColName)?.get(i);
                  foundMatch = true;
                  break;
                }
              }
            }
          } catch (e) {
            console.error(e);
          }
        }
        // if the vid was duplicated, generate a new one
        if (this.duplicateVidCache && !foundMatch &&
              this.duplicateVidCache.valueCounts[this.duplicateVidCache.indexes[newValueIdx]] > 1)
          newVid = null;

        if (!newVid || newVid === '')
          newVid = getNewVid(this.dataFrame!.col(ViDColName)!);

        this.dataFrame!.col(ViDColName)!.set(newValueIdx, newVid, false);

        this.performSingleCellCalculations(newValueIdx, newValue);
      } catch (e) {
        console.error(e);
      }
    }));

    view?.grid && subs.push(view.grid.onCellRender.subscribe((args) => {
      try {
        // color duplicate vid values
        const cell = args.cell;
        if (!cell || !cell.isTableCell || !cell.tableColumn || !this.duplicateVidCache ||
              cell.tableColumn.name !== ViDColName || (cell.tableRowIndex ?? -1) < 0)
          return;

        if (this.duplicateVidCache.valueCounts[this.duplicateVidCache.indexes[cell.tableRowIndex!]] > 1) {
          args.cell.style.backColor =
                DG.Color.setAlpha(DG.Color.getCategoricalColor(this.duplicateVidCache.indexes[cell.tableRowIndex!])
                  , 150);
        }
      } catch (e) {}
    }));

    this.initDesignViewRibbons(view, subs);
    view.parentCall = this.parentCall;
    return view;
  }
}
