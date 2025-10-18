import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {PlateStateManager} from '../../shared/plate-state-manager';
import {PlateFile} from '../../shared/types';
import {Subscription} from 'rxjs';
import './plate-grid-manager.css';

export class PlateGridManager {
  root: HTMLElement;
  grid: DG.Grid;
  private subscriptions: Subscription[] = [];
  private isSelecting: boolean = false;

  constructor(private stateManager: PlateStateManager) {
    this.root = ui.divV([], 'assay-plates--plate-grid-manager');
    this.grid = DG.Grid.create(DG.DataFrame.create());

    this.root.appendChild(this.grid.root);
    this.initGrid();
    this.subscribeToStateChanges();
  }

  private initGrid(): void {
    this.grid.props.allowEdit = false;
    this.grid.props.allowRowSelection = true;
    this.grid.props.showCurrentCellOutline = false;
    this.grid.props.showRowHeader = false;
    this.grid.props.showColumnGridlines = false;
    this.grid.props.showRowGridlines = true;
    this.grid.props.rowHeight = 45;

    this.grid.onCellClick.subscribe((gc: DG.GridCell) => {
      if (!gc || !gc.isTableCell || this.isSelecting) return;

      const tableRowIndex = this.grid.gridRowToTable(gc.gridRow);
      const state = this.stateManager.currentState;

      if (tableRowIndex !== -1 && state && tableRowIndex !== state.activePlateIdx) {
        this.isSelecting = true;
        this.stateManager.selectPlate(tableRowIndex);
        this.isSelecting = false;
      }
    });

    this.grid.onCellRender.subscribe((args: DG.GridCellRenderArgs) => {
      const cell = args.cell;
      const bounds = args.bounds;
      const g = args.g;

      if (cell.isTableCell && (cell.gridColumn.name === '$TemplateName' || cell.gridColumn.name === 'QC')) {
        const value = cell.cell.value;
        g.fillStyle = value ? '#28a745' : '#dc3545';
        const r = Math.min(bounds.width, bounds.height) / 4;
        g.beginPath();
        g.arc(bounds.x + bounds.width / 2, bounds.y + bounds.height / 2, r, 0, 2 * Math.PI);
        g.fill();
        args.preventDefault();
        return;
      }

      if (cell.isColHeader && cell.gridColumn.name === 'Wells Failed') {
        g.save();
        g.strokeStyle = 'var(--grey-5)';
        g.lineWidth = 1.5;
        g.setLineDash([2, 2]);
        g.strokeRect(bounds.x + 8, bounds.y + bounds.height / 2 - 7, 14, 14);
        g.setLineDash([]);
        g.fillStyle = 'var(--grey-6)';
        g.font = '12px Roboto, Roboto Local';
        g.textAlign = 'start';
        g.textBaseline = 'middle';
        g.fillText('Wells Failed', bounds.x + 28, bounds.y + bounds.height / 2);
        g.restore();
        args.preventDefault();
      }
    });
  }

  private renderGrid(): void {
    const state = this.stateManager.currentState;

    if (!state || state.plates.length === 0) {
      this.grid.dataFrame = DG.DataFrame.create();
      return;
    }

    const plateFiles: PlateFile[] = state.plates;

    const barcodes = DG.Column.string('Barcode', plateFiles.length)
      .init((i) => plateFiles[i].plate.barcode ?? `Plate ${i + 1}`);

    const templateMatch = DG.Column.bool('$TemplateName', plateFiles.length)
      .init((i) => i % 2 === 0);

    const qcStatus = DG.Column.bool('QC', plateFiles.length)
      .init((i) => (i + 1) % 3 !== 0);

    const allDynamicKeys = new Set<string>();
    plateFiles.forEach((pf) => {
      pf.commonProperties?.forEach((_: any, key: any) => allDynamicKeys.add(key));
    });

    const dynamicColumns = Array.from(allDynamicKeys).sort().map((key) => {
      const col = DG.Column.string(key, plateFiles.length);
      for (let i = 0; i < plateFiles.length; i++) {
        const value = plateFiles[i].commonProperties?.get(key);
        col.set(i, value !== undefined && value !== null ? String(value) : '');
      }
      return col;
    });

    const experimentalStatus = DG.Column.string('Experimental Status', plateFiles.length)
      .init('Completed');

    const wellsFailed = DG.Column.int('Wells Failed', plateFiles.length)
      .init((i) => i % 3 === 0 ? i % 5 : 0);

    const finalColumns = [
      barcodes,
      templateMatch,
      qcStatus,
      experimentalStatus,
      ...dynamicColumns,
      wellsFailed,
    ];

    this.grid.dataFrame = DG.DataFrame.fromColumns(finalColumns);

    const barcodeCol = this.grid.columns.byName('Barcode');
    if (barcodeCol) barcodeCol.width = 200;

    const templateCol = this.grid.columns.byName('$TemplateName');
    if (templateCol) {
      templateCol.width = 40;
      templateCol.name = 'Template';
    }

    const qcCol = this.grid.columns.byName('QC');
    if (qcCol) qcCol.width = 40;

    if (state.activePlateIdx !== -1 && !this.isSelecting) {
      setTimeout(() => {
        if (this.grid.dataFrame && this.grid.dataFrame.currentRowIdx !== state.activePlateIdx)
          this.grid.dataFrame.currentRowIdx = state.activePlateIdx;
      }, 0);
    }
  }

  private subscribeToStateChanges(): void {
    const sub = this.stateManager.onStateChange$.subscribe((event) => {
      console.log('PlateGridManager: State changed', event);
      this.renderGrid();
    });
    this.subscriptions.push(sub);
  }

  destroy(): void {
    this.subscriptions.forEach((sub) => sub.unsubscribe());
  }
}
