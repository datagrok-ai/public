import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {PlateStateManager} from '../../shared/plate-state-manager';
import {MappingDialogOptions, PlateFile} from '../../shared/types';
import {Subscription} from 'rxjs';

function createPlateGridSkeleton(): HTMLElement {
  const icon = ui.iconFA('table', null, 'Plates Table');
  icon.style.fontSize = '48px';
  icon.style.color = 'var(--grey-3)';

  const message = ui.divText('Import a file to see the list of plates', 'info-message');
  message.style.marginTop = '12px';

  const skeleton = ui.divV([icon, message], 'plate-grid-skeleton');
  skeleton.style.cssText = `
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    width: 100%;
    height: 100%;
    min-height: 150px;
    border: 1px dashed var(--grey-2);
    border-radius: 4px;
    background-color: var(--grey-0);
  `;
  return skeleton;
}

export class PlateGridManager {
  root: HTMLElement;
  grid: DG.Grid;
  private subscriptions: Subscription[] = [];
  private isSelecting: boolean = false;
  private skeletonEl: HTMLElement | null = null;
  private loadingIndicator: HTMLElement | null = null;

  constructor(private stateManager: PlateStateManager) {
    this.root = ui.divV([], 'plate-grid-manager');
    this.grid = DG.Grid.create(DG.DataFrame.create());

    this.buildComponent();
    this.initGrid();
    this.subscribeToStateChanges();
  }

  private buildComponent(): void {
    this.root.style.width = '100%';
    this.root.style.display = 'flex';
    this.root.style.flexDirection = 'column';
    this.grid.root.style.width = '100%';
    this.grid.root.style.flexGrow = '1';
    this.root.appendChild(this.grid.root);
  }

  private initGrid(): void {
    this.grid.props.allowEdit = false;
    this.grid.props.allowRowSelection = true;
    this.grid.props.showCurrentCellOutline = false;
    this.grid.props.showRowHeader = false;
    this.grid.props.showColumnGridlines = false;
    this.grid.props.showRowGridlines = true;
    this.grid.props.rowHeight = 45;
    this.grid.root.style.width = '100%';

    this.grid.onCurrentCellChanged.subscribe((gc: DG.GridCell) => {
      // START OF FIX: Add a null check for the GridCell object
      if (this.isSelecting || !gc || !gc.isTableCell) return;
      // END OF FIX
      const tableRowIndex = this.grid.gridRowToTable(gc.gridRow);
      const state = this.stateManager.currentState;
      if (tableRowIndex !== -1 && state && tableRowIndex !== state.activePlateIdx) {
        this.isSelecting = true;
        this.stateManager.selectPlate(tableRowIndex);
        this.isSelecting = false;
      }
    });


    this.grid.onCellClick.subscribe((gc: DG.GridCell) => {
      // START OF FIX: Add a null check here as well for consistency
      if (!gc || !gc.isTableCell) return;
      // END OF FIX
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
        return;
      }
    });
  }

  private renderGrid(): void {
    const state = this.stateManager.currentState;
    if (!state || state.plates.length === 0) {
      if (!this.skeletonEl) {
        this.skeletonEl = createPlateGridSkeleton();
        this.root.appendChild(this.skeletonEl);
        this.setupDroppableArea(); // <-- ATTACH DROPPABLE LISTENER
      }
      this.skeletonEl.style.display = 'flex';
      this.grid.root.style.display = 'none'; // HIDE the grid
      return;
    }

    if (this.skeletonEl)
      this.skeletonEl.style.display = 'none';


    this.grid.root.style.display = ''; // SHOW the grid by resetting its display style

    const plateFiles: PlateFile[] = state.plates;

    const barcodes = DG.Column.string('Barcode', plateFiles.length)
      .init((i) => plateFiles[i].plate.barcode ?? `Plate ${i + 1}`);

    const templateMatch = DG.Column.bool('$TemplateName', plateFiles.length).init((i) => i % 2 === 0);
    const qcStatus = DG.Column.bool('QC', plateFiles.length).init((i) => (i + 1) % 3 !== 0);

    const allDynamicKeys = new Set<string>();
    plateFiles.forEach((pf) => {
      pf.commonProperties?.forEach((_:any, key:any) => allDynamicKeys.add(key));
    });

    const dynamicColumns = Array.from(allDynamicKeys).sort().map((key) => {
      const col = DG.Column.string(key, plateFiles.length);
      for (let i = 0; i < plateFiles.length; i++) {
        const value = plateFiles[i].commonProperties?.get(key);
        col.set(i, value !== undefined && value !== null ? String(value) : '');
      }
      return col;
    });

    const experimentalStatus = DG.Column.string('Experimental Status', plateFiles.length).init('Completed');
    const wellsFailed = DG.Column.int('Wells Failed', plateFiles.length).init((i) => i % 3 === 0 ? i % 5 : 0);

    const finalColumns = [
      barcodes, templateMatch, qcStatus, experimentalStatus,
      ...dynamicColumns, wellsFailed,
    ];

    this.grid.dataFrame = DG.DataFrame.fromColumns(finalColumns);

    const barcodeCol = this.grid.columns.byName('Barcode');
    if (barcodeCol) barcodeCol.width = 200;
    const templateCol = this.grid.columns.byName('$TemplateName');
    if (templateCol) { templateCol.width = 40; templateCol.name = 'Template'; }
    const qcCol = this.grid.columns.byName('QC');
    if (qcCol) qcCol.width = 40;

    if (state.activePlateIdx !== -1 && !this.isSelecting) {
      setTimeout(() => {
        if (this.grid.dataFrame && this.grid.dataFrame.currentRowIdx !== state.activePlateIdx)
          this.grid.dataFrame.currentRowIdx = state.activePlateIdx;
      }, 0);
    }
  }

  // --- NEW METHODS FOR DROPPABLE FUNCTIONALITY ---

  private setupDroppableArea(): void {
    if (!this.skeletonEl) return;

    // Visual feedback for drag-and-drop
    const highlightColor = 'var(--blue-2)';
    const defaultColor = 'var(--grey-0)';
    const defaultBorder = '1px dashed var(--grey-2)';
    const highlightBorder = '1px solid var(--blue-5)';

    this.skeletonEl.addEventListener('dragenter', () => this.skeletonEl!.style.backgroundColor = highlightColor);
    this.skeletonEl.addEventListener('dragover', (event) => {
      event.preventDefault();
      this.skeletonEl!.style.backgroundColor = highlightColor;
      this.skeletonEl!.style.border = highlightBorder;
    });
    this.skeletonEl.addEventListener('dragleave', () => {
      this.skeletonEl!.style.backgroundColor = defaultColor;
      this.skeletonEl!.style.border = defaultBorder;
    });

    ui.makeDroppable(this.skeletonEl, {
      acceptDrop: (dragObject) => {
        // The dragged object can be a File from the OS or a Datagrok FileInfo
        const file = dragObject instanceof File ? dragObject : (dragObject as any)?.file;
        return !!file && file.name.toLowerCase().endsWith('.csv');
      },
      doDrop: (dragObject) => {
        const file = dragObject instanceof File ? dragObject : (dragObject as any)?.file;
        if (file) {
          this.skeletonEl!.style.backgroundColor = defaultColor;
          this.skeletonEl!.style.border = defaultBorder;
          this.processDroppedFile(file);
        }
      }
    });
  }

  private async processDroppedFile(file: File): Promise<void> {
    this.showLoadingIndicator();
    try {
      const csvString = await new Promise<string>((resolve) => {
        const reader = new FileReader();
        reader.onload = (e) => resolve(e.target!.result as string);
        reader.readAsText(file);
      });
      const df = DG.DataFrame.fromCsv(csvString);
      await this.stateManager.loadDataFrame(df);
    } catch (e: any) {
      grok.shell.error(`Failed to load file: ${e.message}`);
    } finally {
      this.hideLoadingIndicator();
    }
  }

  private showLoadingIndicator(): void {
    if (!this.skeletonEl) return;
    this.loadingIndicator = ui.divV([
      ui.loader(),
      ui.divText('Loading...')
    ], 'loading-indicator');
    this.loadingIndicator.style.cssText = `
      position: absolute;
      top: 0;
      left: 0;
      width: 100%;
      height: 100%;
      background-color: rgba(255, 255, 255, 0.8);
      display: flex;
      flex-direction: column;
      align-items: center;
      justify-content: center;
      z-index: 100;
      border-radius: 4px;
    `;
    this.skeletonEl.appendChild(this.loadingIndicator);
  }

  private hideLoadingIndicator(): void {
    this.loadingIndicator?.remove();
    this.loadingIndicator = null;
  }

  // --- EXISTING METHODS (UNCHANGED) ---

  showMultiPlateMappingDialog(options: MappingDialogOptions): void {
    const dialog = ui.dialog('Apply Field Mappings');
    dialog.root.style.minWidth = '500px';
    const content = ui.divV([], {style: {gap: '20px'}});

    const mappingsSection = ui.divV([], {style: {gap: '8px'}});
    mappingsSection.appendChild(ui.h3('Mappings to Apply'));

    if (options.sourceMappings.size > 0) {
      const mappingsList = ui.divV([], 'mapping-summary-list');
      options.sourceMappings.forEach((oldName, newName) => {
        const row = ui.divH([
          ui.divText(oldName, 'mapping-summary-old-name'),
          ui.iconFA('long-arrow-alt-right', null, 'maps to'),
          ui.divText(newName, 'mapping-summary-new-name'),
          ui.iconFA('times', () => {
            options.onUndo(newName);
            dialog.close();
            grok.shell.info(`Mapping for '${newName}' has been reset. Open "Manage Mappings" again to proceed.`);
          }, 'Undo this mapping definition')
        ], 'mapping-summary-row');
        mappingsList.appendChild(row);
      });
      mappingsSection.appendChild(mappingsList);
    } else {
      mappingsSection.appendChild(
        ui.divText(
          'No mappings defined yet. Use the validation panel to map columns on the active plate.',
          'info-message'
        )
      );
    }

    const platesSection = ui.divV([], {style: {gap: '8px'}});
    platesSection.appendChild(ui.h3('Apply to Plates'));

    const plateCheckboxes: DG.InputBase<boolean>[] = [];
    const plateSelectionHost = ui.divV([], 'plate-selection-host');

    options.allPlates.forEach((plateFile, idx) => {
      // For now, we'll just default to selecting the active plate,
      // since the old mapping indicators have been removed.
      const isSelected = (idx === this.stateManager.currentState?.activePlateIdx);

      const cb = ui.input.bool(plateFile.plate.barcode ?? `Plate ${idx + 1}`, {
        value: isSelected
      });
      plateCheckboxes.push(cb);
      plateSelectionHost.appendChild(cb.root);
    });

    const allPlatesCheckbox = ui.input.bool('Select All', {
      value: false,
      onValueChanged: (v) => {
        plateCheckboxes.forEach((cb) => cb.value = v);
      }
    });
    allPlatesCheckbox.value = plateCheckboxes.every((cb) => cb.value);

    platesSection.appendChild(allPlatesCheckbox.root);
    platesSection.appendChild(plateSelectionHost);

    content.appendChild(mappingsSection);
    if (plateCheckboxes.length > 0)
      content.appendChild(platesSection);

    dialog.add(content);
    dialog.onOK(() => {
      const selectedIndexes = plateCheckboxes
        .map((cb, idx) => cb.value ? idx : -1)
        .filter((idx) => idx !== -1);

      if (options.sourceMappings.size > 0)
        options.onSync(options.sourceMappings, selectedIndexes);
    });

    dialog.show();
  }

  private subscribeToStateChanges(): void {
    const sub = this.stateManager.onStateChange.subscribe((event) => {
      console.log('PlateGridManager: State changed', event);
      this.renderGrid();
    });
    this.subscriptions.push(sub);
  }

  destroy(): void {
    this.subscriptions.forEach((sub) => sub.unsubscribe());
  }
}
