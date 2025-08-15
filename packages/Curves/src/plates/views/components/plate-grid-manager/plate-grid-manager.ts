import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {PlateStateManager} from '../../shared/plate-state-manager';
import {MappingDialogOptions, PlateFile} from '../../shared/types';
import {Subscription} from 'rxjs';

export class PlateGridManager {
  root: HTMLElement;
  grid: DG.Grid;
  private importContainer: HTMLElement;
  private plateIdentifierHost: HTMLElement;
  private replicateIdentifierHost: HTMLElement;
  private plateNumberHost: HTMLElement;
  private fileInput: DG.InputBase<DG.FileInfo | null>;
  private subscriptions: Subscription[] = [];
  private isSelecting: boolean = false;

  constructor(private stateManager: PlateStateManager) {
    this.root = ui.divV([], 'plate-grid-manager');
    this.plateIdentifierHost = ui.div([]);
    this.replicateIdentifierHost = ui.div([]);
    this.plateNumberHost = ui.div([]);

    this.fileInput = this.createFileInput();
    this.importContainer = this.createImportContainer();
    this.grid = DG.Grid.create(DG.DataFrame.create());

    this.buildComponent();
    this.initGrid();
    this.subscribeToStateChanges();
  }

  private buildComponent(): void {
    this.root.style.width = '100%';
    this.root.style.display = 'flex';
    this.root.style.flexDirection = 'column';
    this.root.appendChild(this.importContainer);
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
      if (this.isSelecting || !gc.isTableCell) return;
      const tableRowIndex = this.grid.gridRowToTable(gc.gridRow);
      const state = this.stateManager.currentState;
      if (tableRowIndex !== -1 && state && tableRowIndex !== state.activePlateIdx) {
        this.isSelecting = true;
        this.stateManager.selectPlate(tableRowIndex);
        this.isSelecting = false;
      }
    });

    this.grid.onCellClick.subscribe((gc: DG.GridCell) => {
      if (!gc.isTableCell) return;
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
        // --- FIX: 'g' was used before being assigned. Corrected order. ---
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
      this.grid.dataFrame = DG.DataFrame.create(0);
      return;
    }
    // --- FIX: Cast state.plates to the correct type ---
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

  // ... (rest of the file is unchanged) ...
  private createFileInput(): DG.InputBase<DG.FileInfo | null> {
    const fileInput = ui.input.file('', {
      onValueChanged: async (file: DG.FileInfo | null) => {
        if (!file) return;
        try {
          const df = DG.DataFrame.fromCsv(await file.readAsString());
          await this.stateManager.loadDataFrame(df);
          this.updateIdentifierControls();
        } catch (e: any) {
          grok.shell.error(`Failed to parse CSV: ${e.message}`);
        }
      },
    });

    fileInput.root.classList.add('plate-import-button');
    const fileInputButton = fileInput.root.querySelector('button');
    if (fileInputButton) {
      ui.empty(fileInputButton);
      fileInputButton.appendChild(ui.iconFA('upload'));
    }
    fileInput.root.querySelector('label')?.remove();
    ui.tooltip.bind(fileInput.root, 'Import a plate from a CSV file');
    return fileInput;
  }

  private createImportContainer(): HTMLElement {
    const container = ui.divH([
      ui.divText('Import Plate File'),
      this.plateIdentifierHost,
      this.replicateIdentifierHost,
      this.plateNumberHost,
      this.fileInput.root,
    ], 'plate-import-container');

    container.style.width = '100%';
    return container;
  }

  private updateIdentifierControls(): void {
    this.updatePlateIdentifierControl();
    this.updateReplicateIdentifierControl();
    this.updatePlateNumberControl();
  }

  private updatePlateIdentifierControl(): void {
    ui.empty(this.plateIdentifierHost);
    const df = this.stateManager.sourceDataFrame;
    if (!df) return;

    const choiceInput = ui.input.choice('Plate Index', {
      value: this.stateManager.plateIdentifierColumn,
      items: [null, ...df.columns.names()],
      onValueChanged: (newColumn: string | null) => {
        this.stateManager.plateIdentifierColumn = newColumn;
        this.stateManager.reprocessPlates();
      },
    });
    ui.tooltip.bind(choiceInput.root, 'Select the column that identifies individual plates in the file.');
    this.plateIdentifierHost.appendChild(choiceInput.root);
  }

  private updateReplicateIdentifierControl(): void {
    ui.empty(this.replicateIdentifierHost);
    const df = this.stateManager.sourceDataFrame;
    if (!df) return;

    const choiceInput = ui.input.choice('Replicate Index', {
      value: this.stateManager.replicateIdentifierColumn,
      items: [null, ...df.columns.names()],
      onValueChanged: (newColumn: string | null) => {
        this.stateManager.replicateIdentifierColumn = newColumn;
        this.stateManager.reprocessPlates();
      },
    });
    ui.tooltip.bind(choiceInput.root, 'Optional: Select the column that identifies technical replicates.');
    this.replicateIdentifierHost.appendChild(choiceInput.root);
  }

  private updatePlateNumberControl(): void {
    ui.empty(this.plateNumberHost);
    const df = this.stateManager.sourceDataFrame;
    if (!df) return;

    const choiceInput = ui.input.choice('Plate Number Index', {
      value: this.stateManager.plateNumberColumn,
      items: [null, ...df.columns.names()],
      onValueChanged: (newColumn: string | null) => {
        this.stateManager.plateNumberColumn = newColumn;
        this.stateManager.reprocessPlates();
      },
    });
    ui.tooltip.bind(choiceInput.root, 'Optional: Select the column that identifies the plate number within a barcode.');
    this.plateNumberHost.appendChild(choiceInput.root);
  }

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
      const isCurrentlyMapped = [...options.sourceMappings.keys()].every(
        (key) => plateFile.reconciliationMap.has(key)
      );
      const cb = ui.input.bool(plateFile.plate.barcode ?? `Plate ${idx + 1}`, {
        value: isCurrentlyMapped
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
