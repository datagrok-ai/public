import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {PlateStateManager} from '../../shared/plate-state-manager';
import {MappingDialogOptions} from '../../shared/types';
import {Subscription} from 'rxjs';
import {renderValidationResults} from '../../plates-validation-panel';

export class PlateGridManager {
  root: HTMLElement;
  grid: DG.Grid;
  private importContainer: HTMLElement;
  private plateIdentifierHost: HTMLElement;
  private replicateIdentifierHost: HTMLElement;
  private plateNumberHost: HTMLElement;
  private fileInput: DG.InputBase<DG.FileInfo | null>;
  private subscriptions: Subscription[] = [];

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
    this.root.appendChild(this.importContainer);
    this.root.appendChild(this.grid.root);
  }

  private initGrid(): void {
    // REMOVED: this.grid.root.style.display = 'none';
    this.grid.props.allowEdit = false;
    this.grid.props.allowRowSelection = true;
    this.grid.props.showCurrentRowIndicator = true;
    this.grid.props.showRowHeader = false;

    this.grid.onCurrentCellChanged.subscribe((gc: DG.GridCell) => {
      if (!gc.isTableCell)
        return;

      const state = this.stateManager.currentState;
      const tableRowIndex = this.grid.gridRowToTable(gc.gridRow);

      if (tableRowIndex !== -1 && state && tableRowIndex !== state.activePlateIdx)
        this.stateManager.selectPlate(tableRowIndex);
    });

    this.grid.onCellPrepare((gc) => {
      if (!gc.isTableCell)
        return;

      if (gc.tableColumn?.name === 'Maps' && gc.cell.value > 0) {
        const icon = ui.iconFA('exchange-alt');
        ui.tooltip.bind(icon, `${gc.cell.value} field(s) mapped`);
        const pill = ui.div(icon, 'plate-file-tab__mapped-pill');
        gc.style.element = ui.divH([pill], {style: {justifyContent: 'center', margin: 'auto 0'}});
        gc.customText = '';
      }

      if (gc.tableColumn?.name === 'Issues' && gc.cell.value > 0) {
        const icon = ui.div('', 'plate-file-tab__conflict-dot');
        ui.tooltip.bind(icon, `${gc.cell.value} unresolved fields`);
        gc.style.element = ui.divH([icon], {style: {justifyContent: 'center', margin: 'auto 0'}});
        gc.customText = '';
      }
    });
  }

  private renderGrid(): void {
    const state = this.stateManager.currentState;
    const template = this.stateManager.currentTemplate;

    if (!state || state.plates.length === 0) {
      this.grid.dataFrame = DG.DataFrame.create(0);
      return;
    }

    const barcodes = DG.Column.string('Barcode', state.plates.length)
      .init((i) => state.plates[i].plate.barcode ?? `Plate ${i + 1}`);
    const mappings = DG.Column.int('Mappings', state.plates.length)
      .init((i) => state.plates[i].reconciliationMap.size);
    const conflicts = DG.Column.int('Conflicts', state.plates.length)
      .init((i) => {
        const plateFile = state.plates[i];
        const validation = renderValidationResults(
          ui.div(), plateFile.plate, template, () => {},
          plateFile.reconciliationMap, () => {}
        );
        return validation.conflictCount;
      });

    const df = DG.DataFrame.fromColumns([barcodes, mappings, conflicts]);
    this.grid.dataFrame = df;

    const barcodeCol = this.grid.columns.byName('Barcode');
    if (barcodeCol) barcodeCol.width = 250;
    const mappingsCol = this.grid.columns.byName('Mappings');
    if (mappingsCol) {
      mappingsCol.width = 70;
      mappingsCol.name = 'Maps';
    }
    const conflictsCol = this.grid.columns.byName('Conflicts');
    if (conflictsCol) {
      conflictsCol.width = 70;
      conflictsCol.name = 'Issues';
    }

    if (state.activePlateIdx !== -1 && this.grid.dataFrame.currentRowIdx !== state.activePlateIdx)
      this.grid.dataFrame.currentRowIdx = state.activePlateIdx;
  }

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
    return ui.divH([
      ui.divText('Import Plate File'),
      this.fileInput.root,
      this.plateIdentifierHost,
      this.replicateIdentifierHost,
      this.plateNumberHost,
    ], 'plate-import-container');
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
    const sub = this.stateManager.onStateChange.subscribe(() => this.renderGrid());
    this.subscriptions.push(sub);
  }

  destroy(): void {
    this.subscriptions.forEach((sub) => sub.unsubscribe());
  }
}
