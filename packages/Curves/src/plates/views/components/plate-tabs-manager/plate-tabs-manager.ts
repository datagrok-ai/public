import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {PlateStateManager} from '../../shared/plate-state-manager';
import {MappingDialogOptions, PlateFile, TemplateState} from '../../shared/types';
import {Subscription} from 'rxjs';
import {PlateTemplate} from '../../../plates-crud';
import {renderValidationResults} from '../../plates-validation-panel';

export class PlateTabsManager {
  root: HTMLElement;
  private importContainer: HTMLElement;
  private tabHeaderContainer: HTMLElement;
  private plateIdentifierHost: HTMLElement;
  private replicateIdentifierHost: HTMLElement;
  private plateNumberHost: HTMLElement;
  private fileInput: DG.InputBase<DG.FileInfo | null>;
  private subscriptions: Subscription[] = [];

  constructor(private stateManager: PlateStateManager) {
    this.root = ui.divV([], 'plate-tabs-manager');
    this.tabHeaderContainer = ui.divH([], 'plate-file-tabs-container');
    this.plateIdentifierHost = ui.div([]);
    this.replicateIdentifierHost = ui.div([]);
    this.plateNumberHost = ui.div([]);

    this.fileInput = this.createFileInput();
    this.importContainer = this.createImportContainer();

    this.buildComponent();
    this.subscribeToStateChanges();
  }

  private buildComponent(): void {
    this.root.appendChild(this.importContainer);
    this.root.appendChild(this.tabHeaderContainer);
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

  private renderFileTabs(template: PlateTemplate, state: TemplateState | undefined): void {
    ui.empty(this.tabHeaderContainer);

    if (!state || state.plates.length === 0) {
      this.tabHeaderContainer.style.visibility = 'hidden';
      return;
    }

    this.tabHeaderContainer.style.visibility = 'visible';
    this.tabHeaderContainer.style.display = 'flex';

    state.plates.forEach((plateFile, idx) => {
      const dummyTable = ui.div();
      const validation = renderValidationResults(
        dummyTable,
        plateFile.plate,
        template,
        () => {},
        plateFile.reconciliationMap,
        () => {}
      );

      const hasConflicts = validation.conflictCount > 0;
      const hasMappings = plateFile.reconciliationMap.size > 0;
      const tabLabelText = plateFile.plate.barcode ?? `Plate ${idx + 1}`;

      const headerItems = [];
      if (hasMappings) {
        const mappedIcon = ui.iconFA('exchange-alt');
        const pill = ui.div(mappedIcon, 'plate-file-tab__mapped-pill');
        ui.tooltip.bind(pill, `${plateFile.reconciliationMap.size} field(s) mapped`);
        headerItems.push(pill);
      }

      const labelElement = ui.divText(tabLabelText, 'plate-file-tab__label');
      headerItems.push(labelElement);

      const controlsContainer = ui.divH([], 'plate-file-tab__controls');
      headerItems.push(controlsContainer);

      if (hasConflicts) {
        const conflictDot = ui.div('', 'plate-file-tab__conflict-dot');
        ui.tooltip.bind(conflictDot, `${validation.conflictCount} unresolved fields`);
        controlsContainer.appendChild(conflictDot);
      }

      const clearIcon = ui.iconFA('times', (e: MouseEvent) => {
        e.stopPropagation();
        this.stateManager.removePlate(idx);
      }, 'Remove this plate');
      clearIcon.classList.add('plate-file-tab__clear-icon');
      controlsContainer.appendChild(clearIcon);

      const header = ui.divH(headerItems, 'plate-file-tab__header');
      header.classList.toggle('active', idx === state.activePlateIdx);

      header.onclick = (e) => {
        e.preventDefault();
        e.stopPropagation();
        this.stateManager.selectPlate(idx);
      };

      ui.tooltip.bind(header, () => plateFile.file.name);
      this.tabHeaderContainer.appendChild(header);
    });
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
      const state = this.stateManager.currentState;
      const template = this.stateManager.currentTemplate;
      this.renderFileTabs(template, state);
    });
    this.subscriptions.push(sub);
  }

  destroy(): void {
    this.subscriptions.forEach((sub) => sub.unsubscribe());
  }
}
