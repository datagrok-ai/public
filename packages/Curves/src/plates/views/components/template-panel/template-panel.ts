import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {PlateStateManager} from '../../shared/plate-state-manager';
import {Subscription} from 'rxjs';
import {renderValidationResults} from '../../plates-validation-panel';
import {PlateTemplate, plateTemplates, plateTypes} from '../../../plates-crud';
import {PlateWidget} from '../../../../plate/plate-widget';
/**
 * Creates a consistent two-column form row (Label on left, Input on right).
 * @param label The text for the label.
 * @param input The Datagrok InputBase component.
 * @returns An HTMLElement representing the styled row.
 */
function createFormRow(label: string, input: DG.InputBase<any>): HTMLElement {
  const labelEl = ui.divText(label, 'ui-label');
  // Remove the default label from the input component, as we're providing our own.
  input.root.querySelector('label')?.remove();
  return ui.divH([labelEl, input.root], 'template-panel-form-row');
}

export class TemplatePanel {
  root: HTMLElement;
  private platePropertiesHost: HTMLElement;
  private validationHost: HTMLElement;
  private wellPropsHeaderHost: HTMLElement;
  private subscriptions: Subscription[] = [];

  // Hosts for the new Import section
  private plateIdentifierHost: HTMLElement;
  private replicateIdentifierHost: HTMLElement;
  private plateNumberHost: HTMLElement;

  constructor(
    private stateManager: PlateStateManager,
    private plateWidget: PlateWidget,
    private onManageMappings: () => void
  ) {
    this.root = ui.divV([], 'template-panel');
    this.platePropertiesHost = ui.divV([]);
    this.validationHost = ui.divV([]);
    this.wellPropsHeaderHost = ui.div();

    // Initialize new hosts
    this.plateIdentifierHost = ui.div([]);
    this.replicateIdentifierHost = ui.div([]);
    this.plateNumberHost = ui.div([]);

    this.buildPanel();
    this.subscribeToStateChanges();
  }

  private buildPanel(): void {
    const importSection = this.createImportSection(); // New section
    const templateSection = this.createTemplateSection();
    const platePropsSection = this.createCollapsiblePanel(
      ui.h2('Plate Properties'),
      this.platePropertiesHost,
      true
    );
    const wellPropsSection = this.createCollapsiblePanel(
      this.wellPropsHeaderHost,
      this.validationHost,
      true
    );

    this.root.appendChild(importSection); // Add new section to the top
    this.root.appendChild(templateSection);
    this.root.appendChild(platePropsSection);
    this.root.appendChild(wellPropsSection);
  }


  private createImportSection(): HTMLElement {
    const fileInput = this.createFileInput();
    const indexControls = ui.divV([
      this.plateIdentifierHost,
      this.replicateIdentifierHost,
      this.plateNumberHost,
    ], 'left-panel-section-content');

    const content = ui.divV([
      fileInput.root,
      indexControls
    ]);

    // Initially populate the controls (they will be empty)
    this.updateIdentifierControls();

    return this.createCollapsiblePanel(ui.h2('Import'), content, true);
  }

  private createFileInput(): DG.InputBase<DG.FileInfo | null> {
    const fileInput = ui.input.file('', {
      onValueChanged: async (file: DG.FileInfo | null) => {
        if (!file) return;
        try {
          const df = DG.DataFrame.fromCsv(await file.readAsString());
          await this.stateManager.loadDataFrame(df);
          // After loading the data, update the dropdowns with the new columns
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
      fileInputButton.appendChild(ui.divText('Import Plate File'));
    }
    fileInput.root.querySelector('label')?.remove();
    ui.tooltip.bind(fileInput.root, 'Import a plate from a CSV file');
    return fileInput;
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

    const choiceInput = ui.input.choice('', { // Label is now empty
      value: this.stateManager.plateIdentifierColumn,
      items: [null, ...df.columns.names()],
      onValueChanged: (newColumn: string | null) => {
        this.stateManager.plateIdentifierColumn = newColumn;
        this.stateManager.reprocessPlates();
      },
    });
    ui.tooltip.bind(choiceInput.root, 'Select the column that identifies individual plates in the file.');
    // Use the new helper to create a styled row
    this.plateIdentifierHost.appendChild(createFormRow('Plate Index', choiceInput));
  }

  private updateReplicateIdentifierControl(): void {
    ui.empty(this.replicateIdentifierHost);
    const df = this.stateManager.sourceDataFrame;
    if (!df) return;

    const choiceInput = ui.input.choice('', { // Label is now empty
      value: this.stateManager.replicateIdentifierColumn,
      items: [null, ...df.columns.names()],
      onValueChanged: (newColumn: string | null) => {
        this.stateManager.replicateIdentifierColumn = newColumn;
        this.stateManager.reprocessPlates();
      },
    });
    ui.tooltip.bind(choiceInput.root, 'Optional: Select the column that identifies technical replicates.');
    // Use the new helper to create a styled row
    this.replicateIdentifierHost.appendChild(createFormRow('Replicate Index', choiceInput));
  }

  private updatePlateNumberControl(): void {
    ui.empty(this.plateNumberHost);
    const df = this.stateManager.sourceDataFrame;
    if (!df) return;

    const choiceInput = ui.input.choice('', { // Label is now empty
      value: this.stateManager.plateNumberColumn,
      items: [null, ...df.columns.names()],
      onValueChanged: (newColumn: string | null) => {
        this.stateManager.plateNumberColumn = newColumn;
        this.stateManager.reprocessPlates();
      },
    });
    ui.tooltip.bind(choiceInput.root, 'Optional: Select the column that identifies the plate number.');
    // Use the new helper to create a styled row
    this.plateNumberHost.appendChild(createFormRow('Plate Number', choiceInput));
  }

  // --- EXISTING METHODS (UNCHANGED OR SLIGHTLY MODIFIED) ---


  private createTemplateSection(): HTMLElement {
    const templateIcon = ui.iconFA('file-alt', null, 'Template-defined properties');
    templateIcon.classList.add('legend-icon', 'legend-icon-template');

    const templateHeader = ui.divH(
      [ui.h2('Template'), templateIcon],
      {style: {alignItems: 'center', gap: '8px', flexGrow: '1'}}
    );

    const plateTypeSelector = ui.input.choice('', { // Label is now empty
      value: this.stateManager.currentPlateType.name,
      items: plateTypes.map((pt) => pt.name),
      onValueChanged: (v) => {
        const plateType = plateTypes.find((pt) => pt.name === v)!;
        this.stateManager.setPlateType(plateType);
        this.stateManager.setTemplate(this.stateManager.currentTemplate);
      },
    });

    const plateTemplateSelector = ui.input.choice('', { // Label is now empty
      value: this.stateManager.currentTemplate.name,
      items: plateTemplates.map((pt) => pt.name),
      onValueChanged: (v) => {
        const template = plateTemplates.find((pt) => pt.name === v)!;
        this.stateManager.setTemplate(template);
      },
    });

    // Create the rows using the new helper function
    const plateTypeRow = createFormRow('Plate Type', plateTypeSelector);
    const templateRow = createFormRow('Template', plateTemplateSelector);

    const templateContent = ui.divV(
      [plateTypeRow, templateRow], // Add the new rows to the content
      'left-panel-section-content'
    );

    return this.createCollapsiblePanel(templateHeader, templateContent, true);
  }

  private createCollapsiblePanel(
    header: HTMLElement,
    content: HTMLElement,
    expanded: boolean = true
  ): HTMLElement {
    // ... (this.createCollapsiblePanel method remains completely unchanged) ...
    const icon = ui.iconFA(expanded ? 'chevron-down' : 'chevron-right', () => {
      const isExpanded = content.style.display !== 'none';
      content.style.display = isExpanded ? 'none' : 'block';
      icon.classList.toggle('fa-chevron-down', !isExpanded);
      icon.classList.toggle('fa-chevron-right', isExpanded);
    });
    icon.style.marginRight = '8px';
    icon.style.color = 'var(--grey-4)';
    icon.style.cursor = 'pointer';

    header.style.cursor = 'pointer';
    header.onclick = () => icon.click();

    const headerContainer = ui.divH([icon, header]);
    headerContainer.style.alignItems = 'center';

    content.style.display = expanded ? 'block' : 'none';
    content.style.paddingLeft = '24px';
    content.classList.add('left-panel-section-content');

    return ui.divV([headerContainer, content], 'left-panel-section');
  }

  private subscribeToStateChanges(): void {
    const sub = this.stateManager.onStateChange.subscribe((event) => {
      // The identifier controls are now updated by the file input directly,
      // so no changes are needed here for them.
      this.updatePanelContent();
    });
    this.subscriptions.push(sub);
  }

  private updatePanelContent(): void {
    // ... (this.updatePanelContent method remains completely unchanged) ...
    const template = this.stateManager.currentTemplate;
    const activePlate = this.stateManager.activePlate;
    const state = this.stateManager.currentState;

    ui.empty(this.platePropertiesHost);
    ui.empty(this.validationHost);
    ui.empty(this.wellPropsHeaderHost);

    if (activePlate) {
      const handleMapping = (sourceColumn: string, targetProperty: string) => {
        const currentState = this.stateManager.currentState;
        if (currentState) {
          this.stateManager.remapProperty(
            currentState.activePlateIdx,
            targetProperty,
            sourceColumn
          );
        }
      };

      const handleUndo = (targetProperty: string) => {
        const currentState = this.stateManager.currentState;
        if (currentState)
          this.stateManager.undoMapping(currentState.activePlateIdx, targetProperty);
      };

      renderValidationResults(
        this.validationHost,
        activePlate.plate,
        template,
        handleMapping,
        activePlate.reconciliationMap,
        handleUndo
      );

      const plateProperties = template.plateProperties
        .filter((p) => p && p.name && p.type)
        .map((p) => DG.Property.js(p.name!, p.type! as DG.TYPE));

      if (plateProperties.length > 0) {
        const form = ui.input.form(activePlate.plate.details || {}, plateProperties);
        this.platePropertiesHost.appendChild(form);
      }
    } else {
      const templateProperties = template.plateProperties
        .filter((p) => p && p.name && p.type)
        .map((p) => DG.Property.js(p.name!, p.type! as DG.TYPE));

      if (templateProperties.length > 0) {
        const form = ui.input.form({}, templateProperties);
        this.platePropertiesHost.appendChild(form);
      }

      this.validationHost.appendChild(
        ui.divText('Import a CSV file to see validation results.', 'info-message')
      );
    }

    this.updateWellPropsHeader(template, state);
  }

  private updateWellPropsHeader(template: PlateTemplate, state: any): void {
    // ... (this.updateWellPropsHeader method remains completely unchanged) ...
    const manageMappingsButton = ui.button(
      ui.iconFA('exchange-alt'),
      this.onManageMappings,
      'Synchronize field mappings across multiple plates'
    );
    manageMappingsButton.classList.add('manage-mappings-button');

    if (state && state.plates.length > 0) {
      const badges = ui.divH([]);
      let totalConflicts = 0;
      let totalMappings = 0;

      state.plates.forEach((pf: any) => {
        const validation = renderValidationResults(
          ui.div(),
          pf.plate,
          template,
          () => {},
          pf.reconciliationMap,
          () => {}
        );
        totalConflicts += validation.conflictCount;
        totalMappings += pf.reconciliationMap.size;
      });

      if (totalConflicts > 0) {
        const conflictBadge = ui.span([`${totalConflicts}`], 'ui-badge-red');
        ui.tooltip.bind(conflictBadge, `${totalConflicts} unresolved fields across all plates`);
        badges.appendChild(conflictBadge);
      }

      if (totalMappings > 0) {
        const mappingBadge = ui.span([`${totalMappings}`], 'ui-badge-blue');
        ui.tooltip.bind(mappingBadge, `${totalMappings} total mappings applied across all plates`);
        badges.appendChild(mappingBadge);
      }

      if (badges.children.length > 0)
        manageMappingsButton.prepend(badges);
    }

    const header = ui.h2('Well Properties');
    const headerContainer = ui.divH(
      [header, ui.div([manageMappingsButton], {style: {marginLeft: 'auto'}})],
      'space-between-center'
    );

    this.wellPropsHeaderHost.appendChild(headerContainer);
  }

  destroy(): void {
    this.subscriptions.forEach((sub) => sub.unsubscribe());
  }
}
