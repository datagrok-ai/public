/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {PlateStateManager} from '../../shared/plate-state-manager';
import {Subscription} from 'rxjs';
import {PlateTemplate, plateTemplates, plateTypes} from '../../../plates-crud';
import {PlateWidget} from '../../../../plate/plate-widget';
import {renderMappingEditor, TargetProperty} from '../mapping-editor/mapping-editor';
import {renderValidationResults} from '../../plates-validation-panel';
import {MAPPING_SCOPES} from '../../shared/scopes';

function createFormRow(label: string, input: DG.InputBase<any>): HTMLElement {
  const labelEl = ui.divText(label, 'ui-label');
  input.root.querySelector('label')?.remove();
  return ui.divH([labelEl, input.root], 'template-panel-form-row');
}

export class TemplatePanel {
  root: HTMLElement;
  private platePropertiesHost: HTMLElement;
  private validationHost: HTMLElement;
  private wellPropsHeaderHost: HTMLElement;
  private subscriptions: Subscription[] = [];
  private identifierControlsHost: HTMLElement;

  constructor(
  private stateManager: PlateStateManager,
  private plateWidget: PlateWidget,
  private onManageMappings: () => void
  ) {
    this.root = ui.divV([], 'template-panel');
    this.platePropertiesHost = ui.divV([]);
    this.validationHost = ui.divV([]);
    this.wellPropsHeaderHost = ui.div();
    this.identifierControlsHost = ui.divV([], 'identifier-controls-host');
    this.buildPanel();
    this.subscribeToStateChanges();
  }
  private updateWellPropsHeader(template: PlateTemplate, state: any): void {
    const header = ui.h2('Template Properties');
    const headerContainer = ui.divH(
      [header],
      'space-between-center'
    );
    this.wellPropsHeaderHost.appendChild(headerContainer);
  }

  private buildPanel(): void {
    const importSection = this.createImportSection();
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
    this.root.appendChild(importSection);
    this.root.appendChild(templateSection);
    this.root.appendChild(platePropsSection);
    this.root.appendChild(wellPropsSection);
  }

  private createImportSection(): HTMLElement {
    const fileInput = this.createFileInput();
    // START of CHANGE: Use the new single host for identifier controls
    const content = ui.divV([
      fileInput.root,
      this.identifierControlsHost
    ]);
    // END of CHANGE
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
          this.updateIdentifierControls(); // This will now get the auto-detected column
        } catch (e: any) {
          grok.shell.error(`Failed to parse CSV: ${e.message}`);
        }
      },
    });
    fileInput.root.classList.add('plate-import-input');
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

  // START of CHANGE: Replace the three update...Control methods with a single dynamic one.
  private updateIdentifierControls(): void {
    ui.empty(this.identifierControlsHost);
    const df = this.stateManager.sourceDataFrame;

    if (!df) return;

    // START of CHANGE: Logic to prevent duplicate index columns
    const allDfColumns = df.columns.names();
    const selectedIdentifiers = this.stateManager.identifierColumns.filter((c): c is string => c !== null);

    this.stateManager.identifierColumns.forEach((colName, index) => {
      // For each dropdown, the available choices are all columns MINUS those already selected in OTHER dropdowns.
      const otherSelected = selectedIdentifiers.filter((c) => c !== colName);
      const availableColumns = allDfColumns.filter((c) => !otherSelected.includes(c));

      const choiceInput = ui.input.choice<string | null>('', {
        value: colName,
        items: [null, ...availableColumns], // Use the filtered list
        onValueChanged: (newColumn) => {
          this.stateManager.setIdentifierColumn(index, newColumn);
        },
      });
      ui.tooltip.bind(choiceInput.root, 'Select a column to identify unique plates.');

      // START of CHANGE: Improved remove button UI
      const removeBtn = ui.button(ui.iconFA('trash-alt'), () => {
        this.stateManager.removeIdentifierColumn(index);
      }, 'Remove this index column');
      removeBtn.classList.add('curves-icon-button', 'curves-remove-button');
      // END of CHANGE

      const formRow = createFormRow(`Plate Index ${index + 1}`, choiceInput);
      const controlRow = ui.divH([formRow, removeBtn], {style: {alignItems: 'center', flexWrap: 'nowrap'}});
      (formRow.lastChild as HTMLElement).style.flexGrow = '1';

      this.identifierControlsHost.appendChild(controlRow);
    });

    // START of CHANGE: Improved add button UI
    const addBtn = ui.button(ui.iconFA('plus'), () => {
      this.stateManager.addIdentifierColumn();
    }, 'Add another index column');
    addBtn.classList.add('curves-icon-button', 'curves-add-button');
    // END of CHANGE

    this.identifierControlsHost.appendChild(ui.div(addBtn, {style: {marginTop: '8px'}}));
  }
  // END of CHANGE

  private createTemplateSection(): HTMLElement {
    const templateIcon = ui.iconFA('file-alt', null, 'Template-defined properties');
    templateIcon.classList.add('legend-icon', 'legend-icon-template');
    const templateHeader = ui.divH(
      [ui.h2('Template'), templateIcon],
      {style: {alignItems: 'center', gap: '8px', flexGrow: '1'}}
    );
    const plateTypeSelector = ui.input.choice('', {
      value: this.stateManager.currentPlateType.name,
      items: plateTypes.map((pt) => pt.name),
      onValueChanged: (v) => {
        const plateType = plateTypes.find((pt) => pt.name === v)!;
        this.stateManager.setPlateType(plateType);
        this.stateManager.setTemplate(this.stateManager.currentTemplate);
      },
    });
    const plateTemplateSelector = ui.input.choice('', {
      value: this.stateManager.currentTemplate.name,
      items: plateTemplates.map((pt) => pt.name),
      onValueChanged: (v) => {
        const template = plateTemplates.find((pt) => pt.name === v)!;
        this.stateManager.setTemplate(template);
      },
    });
    const plateTypeRow = createFormRow('Plate Type', plateTypeSelector);
    const templateRow = createFormRow('Template', plateTemplateSelector);
    const templateContent = ui.divV(
      [plateTypeRow, templateRow],
      'left-panel-section-content'
    );
    return this.createCollapsiblePanel(templateHeader, templateContent, true);
  }

  private createCollapsiblePanel(
    header: HTMLElement,
    content: HTMLElement,
    expanded: boolean = true
  ): HTMLElement {
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
      if (event.type === 'identifier-changed')
        this.updateIdentifierControls();
      else
        this.updatePanelContent();
    });
    this.subscriptions.push(sub);
  }

  private createPropertyInput(prop: any, currentValue: any): DG.InputBase {
    if (prop.choices) {
      let choicesList: string[];
      if (typeof prop.choices === 'string') {
        try {
          choicesList = JSON.parse(prop.choices);
        } catch {
          choicesList = [prop.choices];
        }
      } else if (Array.isArray(prop.choices)) {
        choicesList = prop.choices;
      } else {
        choicesList = [];
      }

      return ui.input.choice(prop.name, {
        items: choicesList,
        value: currentValue || choicesList[0]
      });
    }

    // For numeric properties with min/max
    if (prop.type === DG.TYPE.FLOAT || prop.type === DG.TYPE.INT) {
      const input = prop.type === DG.TYPE.INT ?
        ui.input.int(prop.name, {value: currentValue || 0}) :
        ui.input.float(prop.name, {value: currentValue || 0});

      if (prop.min !== undefined || prop.max !== undefined) {
        input.addValidator((valueStr: string) => {
          const value = parseFloat(valueStr);
          if (isNaN(value))
            return 'Invalid number';
          if (prop.min !== undefined && value < prop.min)
            return `Minimum value is ${prop.min}`;
          if (prop.max !== undefined && value > prop.max)
            return `Maximum value is ${prop.max}`;
          return null;
        });
      }

      return input;
    }

    // For boolean properties
    if (prop.type === DG.TYPE.BOOL)
      return ui.input.bool(prop.name, {value: currentValue || false});


    // Default to string input
    return ui.input.string(prop.name, {value: currentValue || ''});
  }

  private updatePanelContent(): void {
    const template = this.stateManager.currentTemplate;
    const activePlate = this.stateManager.activePlate;
    const state = this.stateManager.currentState;
    ui.empty(this.platePropertiesHost);
    ui.empty(this.wellPropsHeaderHost);


    if (activePlate) {
      const handleMapping = (targetProperty: string, sourceColumn: string) => {
        const currentState = this.stateManager.currentState;
        if (currentState)
          this.stateManager.remapScopedProperty(currentState.activePlateIdx, MAPPING_SCOPES.TEMPLATE, targetProperty, sourceColumn);
      };

      const handleUndo = (targetProperty: string) => {
        const currentState = this.stateManager.currentState;
        if (currentState)
          this.stateManager.undoScopedMapping(currentState.activePlateIdx, MAPPING_SCOPES.TEMPLATE, targetProperty);
      };
      const identifierColumns = this.stateManager.identifierColumns.filter((c): c is string => c !== null);
      const allSourceColumns = activePlate.plate.data.columns.names();
      const availableSourceColumns = allSourceColumns.filter((c) => !identifierColumns.includes(c));

      const currentMappings = activePlate.plate.getScopedAliases(MAPPING_SCOPES.TEMPLATE);

      const MOCKED_REQUIRED_TEMPLATE_FIELDS = ['Target', 'Assay Format'];
      const templateProps: TargetProperty[] = template.wellProperties
        .filter((p) => p && p.name)
        .map((p) => ({
          name: p.name!,
          required: MOCKED_REQUIRED_TEMPLATE_FIELDS.includes(p.name!)
        }));

      renderMappingEditor(this.validationHost, {
        targetProperties: templateProps,
        sourceColumns: availableSourceColumns, // Use the filtered list here
        mappings: currentMappings,
        onMap: handleMapping,
        onUndo: handleUndo,
      });

      // Create proper inputs for plate properties
      const platePropertyInputs: HTMLElement[] = [];
      for (const prop of template.plateProperties) {
        if (!prop || !prop.name || !prop.type) continue;

        const currentValue = activePlate.plate.details?.[prop.name!];
        const input = this.createPropertyInput(prop, currentValue);

        // Update the plate details when value changes
        input.onChanged.subscribe(() => {
          if (!activePlate.plate.details)
            activePlate.plate.details = {};
          activePlate.plate.details[prop.name!] = input.value;
        });

        platePropertyInputs.push(createFormRow(prop.name!, input));
      }

      if (platePropertyInputs.length > 0) {
        const form = ui.divV(platePropertyInputs, 'template-plate-properties-form');
        this.platePropertiesHost.appendChild(form);
      }
    } else {
      renderMappingEditor(this.validationHost, {
        targetProperties: [],
        sourceColumns: [],
        mappings: new Map(),
        onMap: () => {},
        onUndo: () => {},
      });

      // Create proper inputs even when no plate is active
      const platePropertyInputs: HTMLElement[] = [];
      for (const prop of template.plateProperties) {
        if (!prop || !prop.name || !prop.type) continue;

        const input = this.createPropertyInput(prop, null);
        platePropertyInputs.push(createFormRow(prop.name!, input));
      }

      if (platePropertyInputs.length > 0) {
        const form = ui.divV(platePropertyInputs, 'template-plate-properties-form');
        this.platePropertiesHost.appendChild(form);
      }
    }

    this.updateWellPropsHeader(template, state);
  }


  destroy(): void {
    this.subscriptions.forEach((sub) => sub.unsubscribe());
  }
}
