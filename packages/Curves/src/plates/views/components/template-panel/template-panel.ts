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
    this.plateIdentifierHost = ui.div([]);
    this.replicateIdentifierHost = ui.div([]);
    this.plateNumberHost = ui.div([]);
    this.buildPanel();
    this.subscribeToStateChanges();
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
    const indexControls = ui.divV([
      this.plateIdentifierHost,
      this.replicateIdentifierHost,
      this.plateNumberHost,
    ], 'left-panel-section-content');
    const content = ui.divV([
      fileInput.root,
      indexControls
    ]);
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
          this.updateIdentifierControls();
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

  private updateIdentifierControls(): void {
    this.updatePlateIdentifierControl();
    this.updateReplicateIdentifierControl();
    this.updatePlateNumberControl();
  }

  private updatePlateIdentifierControl(): void {
    ui.empty(this.plateIdentifierHost);
    const df = this.stateManager.sourceDataFrame;
    if (!df) return;
    const choiceInput = ui.input.choice('', {
      value: this.stateManager.plateIdentifierColumn,
      items: [null, ...df.columns.names()],
      onValueChanged: (newColumn: string | null) => {
        this.stateManager.plateIdentifierColumn = newColumn;
        this.stateManager.reprocessPlates();
      },
    });
    ui.tooltip.bind(choiceInput.root, 'Select the column that identifies individual plates in the file.');
    this.plateIdentifierHost.appendChild(createFormRow('Plate Index', choiceInput));
  }

  private updateReplicateIdentifierControl(): void {
    ui.empty(this.replicateIdentifierHost);
    const df = this.stateManager.sourceDataFrame;
    if (!df) return;
    const choiceInput = ui.input.choice('', {
      value: this.stateManager.replicateIdentifierColumn,
      items: [null, ...df.columns.names()],
      onValueChanged: (newColumn: string | null) => {
        this.stateManager.replicateIdentifierColumn = newColumn;
        this.stateManager.reprocessPlates();
      },
    });
    ui.tooltip.bind(choiceInput.root, 'Optional: Select the column that identifies technical replicates.');
    this.replicateIdentifierHost.appendChild(createFormRow('Replicate Index', choiceInput));
  }

  private updatePlateNumberControl(): void {
    ui.empty(this.plateNumberHost);
    const df = this.stateManager.sourceDataFrame;
    if (!df) return;
    const choiceInput = ui.input.choice('', {
      value: this.stateManager.plateNumberColumn,
      items: [null, ...df.columns.names()],
      onValueChanged: (newColumn: string | null) => {
        this.stateManager.plateNumberColumn = newColumn;
        this.stateManager.reprocessPlates();
      },
    });
    ui.tooltip.bind(choiceInput.root, 'Optional: Select the column that identifies the plate number.');
    this.plateNumberHost.appendChild(createFormRow('Plate Number', choiceInput));
  }

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

      const sourceColumns = activePlate.plate.data.columns.names();
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
        sourceColumns: sourceColumns,
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

  private updateWellPropsHeader(template: PlateTemplate, state: any): void {
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
        const plateMappings = new Map<string, string>();
        for (const layer of pf.plate.getLayerNames()) {
          const aliases = pf.plate.getAliases(layer);
          for (const alias of aliases)
            plateMappings.set(alias, layer);
        }

        const validation = renderValidationResults(
          ui.div(),
          pf.plate,
          template,
          () => {},
          plateMappings,
          () => {}
        );
        totalConflicts += validation.conflictCount;
        totalMappings += plateMappings.size;
      });

      if (totalConflicts > 0) {
        const conflictBadge = ui.span([`${totalConflicts}`], 'ui-badge-red');
        ui.tooltip.bind(conflictBadge, `${totalConflicts} unresolved template fields across all plates`);
        badges.appendChild(conflictBadge);
      }
      if (totalMappings > 0) {
        const mappingBadge = ui.span([`${totalMappings}`], 'ui-badge-blue');
        ui.tooltip.bind(mappingBadge, `${totalMappings} total template mappings applied across all plates`);
        badges.appendChild(mappingBadge);
      }
      if (badges.children.length > 0)
        manageMappingsButton.prepend(badges);
    }

    const header = ui.h2('Template Properties');
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
