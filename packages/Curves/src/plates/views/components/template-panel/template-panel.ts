/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {PlateStateManager} from '../../shared/plate-state-manager';
import {Subscription} from 'rxjs';
import {plateTemplates, plateTypes} from '../../../plates-crud';
import {MAPPING_SCOPES} from '../../shared/scopes';
import './template-panel.css';
import {PlateWidget} from '../../../../plate/plate-widget/plate-widget';
import {TargetProperty} from '../mapping-editor/mapping-editor';
function createFormRow(label: string | HTMLElement, input: DG.InputBase<any> | HTMLElement): HTMLElement {
  const labelEl = typeof label === 'string' ? ui.divText(label, 'ui-label') : label;
  const inputEl = (input instanceof DG.InputBase) ? input.root : input;
  inputEl.querySelector('label')?.remove();
  return ui.divH([labelEl, inputEl], 'assay-plates--form-row');
}

export class TemplatePanel {
  root: HTMLElement;
  private platePropertiesHost: HTMLElement;
  private wellPropertiesHost: HTMLElement;
  private subscriptions: Subscription[] = [];
  private identifierControlsHost: HTMLElement;

  constructor(
    private stateManager: PlateStateManager,
    private plateWidget: PlateWidget,
    private onManageMappings: () => void
  ) {
    this.root = ui.divV([], 'assay-plates--template-panel');
    this.platePropertiesHost = ui.divV([]);
    this.wellPropertiesHost = ui.divV([]);
    this.identifierControlsHost = ui.divV([], 'identifier-controls-host');
    this.buildPanel();
    this.subscribeToStateChanges();
  }

  private createPropertyInput(prop: any, currentValue: any, isRequired: boolean): DG.InputBase {
    if (prop.choices) {
      let choices = prop.choices;
      if (typeof choices === 'string') {
        try {
          choices = JSON.parse(choices);
        } catch (e) {
          console.warn('Failed to parse choices:', choices);
          choices = [];
        }
      }
      return ui.input.choice(prop.name || '', {
        value: currentValue || (Array.isArray(choices) ? choices[0] : null),
        items: Array.isArray(choices) ? choices : [choices],
      });
    }

    const property = DG.Property.fromOptions(prop);
    const input = DG.InputBase.forProperty(property);
    if (currentValue !== undefined)
      input.value = currentValue;

    input.nullable = !isRequired;
    return input;
  }

  private buildPanel(): void {
    const importSection = this.createImportSection();
    const templateSection = this.createTemplateSection();
    const platePropsSection = this.createCollapsiblePanel(
      ui.h2('Plate Properties'), this.platePropertiesHost, true);
    const wellPropsSection = this.createCollapsiblePanel(
      ui.h2('Well Properties'), this.wellPropertiesHost, true);

    this.root.append(importSection, templateSection, platePropsSection, wellPropsSection);
  }

  private createImportSection(): HTMLElement {
    const fileInput = this.createFileInput();
    const content = ui.divV([fileInput.root, this.identifierControlsHost]);
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
    fileInput.root.classList.add('assay-plates--plate-import-input');
    const fileInputButton = fileInput.root.querySelector('button');
    if (fileInputButton) {
      ui.empty(fileInputButton);
      fileInputButton.append(ui.iconFA('upload'), ui.divText('Import Plate File'));
    }
    fileInput.root.querySelector('label')?.remove();
    ui.tooltip.bind(fileInput.root, 'Import a plate from a CSV file');
    return fileInput;
  }

  private updateIdentifierControls(): void {
    ui.empty(this.identifierControlsHost);
    const df = this.stateManager.sourceDataFrame;
    if (!df) return;

    const allDfColumns = df.columns.names();
    const selectedIdentifiers = this.stateManager.identifierColumns.filter((c): c is string => c !== null);

    this.stateManager.identifierColumns.forEach((colName, index) => {
      const otherSelected = selectedIdentifiers.filter((c) => c !== colName);
      const availableColumns = allDfColumns.filter((c) => !otherSelected.includes(c));

      const choiceInput = ui.input.choice<string | null>(`Plate Index ${index + 1}`, {
        value: colName,
        items: [null, ...availableColumns],
        onValueChanged: (newColumn) => this.stateManager.setIdentifierColumn(index, newColumn),
      });
      ui.tooltip.bind(choiceInput.root, 'Select a column to identify unique plates.');

      const removeBtn = ui.button(ui.iconFA('trash-alt'), () => {
        this.stateManager.removeIdentifierColumn(index);
      }, 'Remove this index column');
      removeBtn.classList.add('assay-plates--icon-button', 'assay-plates--remove-button');

      const formRow = createFormRow(`Plate Index ${index + 1}`, choiceInput);
      const controlRow = ui.divH([formRow, removeBtn], 'assay-plates--identifier-control-row');
      this.identifierControlsHost.appendChild(controlRow);
    });

    const addBtn = ui.button(ui.iconFA('plus'), () => this.stateManager.addIdentifierColumn(), 'Add another index column');
    addBtn.classList.add('assay-plates--icon-button', 'assay-plates--add-button');
    this.identifierControlsHost.appendChild(ui.div(addBtn, 'assay-plates--add-identifier-container'));
  }

  private createTemplateSection(): HTMLElement {
    const templateHeader = ui.h2('Template');
    const plateTypeSelector = ui.input.choice('Plate Type', {
      value: this.stateManager.currentPlateType.name,
      items: plateTypes.map((pt) => pt.name),
      onValueChanged: (v) => {
        const plateType = plateTypes.find((pt) => pt.name === v)!;
        this.stateManager.setPlateType(plateType);
        this.stateManager.setTemplate(this.stateManager.currentTemplate);
      },
    });
    const plateTemplateSelector = ui.input.choice('Template', {
      value: this.stateManager.currentTemplate.name,
      items: plateTemplates.map((pt) => pt.name),
      onValueChanged: (v) => {
        const template = plateTemplates.find((pt) => pt.name === v)!;
        this.stateManager.setTemplate(template);
      },
    });
    const templateContent = ui.divV([
      createFormRow('Plate Type', plateTypeSelector),
      createFormRow('Template', plateTemplateSelector),
    ]);
    return this.createCollapsiblePanel(templateHeader, templateContent, true);
  }

  private createCollapsiblePanel(header: HTMLElement, content: HTMLElement, expanded: boolean = true): HTMLElement {
    const icon = ui.iconFA(expanded ? 'chevron-down' : 'chevron-right', () => {
      const isExpanded = content.style.display !== 'none';
      content.style.display = isExpanded ? 'none' : 'block';
      icon.classList.toggle('fa-chevron-down', !isExpanded);
      icon.classList.toggle('fa-chevron-right', isExpanded);
    });
    icon.classList.add('assay-plates--collapsible-icon');

    header.onclick = () => icon.click();
    const headerContainer = ui.divH([icon, header], 'assay-plates--collapsible-header-container');

    content.style.display = expanded ? 'block' : 'none';
    content.classList.add('assay-plates--collapsible-content');

    return ui.divV([headerContainer, content], 'assay-plates--left-panel-section');
  }

  private subscribeToStateChanges(): void {
    const sub = this.stateManager.onStateChange$.subscribe((event) => {
      if (event.type === 'identifier-changed')
        this.updateIdentifierControls();
      else
        this.updatePanelContent();
    });
    this.subscriptions.push(sub);
  }

  private updatePanelContent(): void {
    const template = this.stateManager.currentTemplate;
    const activePlate = this.stateManager.activePlate;
    const state = this.stateManager.currentState;
    const requiredPropIds = new Set(template.required_props.map((tuple) => tuple[0]));

    ui.empty(this.platePropertiesHost);
    ui.empty(this.wellPropertiesHost);

    const platePropertyInputs: HTMLElement[] = [];
    for (const prop of template.plateProperties) {
      if (!prop || !prop.name || !prop.type) continue;

      const currentValue = activePlate?.plate.details?.[prop.name!];
      const isRequired = requiredPropIds.has(prop.id!);
      const input = this.createPropertyInput(prop, currentValue, isRequired);
      if (currentValue === undefined && input.value !== null && input.value !== undefined && activePlate) {
        if (!activePlate.plate.details) activePlate.plate.details = {};
        activePlate.plate.details[prop.name!] = input.value;
      }

      input.onChanged.subscribe(() => {
        if (activePlate) {
          if (!activePlate.plate.details) activePlate.plate.details = {};
          activePlate.plate.details[prop.name!] = input.value;
        }
      });

      const label = ui.label(prop.name!);
      if (isRequired) {
        label.classList.add('assay-plates--required-field');
        const asterisk = ui.span([' *'], 'assay-plates--required-asterisk');
        asterisk.style.color = 'red';
        asterisk.style.marginLeft = '2px';
        label.appendChild(asterisk);
      }
      platePropertyInputs.push(createFormRow(label, input));
    }
    if (platePropertyInputs.length > 0)
      this.platePropertiesHost.appendChild(ui.divV(platePropertyInputs));

    if (activePlate && state) {
      const handleMapping = (target: string, source: string) => this.stateManager.setMapping(state.activePlateIdx, MAPPING_SCOPES.TEMPLATE, target, source);
      const handleUndo = (target: string) => this.stateManager.removeMapping(state.activePlateIdx, MAPPING_SCOPES.TEMPLATE, target);
      const allSourceColumns = activePlate.plate.data.columns.names();
      const currentMappings = this.stateManager.getMappings(state.activePlateIdx, MAPPING_SCOPES.TEMPLATE);
      const templateProps: TargetProperty[] = template.wellProperties
        .filter((p) => p && p.name).map((p) => ({
          name: p.name!,
          required: requiredPropIds.has(p.id!)
        }));

      this.renderWellPropertyMappings(this.wellPropertiesHost, {
        targetProperties: templateProps, sourceColumns: allSourceColumns,
        mappings: currentMappings, onMap: handleMapping, onUndo: handleUndo,
      });
    } else {
      this.wellPropertiesHost.appendChild(ui.divText('Import a data file to map well properties.', 'assay-plates--info-message'));
    }
  }

  private createDynamicMappingRow(sourceColumns: string[], onMap: (target: string, source: string) => void, onCancel: () => void): HTMLElement {
    let propName: string | null = null; let sourceCol: string | null = null;
    const applyMapping = () => { if (propName && sourceCol) onMap(propName, sourceCol); };

    const propInput = ui.input.string('', {placeholder: 'Property name...'});
    propInput.onChanged.subscribe(() => { propName = propInput.value; });
    propInput.input.addEventListener('keydown', (e) => { if (e.key === 'Enter') { e.preventDefault(); applyMapping(); } });

    const colChoice = ui.input.choice('', {value: null, items: [null, ...sourceColumns], nullable: true, onValueChanged: (v: string|null) => { sourceCol = v; applyMapping(); }});
    const cancelBtn = ui.iconFA('times', onCancel, 'Cancel');
    cancelBtn.classList.add('assay-plates--mapping-cancel-icon');

    const inputContainer = ui.divH([colChoice.root, cancelBtn], 'assay-plates--mapping-input-container');
    const row = createFormRow(propInput.root, inputContainer);
    (propInput.input as HTMLElement).focus();
    return row;
  }

  private renderWellPropertyMappings(host: HTMLElement, options: any): void {
    const {targetProperties, sourceColumns, mappings, onMap, onUndo} = options;
    ui.empty(host);

    if (sourceColumns.length === 0) {
      host.appendChild(ui.divText('No columns available in the source file to map.', 'assay-plates--info-message'));
      return;
    }

    const allPropsMap = new Map<string, TargetProperty>();
    targetProperties.forEach((p: TargetProperty) => allPropsMap.set(p.name, p));
    mappings.forEach((_: string, targetCol: string) => {
      if (!allPropsMap.has(targetCol)) allPropsMap.set(targetCol, {name: targetCol, required: false});
    });

    const allTargetProps = Array.from(allPropsMap.values()).sort((a, b) => {
      if (a.required !== b.required) return a.required ? -1 : 1;
      return a.name.localeCompare(b.name);
    });

    allTargetProps.forEach((prop) => {
      const mappedSource = mappings.get(prop.name);
      const propNameEl = ui.span([prop.name]);
      if (prop.required) {
        const asterisk = ui.span([' *'], 'assay-plates--required-asterisk');
        asterisk.style.color = 'red';
        asterisk.style.marginLeft = '2px';
        propNameEl.appendChild(asterisk);
      }
      const choiceControl = ui.input.choice('', {
        value: mappedSource || null, items: [null, ...sourceColumns], nullable: true,
        onValueChanged: (v: string | null) => { if (v) onMap(prop.name, v); else if (mappedSource) onUndo(prop.name); },
      });
      const rightCell = ui.divH([choiceControl.root], 'assay-plates--mapping-input-container');
      if (mappedSource) {
        const undoIcon = ui.iconFA('times', () => onUndo(prop.name), 'Undo mapping');
        undoIcon.classList.add('assay-plates--mapping-undo-icon');
        rightCell.appendChild(undoIcon);
      }
      host.appendChild(createFormRow(propNameEl, rightCell));
    });

    const addRowHost = ui.divH([], 'assay-plates--add-row-container');
    const addBtn = ui.button(ui.iconFA('plus'), () => {
      const newRow = this.createDynamicMappingRow(sourceColumns, onMap, () => newRow.remove());
      host.insertBefore(newRow, addRowHost);
    }, 'Add new property mapping');
    addBtn.classList.add('assay-plates--icon-button', 'assay-plates--add-button');
    addRowHost.appendChild(addBtn);
    host.appendChild(addRowHost);
  }

  destroy(): void {
    this.subscriptions.forEach((sub) => sub.unsubscribe());
  }
}
