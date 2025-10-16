/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {PlateStateManager} from '../../shared/plate-state-manager';
import {Subscription} from 'rxjs';
import {PlateTemplate, plateTemplates, plateTypes} from '../../../plates-crud';
// import {PlateWidget} from '../../../../plate/plate-widget';
import {renderMappingEditor, TargetProperty} from '../mapping-editor/mapping-editor';
import {MAPPING_SCOPES} from '../../shared/scopes';
import './template-panel-and-mapping.css';
import { PlateWidget } from '../../../../plate/plate-widget/plate-widget';

function createFormRow(label: string, input: DG.InputBase<any>): HTMLElement {
  const labelEl = ui.divText(label, 'ui-label');
  input.root.querySelector('label')?.remove();
  return ui.divH([labelEl, input.root], 'assay-plates--form-row');
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
    this.root = ui.divV([], 'assay-plates--template-panel');
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
    const content = ui.divV([
      fileInput.root,
      this.identifierControlsHost
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
    fileInput.root.classList.add('assay-plates--plate-import-input');
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
    ui.empty(this.identifierControlsHost);
    const df = this.stateManager.sourceDataFrame;
    if (!df) return;

    const allDfColumns = df.columns.names();
    const selectedIdentifiers = this.stateManager.identifierColumns.filter((c): c is string => c !== null);

    this.stateManager.identifierColumns.forEach((colName, index) => {
      const otherSelected = selectedIdentifiers.filter((c) => c !== colName);
      const availableColumns = allDfColumns.filter((c) => !otherSelected.includes(c));

      const choiceInput = ui.input.choice<string | null>('', {
        value: colName,
        items: [null, ...availableColumns],
        onValueChanged: (newColumn) => {
          this.stateManager.setIdentifierColumn(index, newColumn);
        },
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

    const addBtn = ui.button(ui.iconFA('plus'), () => {
      this.stateManager.addIdentifierColumn();
    }, 'Add another index column');
    addBtn.classList.add('assay-plates--icon-button', 'assay-plates--add-button');
    this.identifierControlsHost.appendChild(ui.div(addBtn, 'assay-plates--add-identifier-container'));
  }

  private createTemplateSection(): HTMLElement {
    const templateIcon = ui.iconFA('file-alt', null, 'Template-defined properties');
    templateIcon.classList.add('legend-icon', 'legend-icon-template');
    const templateHeader = ui.divH(
      [ui.h2('Template'), templateIcon],
      'assay-plates--collapsible-header'
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
      'assay-plates--left-panel-section-content'
    );
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
    content.classList.add('assay-plates--left-panel-section-content', 'assay-plates--collapsible-content');

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

  private createPropertyInput(prop: any, currentValue: any): DG.InputBase {
  // Check if this is a choices property
    if (prop.choices) {
      let choices = prop.choices;
      if (typeof choices === 'string') {
        try {
          choices = JSON.parse(choices);
        } catch (e) {
          console.warn('Failed to parse choices:', choices);
        }
      }

      return ui.input.choice('', {
        value: currentValue || (Array.isArray(choices) ? choices[0] : null),
        items: Array.isArray(choices) ? choices : [choices],
      });
    }

    const property = DG.Property.fromOptions({
      ...prop,
      choices: prop.choices ?
        (typeof prop.choices === 'string' ? JSON.parse(prop.choices) : prop.choices) :
        undefined
    });

    return ui.input.forProperty(property, null, {value: currentValue});
  }

  private updatePanelContent(): void {
    const template = this.stateManager.currentTemplate;
    const activePlate = this.stateManager.activePlate;
    const state = this.stateManager.currentState;
    ui.empty(this.platePropertiesHost);
    ui.empty(this.wellPropsHeaderHost);

    if (activePlate && state) {
      const handleMapping = (targetProperty: string, sourceColumn: string) => {
        this.stateManager.setMapping(state.activePlateIdx, MAPPING_SCOPES.TEMPLATE, targetProperty, sourceColumn);
      };

      const handleUndo = (targetProperty: string) => {
        this.stateManager.removeMapping(state.activePlateIdx, MAPPING_SCOPES.TEMPLATE, targetProperty);
      };

      const identifierColumns = this.stateManager.identifierColumns.filter((c): c is string => c !== null);
      const allSourceColumns = activePlate.plate.data.columns.names();
      const availableSourceColumns = allSourceColumns.filter((c) => !identifierColumns.includes(c));
      const currentMappings = this.stateManager.getMappings(state.activePlateIdx, MAPPING_SCOPES.TEMPLATE);

      const MOCKED_REQUIRED_TEMPLATE_FIELDS = ['Target', 'Assay Format'];
      const templateProps: TargetProperty[] = template.wellProperties
        .filter((p) => p && p.name)
        .map((p) => ({
          name: p.name!,
          required: MOCKED_REQUIRED_TEMPLATE_FIELDS.includes(p.name!)
        }));

      renderMappingEditor(this.validationHost, {
        targetProperties: templateProps,
        sourceColumns: availableSourceColumns,
        mappings: currentMappings,
        onMap: handleMapping,
        onUndo: handleUndo,
      });

      const platePropertyInputs: HTMLElement[] = [];
      for (const prop of template.plateProperties) {
        if (!prop || !prop.name || !prop.type) continue;
        const currentValue = activePlate.plate.details?.[prop.name!];
        const input = this.createPropertyInput(prop, currentValue);
        input.onChanged.subscribe(() => {
          if (!activePlate.plate.details)
            activePlate.plate.details = {};
          activePlate.plate.details[prop.name!] = input.value;
        });
        platePropertyInputs.push(createFormRow(prop.name!, input));
      }

      if (platePropertyInputs.length > 0) {
        const form = ui.divV(platePropertyInputs);
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

      const platePropertyInputs: HTMLElement[] = [];
      for (const prop of template.plateProperties) {
        if (!prop || !prop.name || !prop.type) continue;
        const input = this.createPropertyInput(prop, null);
        platePropertyInputs.push(createFormRow(prop.name!, input));
      }

      if (platePropertyInputs.length > 0) {
        const form = ui.divV(platePropertyInputs);
        this.platePropertiesHost.appendChild(form);
      }
    }
    this.updateWellPropsHeader(template, state);
  }

  destroy(): void {
    this.subscriptions.forEach((sub) => sub.unsubscribe());
  }
}
