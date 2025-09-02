import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {PlateStateManager} from '../../shared/plate-state-manager';
import {Subscription} from 'rxjs';
// import {renderValidationResults} from '../../plates-validation-panel';
import {PlateTemplate, plateTemplates, plateTypes} from '../../../plates-crud';
import {PlateWidget} from '../../../../plate/plate-widget';
import {renderMappingEditor, TargetProperty} from '../mapping-editor/mapping-editor';
import {renderValidationResults} from '../../plates-validation-panel';
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

    const choiceInput = ui.input.choice('', { // Label is now empty
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
          this.stateManager.remapProperty(currentState.activePlateIdx, targetProperty, sourceColumn);
      };

      const handleUndo = (targetProperty: string) => {
        const currentState = this.stateManager.currentState;
        if (currentState)
          this.stateManager.undoMapping(currentState.activePlateIdx, targetProperty);
      };

      const sourceColumns = activePlate.plate.data.columns.names();

      // Define required fields for both existing and new analyses
      const DRC_FIELDS: TargetProperty[] = [
        {name: 'Activity', required: true},
        {name: 'Concentration', required: true},
        {name: 'SampleID', required: true},
      ];
      const DOSE_RATIO_FIELDS: TargetProperty[] = [
        {name: 'Percent_Inhibition', required: true},
        {name: 'Agonist_Concentration_M', required: true},
        {name: 'Antagonist_Concentration_M', required: true},
        // These are useful for titles/grouping but not strictly essential for calculation
        {name: 'Antagonist_ID', required: false},
        {name: 'Agonist_ID', required: false},
      ];
      const MOCKED_REQUIRED_TEMPLATE_FIELDS = ['Target', 'Assay Format'];
      const templateProps: TargetProperty[] = template.wellProperties
        .filter((p) => p && p.name)
        .map((p) => ({name: p.name!, required: MOCKED_REQUIRED_TEMPLATE_FIELDS.includes(p.name!)}));


      // Combine all properties into a single list for the mapping editor
      const allPropsMap = new Map<string, TargetProperty>();
      [...templateProps, ...DRC_FIELDS, ...DOSE_RATIO_FIELDS].forEach((prop) => {
        const existing = allPropsMap.get(prop.name);
        if (existing)
          existing.required = existing.required || prop.required;
        else
          allPropsMap.set(prop.name, {...prop});
      });
      const allTargetProps = Array.from(allPropsMap.values());

      renderMappingEditor(this.validationHost, {
        targetProperties: allTargetProps, // Pass the corrected list
        sourceColumns: sourceColumns,
        mappings: activePlate.reconciliationMap,
        onMap: handleMapping,
        onUndo: handleUndo,
      });


      // Plate properties logic remains the same...
      const plateProperties = template.plateProperties
        .filter((p) => p && p.name && p.type)
        .map((p) => DG.Property.js(p.name!, p.type! as DG.TYPE));
      if (plateProperties.length > 0)
        this.platePropertiesHost.appendChild(ui.input.form(activePlate.plate.details || {}, plateProperties));
    } else {
      // Handle the "no plate" case by calling the renderer with empty data
      renderMappingEditor(this.validationHost, {
        targetProperties: [],
        sourceColumns: [],
        mappings: new Map(),
        onMap: () => {},
        onUndo: () => {},
      });

      // Plate properties logic remains the same...
      const templateProperties = template.plateProperties
        .filter((p) => p && p.name && p.type)
        .map((p) => DG.Property.js(p.name!, p.type! as DG.TYPE));
      if (templateProperties.length > 0)
        this.platePropertiesHost.appendChild(ui.input.form({}, templateProperties));
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
    // No need to destroy a mapping editor instance anymore
  }
}
