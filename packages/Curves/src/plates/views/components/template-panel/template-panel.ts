import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {PlateStateManager} from '../../shared/plate-state-manager';
import {Subscription} from 'rxjs';
import {renderValidationResults} from '../../plates-validation-panel';
import {PlateTemplate, plateTemplates, plateTypes} from '../../../plates-crud';
import {PlateWidget} from '../../../../plate/plate-widget';
import {toExcelPosition} from '../../../../plate/utils';

export class TemplatePanel {
  root: HTMLElement;
  private platePropertiesHost: HTMLElement;
  private validationHost: HTMLElement;
  private wellPropsHeaderHost: HTMLElement;
  private roleAssignmentSection: HTMLElement;
  private subscriptions: Subscription[] = [];

  constructor(
    private stateManager: PlateStateManager,
    private plateWidget: PlateWidget,
    private onManageMappings: () => void
  ) {
    this.root = ui.divV([], 'template-panel');
    this.platePropertiesHost = ui.divV([]);
    this.validationHost = ui.divV([]);
    this.wellPropsHeaderHost = ui.div();
    this.roleAssignmentSection = this.createRoleAssignmentSection();

    this.buildPanel();
    this.subscribeToStateChanges();
  }

  private buildPanel(): void {
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
    const roleSection = this.createCollapsiblePanel(
      ui.h2('Role Assignment'),
      this.roleAssignmentSection,
      true
    );

    this.root.appendChild(templateSection);
    this.root.appendChild(platePropsSection);
    this.root.appendChild(wellPropsSection);
    this.root.appendChild(roleSection);
  }

  private createTemplateSection(): HTMLElement {
    const templateIcon = ui.iconFA('file-alt', null, 'Template-defined properties');
    templateIcon.classList.add('legend-icon', 'legend-icon-template');

    const templateHeader = ui.divH(
      [ui.h2('Template'), templateIcon],
      {style: {alignItems: 'center', gap: '8px', flexGrow: '1'}}
    );

    const plateTypeSelector = ui.input.choice('Plate Type', {
      value: this.stateManager.currentPlateType.name,
      items: plateTypes.map((pt) => pt.name),
      onValueChanged: (v) => {
        const plateType = plateTypes.find((pt) => pt.name === v)!;
        this.stateManager.setPlateType(plateType);
        this.stateManager.setTemplate(this.stateManager.currentTemplate);
      }
    });

    plateTypeSelector.root.style.display = 'flex';
    plateTypeSelector.root.style.justifyContent = 'space-between';
    plateTypeSelector.root.style.alignItems = 'center';

    const plateTemplateSelector = ui.input.choice('Template', {
      value: this.stateManager.currentTemplate.name,
      items: plateTemplates.map((pt) => pt.name),
      onValueChanged: (v) => {
        const template = plateTemplates.find((pt) => pt.name === v)!;
        this.stateManager.setTemplate(template);
      }
    });

    plateTemplateSelector.root.style.display = 'flex';
    plateTemplateSelector.root.style.justifyContent = 'space-between';
    plateTemplateSelector.root.style.alignItems = 'center';

    const templateContent = ui.divV(
      [plateTypeSelector.root, plateTemplateSelector.root],
      'left-panel-section-content'
    );

    return this.createCollapsiblePanel(templateHeader, templateContent, true);
  }

  private createRoleAssignmentSection(): HTMLElement {
    const selectionStatusDiv = ui.divText('Drag on the plate to select wells.');
    selectionStatusDiv.style.marginTop = '10px';

    const roles = ['Control', 'Buffer', 'Assay Reagent', 'Sample'];
    const roleInput = ui.input.multiChoice('Roles', {items: roles, value: []});

    const assignButton = ui.button('Assign Roles', async () => {
      const selectedRoles = roleInput.value;
      const plate = this.plateWidget.plate;
      const selection = plate.data.selection;

      if (selectedRoles?.length === 0 || selection.trueCount === 0) return;

      const roleCol = plate.data.columns.getOrCreate('Role', DG.COLUMN_TYPE.STRING);
      const selectedIndexes = selection.getSelectedIndexes();

      for (const i of selectedIndexes)
        roleCol.set(i, selectedRoles?.join(',') ?? '');


      this.plateWidget._colorColumn = roleCol;
      this.plateWidget.grid.invalidate();
      this.plateWidget.updateRoleSummary();

      selection.setAll(false, true);
    });
    assignButton.disabled = true;

    const updateAssignButtonState = () => {
      const selectionCount = this.plateWidget.plate.data.selection.trueCount;
      const rolesCount = roleInput.value?.length ?? 0;
      assignButton.disabled = selectionCount === 0 || rolesCount === 0;
    };

    this.plateWidget.plate.data.selection.onChanged.subscribe(() => {
      const selectionCount = this.plateWidget.plate.data.selection.trueCount;

      if (selectionCount > 0) {
        selectionStatusDiv.textContent = `${selectionCount} wells selected`;
        const selectedPositions = Array.from(
          this.plateWidget.plate.data.selection.getSelectedIndexes()
        ).map((idx) => {
          const [row, col] = this.plateWidget.plate.rowIndexToExcel(idx);
          return toExcelPosition(row, col);
        });
        ui.tooltip.bind(selectionStatusDiv, () => ui.divText(selectedPositions.join(', ')));
      } else {
        selectionStatusDiv.textContent = 'Drag on the plate to select wells.';
        ui.tooltip.bind(selectionStatusDiv, () => null);
      }
      updateAssignButtonState();
    });

    roleInput.onInput.subscribe(() => updateAssignButtonState());

    return ui.divV(
      [selectionStatusDiv, roleInput.root, assignButton],
      'left-panel-section-content'
    );
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

    // Clear content
    ui.empty(this.platePropertiesHost);
    ui.empty(this.validationHost);
    ui.empty(this.wellPropsHeaderHost);

    const styleFormInputs = (formElement: HTMLElement) => {
      for (const child of Array.from(formElement.children)) {
        if (child.classList.contains('d4-input-base')) {
          (child as HTMLElement).style.display = 'flex';
          (child as HTMLElement).style.justifyContent = 'space-between';
          (child as HTMLElement).style.alignItems = 'center';
        }
      }
      return formElement;
    };

    if (activePlate) {
      const handleMapping = (currentColName: string, templatePropName: string) => {
        const currentState = this.stateManager.currentState;
        if (currentState) {
          this.stateManager.applyMapping(
            currentState.activePlateIdx,
            templatePropName,
            currentColName
          );
        }
      };

      const handleUndo = (mappedField: string) => {
        const currentState = this.stateManager.currentState;
        if (currentState)
          this.stateManager.undoMapping(currentState.activePlateIdx, mappedField);
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
        this.platePropertiesHost.appendChild(styleFormInputs(form));
      }
    } else {
      const templateProperties = template.plateProperties
        .filter((p) => p && p.name && p.type)
        .map((p) => DG.Property.js(p.name!, p.type! as DG.TYPE));

      if (templateProperties.length > 0) {
        const form = ui.input.form({}, templateProperties);
        this.platePropertiesHost.appendChild(styleFormInputs(form));
      }


      this.validationHost.appendChild(
        ui.divText('Import a CSV file to see validation results.', 'info-message')
      );
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
  }
}

