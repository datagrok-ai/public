import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {ClinicalCaseViewBase} from '../model/ClinicalCaseViewBase';
import {ACT_TRT_ARM, PLANNED_TRT_ARM} from '../constants/columns-constants';
import {studies} from '../utils/app-utils';
import {calculateLBControlColumns} from '../data-preparation/data-preparation';

interface TreatmentControlRow {
  treatment: string;
  control: string;
}

export class ConfigurationView extends ClinicalCaseViewBase {
  private tableRows: TreatmentControlRow[] = [{treatment: '', control: ''}];
  private armChoices: string[] = [];
  private inputs: DG.InputBase[][] = []; // [row][cell] where cell 0=treatment, 1=control

  constructor(name: string, studyId: string) {
    super(name, studyId);
    this.name = name;
  }

  createView(): void {
    // Get ARM column categories from DM domain
    const dmDomain = studies[this.studyId].domains.dm;
    if (dmDomain) {
      // Try to get categories from PLANNED_TRT_ARM first, then ACT_TRT_ARM
      const armCol = dmDomain.col(PLANNED_TRT_ARM) || dmDomain.col(ACT_TRT_ARM);
      if (armCol?.categories)
        this.armChoices = armCol.categories;
    }

    // Load existing configuration if available
    const study = studies[this.studyId];
    if (study.treatmentAndControlConfig?.length > 0)
      this.tableRows = study.treatmentAndControlConfig.map((tc) => ({...tc}));

    // Create accordion
    const acc = ui.accordion();
    acc.addPane('Treatment and control', () => this.createTreatmentControlPane(), true);

    this.root.className = 'grok-view ui-box';
    this.root.append(acc.root);
  }

  private createTreatmentControlPane(): HTMLElement {
    const container = ui.divV([], {style: {padding: '10px'}});

    // Create table container
    const tableContainer = ui.div([], {style: {marginBottom: '10px'}});
    this.refreshTable(tableContainer);

    // Create controls container (plus icon and Apply button)
    const controlsContainer = ui.divH([], {style: {marginTop: '10px', gap: '10px', alignItems: 'center'}});

    // Plus icon button to add row
    const addRowButton = ui.button(ui.icons.add(() => {}), () => {
      // Save current values before adding new row
      this.saveCurrentTableValues();
      // Add row to bottom
      this.addTableRow();
      this.refreshTable(tableContainer);
    });
    addRowButton.title = 'Add row';

    // Apply button
    const applyButton = ui.button('Apply', () => {
      this.applyConfiguration();
    });

    controlsContainer.append(addRowButton);
    controlsContainer.append(applyButton);

    container.append(tableContainer);
    container.append(controlsContainer);

    return container;
  }

  private saveCurrentTableValues(): void {
    // Save values from stored input references
    this.inputs.forEach((rowInputs, rowIndex) => {
      if (rowIndex >= this.tableRows.length)
        return;

      if (rowInputs.length >= 2) {
        // Get treatment input value (cell 0)
        const treatmentInput = rowInputs[0];
        if (treatmentInput)
          this.tableRows[rowIndex].treatment = (treatmentInput.value as string) || '';

        // Get control input value (cell 1)
        const controlInput = rowInputs[1];
        if (controlInput)
          this.tableRows[rowIndex].control = (controlInput.value as string) || '';
      }
    });
  }

  private refreshTable(container: HTMLElement): void {
    // Save current values before refreshing
    this.saveCurrentTableValues();
    ui.empty(container);

    // Clear and rebuild inputs array
    this.inputs = [];

    const table = ui.table(this.tableRows, (row: TreatmentControlRow, idx: number) => {
      const index = idx - 1;
      // Treatment input (cell 0)
      const treatmentInput = ui.input.choice('', {
        value: row.treatment,
        items: this.armChoices,
        nullable: false,
      });
      treatmentInput.onChanged.subscribe((value: string) => {
        this.tableRows[index].treatment = value || '';
      });

      // Control input (cell 1)
      const controlInput = ui.input.choice('', {
        value: row.control,
        items: this.armChoices,
        nullable: false,
      });
      controlInput.onChanged.subscribe((value: string) => {
        this.tableRows[index].control = value || '';
      });

      // Store input references: [row][cell]
      if (!this.inputs[index])
        this.inputs[index] = [];
      this.inputs[index][0] = treatmentInput;
      this.inputs[index][1] = controlInput;

      return [treatmentInput.root, controlInput.root];
    }, ['Treatment', 'Control']);

    container.appendChild(table);
  }

  private addTableRow(): void {
    this.tableRows.push({treatment: '', control: ''});
  }

  private applyConfiguration(): void {
    // Save current values from the table before applying
    this.saveCurrentTableValues();

    // Validate that all rows have both treatment and control values
    const invalidRows = this.tableRows.filter((row) => !row.treatment || !row.control);
    if (invalidRows.length > 0) {
      grok.shell.warning('Please fill in both Treatment and Control for all rows.');
      return;
    }

    // Save to study configuration
    const study = studies[this.studyId];
    study.treatmentAndControlConfig = this.tableRows.map((row) => ({
      treatment: row.treatment,
      control: row.control,
    }));

    // Calculate control columns for LB domain if it exists
    if (study.domains.lb)
      calculateLBControlColumns(study.domains.lb, study.treatmentAndControlConfig);

    grok.shell.info(`Configuration saved: ${this.tableRows.length} treatment/control pair(s)`);
  }
}

