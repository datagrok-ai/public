import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {ACT_TRT_ARM, PLANNED_TRT_ARM} from '../constants/columns-constants';
import {studies} from '../utils/app-utils';
import {calculateLBControlColumns} from '../data-preparation/data-preparation';

interface TreatmentControlRow {
  treatment: string;
  control: string;
}

export class ConfigurationView extends DG.ViewBase {
  private tableRows: TreatmentControlRow[] = [{treatment: '', control: ''}];
  private armChoices: string[] = [];
  private inputs: DG.InputBase[][] = [];
  private studyId: string;
  loaded = false;

  constructor(name: string, studyId: string) {
    super();
    this.name = name;
    this.studyId = studyId;
  }

  load(): void {
    if (this.loaded)
      return;
    this.loaded = true;

    const dmDomain = studies[this.studyId].domains.dm;
    if (dmDomain) {
      const armCol = dmDomain.col(PLANNED_TRT_ARM) || dmDomain.col(ACT_TRT_ARM);
      if (armCol?.categories)
        this.armChoices = armCol.categories;
    }

    const study = studies[this.studyId];
    if (study.treatmentAndControlConfig?.length > 0)
      this.tableRows = study.treatmentAndControlConfig.map((tc) => ({...tc}));

    const acc = ui.accordion();
    acc.addPane('Treatment and control', () => this.createTreatmentControlPane(), true);

    this.root.className = 'grok-view ui-box';
    this.root.append(acc.root);
  }

  private createTreatmentControlPane(): HTMLElement {
    const container = ui.divV([], {style: {padding: '10px'}});

    const tableContainer = ui.div([], {style: {marginBottom: '10px'}});
    this.refreshTable(tableContainer);

    const controlsContainer = ui.divH([], {style: {marginTop: '10px', gap: '10px', alignItems: 'center'}});

    const addRowButton = ui.button(ui.icons.add(() => {}), () => {
      this.saveCurrentTableValues();
      this.addTableRow();
      this.refreshTable(tableContainer);
    });
    addRowButton.title = 'Add row';

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
    this.inputs.forEach((rowInputs, rowIndex) => {
      if (rowIndex >= this.tableRows.length)
        return;

      if (rowInputs.length >= 2) {
        const treatmentInput = rowInputs[0];
        if (treatmentInput)
          this.tableRows[rowIndex].treatment = (treatmentInput.value as string) || '';

        const controlInput = rowInputs[1];
        if (controlInput)
          this.tableRows[rowIndex].control = (controlInput.value as string) || '';
      }
    });
  }

  private refreshTable(container: HTMLElement): void {
    this.saveCurrentTableValues();
    ui.empty(container);

    this.inputs = [];

    const table = ui.table(this.tableRows, (row: TreatmentControlRow, idx: number) => {
      const index = idx - 1;
      const treatmentInput = ui.input.choice('', {
        value: row.treatment,
        items: this.armChoices,
        nullable: false,
      });
      treatmentInput.onChanged.subscribe(() => {
        this.tableRows[index].treatment = treatmentInput.value || '';
      });

      const controlInput = ui.input.choice('', {
        value: row.control,
        items: this.armChoices,
        nullable: false,
      });
      controlInput.onChanged.subscribe(() => {
        this.tableRows[index].control = treatmentInput.value || '';
      });

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
    this.saveCurrentTableValues();

    const invalidRows = this.tableRows.filter((row) => !row.treatment || !row.control);
    if (invalidRows.length > 0) {
      grok.shell.warning('Please fill in both Treatment and Control for all rows.');
      return;
    }

    const study = studies[this.studyId];
    study.treatmentAndControlConfig = this.tableRows.map((row) => ({
      treatment: row.treatment,
      control: row.control,
    }));

    if (study.domains.lb)
      calculateLBControlColumns(study.domains.lb, study.treatmentAndControlConfig);

    grok.shell.info(`Configuration saved: ${this.tableRows.length} treatment/control pair(s)`);
  }
}
