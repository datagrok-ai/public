/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import { MolTrackProp } from '../utils/constants';
import { createPropertySection } from '../utils/view-utils';
import { fetchSchema, registerAssays } from '../package';
import { RegistrationViewBase } from './registration-view-base';

export class AssayRegistrationView extends RegistrationViewBase {
  assayInputs: DG.InputBase[] = [];
  assayResultInputs: DG.InputBase[] = [];

  constructor() {
    super('Register a new assay');
    this.path = 'Assay';
    this.buildUI();
  }

  protected async handleRegisterClick() {
    const assayValues = this.collectNonEmptyInputValues(this.assayInputs);
    const assayResultValues = this.assayResultInputs.map((inp) => ({
      name: inp.property.name,
      required: inp.value,
    }));

    assayValues['assay_result_properties'] = assayResultValues;

    try {
      const response = await registerAssays(JSON.stringify([assayValues]));
      const { status } = JSON.parse(response);
      this.showRegistrationMessage(status, assayValues.name);
    } catch (err) {
      this.showRegistrationMessage('error', (err as Error).message);
    }
  }

  private showRegistrationMessage(status: string, assayName?: string) {
    const isSuccess = status === 'success';
    const message = isSuccess ?
      `Assay ${assayName ?? ''} successfully registered!` :
      `Registration failed${assayName ? ` for ${assayName}` : ''}.`;

    this.showMessage(isSuccess, message, '');
  }

  private async fetchProperties(
    entityType: 'ASSAY' | 'ASSAY_RESULT',
    options?: {
    append?: MolTrackProp[],
    overrideValueType?: string
  },
  ): Promise<MolTrackProp[]> {
    const schema = await fetchSchema();
    const parsed = JSON.parse(schema);
    let props = parsed.filter((p: any) => p.entity_type === entityType);

    if (options?.overrideValueType)
      props = props.map((p: any) => ({ ...p, value_type: options.overrideValueType }));

    if (options?.append)
      props = [...options.append, ...props];

    return props;
  }

  private async fetchAssayProperties(): Promise<MolTrackProp[]> {
    return this.fetchProperties('ASSAY', { append: [{ name: 'name', value_type: 'string', nullable: false }] });
  }

  private async fetchAssayResultProperties(): Promise<MolTrackProp[]> {
    return this.fetchProperties('ASSAY_RESULT', { overrideValueType: 'bool' });
  }

  private async buildUI() {
    this.registerButton = ui.bigButton('REGISTER', () => this.handleRegisterClick());
    const clearAllButton = ui.iconFA('eraser', () => this.clearInputs(this.assayInputs), 'Clear all');
    this.view.setRibbonPanels([[this.registerButton, clearAllButton]]);

    const [assayProps, assayResultProps] = await Promise.all([
      this.fetchAssayProperties(),
      this.fetchAssayResultProperties(),
    ]);

    const { section: assaySection, inputs, formBackingObject } = await createPropertySection(
      'Assay properties',
      async () => assayProps,
      this.convertToDGProperty.bind(this),
      {
        initiallyOpen: true,
        onValidationChange: (invalid) => {
          this.invalidForm = invalid;
          this.registerButton?.classList.toggle('dim', this.invalidForm);
        },
      },
    );
    this.assayInputs = inputs;
    this.formBackingObject = { ...formBackingObject };

    const settings = ui.div();
    const df = DG.DataFrame.fromColumns(
      assayResultProps.map((p) => DG.Column.fromStrings(p.name, [p.name])),
    );

    const columnsEditor = ui.input.columns('Assay result properties', {
      table: df,
      value: df.columns.toList(),
      tooltipText: 'Select which assay result properties you want to configure below',
    });

    const updateSettings = () => {
      const chosen = columnsEditor.value.map((c) => c.name);
      ui.empty(settings);
      this.assayResultInputs = assayResultProps
        .filter((p) => chosen.includes(p.name))
        .map((p) => DG.InputBase.forProperty(this.convertToDGProperty(p)));

      if (this.assayResultInputs.length > 0) {
        const renderedTable = ui.table(
          this.assayResultInputs,
          (input: DG.InputBase, i: number) => [input.property.name, input.input],
          ['Property', 'Required'],
        );
        settings.appendChild(renderedTable);
      }
    };

    updateSettings();
    columnsEditor.onChanged.subscribe(updateSettings);

    const assayResultAccordion = ui.accordion();
    assayResultAccordion.addPane(
      'Assay result properties',
      () => ui.divV([columnsEditor.root, settings]),
      true,
    );

    const container = ui.divV([
      this.messageContainer,
      ui.divV([assaySection, assayResultAccordion.root], 'moltrack-top-sections'),
    ], 'moltrack-register-container');

    this.view.root.append(container);
  }
}
