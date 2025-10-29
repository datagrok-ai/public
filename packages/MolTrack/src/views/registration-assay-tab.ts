/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import { MolTrackProp } from '../utils/constants';
import { createPropertySection } from '../utils/view-utils';
import { buildPropertyOptions } from '../utils/utils';
import { fetchSchema, registerAssays } from '../package';

let openedView: DG.ViewBase | null = null;

export class AssayRegistrationView {
  view: DG.View;
  assayInputs: DG.InputBase[] = [];
  assayResultInputs: DG.InputBase[] = [];
  formBackingObject: Record<string, any> = {};

  title = 'Register a new assay';
  invalidForm = true;
  path = 'Assay';

  private messageContainer = ui.div([], 'moltrack-info-container');
  private registerButton?: HTMLButtonElement;

  constructor() {
    this.view = DG.View.create();
    const titleText = ui.divText(this.title, 'moltrack-title');
    this.messageContainer.append(titleText);
    this.buildUI();
  }

  private convertToDGProperty(prop: MolTrackProp): DG.Property {
    return DG.Property.fromOptions(buildPropertyOptions(prop));
  }

  private resetInputs() {
    this.assayInputs.forEach((inp) => inp.value = null);
  }

  private async handleRegisterClick() {
    const assayValues = this.collectNonEmptyInputValues();
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

  private collectNonEmptyInputValues(): Record<string, any> {
    return Object.fromEntries(
      this.assayInputs
        .filter((inp) => inp.value != null && inp.value !== '')
        .map((inp) => [inp.property.name, inp.stringValue]),
    );
  }

  private showRegistrationMessage(status: string, assayName?: string) {
    const isSuccess = status === 'success';
    const message = isSuccess ?
      `Assay ${assayName ?? ''} successfully registered!` :
      `Registration failed${assayName ? ` for ${assayName}` : ''}.`;

    const infoDiv = ui.info('', message, true);
    const bar = infoDiv.querySelector('.grok-info-bar') as HTMLElement;

    if (bar) {
      bar.classList.toggle('moltrack-bar-success', isSuccess);
      bar.classList.toggle('moltrack-bar-error', !isSuccess);
    }

    ui.empty(this.messageContainer);
    this.messageContainer.append(infoDiv);
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
    const clearAllButton = ui.iconFA('eraser', () => this.resetInputs(), 'Clear all');
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

    this.assayResultInputs = assayResultProps.map((p) =>
      DG.InputBase.forProperty(this.convertToDGProperty(p)),
    );

    const assayResultAccordion = ui.accordion();
    assayResultAccordion.addPane('Assay result properties', () =>
      ui.wideForm(this.assayResultInputs), true);

    const container = ui.divV([
      this.messageContainer,
      ui.divV([assaySection, assayResultAccordion.root], 'moltrack-top-sections'),
    ], 'moltrack-register-container');

    this.view.root.append(container);
  }

  show() {
    openedView?.close();
    openedView = this.view;
    grok.shell.addPreview(this.view);
  }
}
