/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MolTrackProp} from '../utils/constants';
import {fetchSchema, registerAssays} from '../package';
import {RegistrationViewBase} from './registration-view-base';

export class AssayRegistrationView extends RegistrationViewBase {
  assayInputs: DG.InputBase[] = [];
  selectedAssayResultColumns: string[] = [];
  selectedAssayResultRequiredColumns: string[] = [];

  constructor() {
    const entityType = 'Assay';
    const entityTypeLower = entityType.toLowerCase();

    super(`Register a new ${entityTypeLower}`);
    this.path = entityType;
    this.view.name = `Register an ${entityTypeLower}`;
    this.buildUI();
  }
  protected async handleRegisterClick() {
    const assayValues = this.collectNonEmptyInputValues(this.assayInputs);

    const assayResultValues: { [key: string]: any }[] = this.selectedAssayResultColumns.map((colName) => {
      const isRequired = this.selectedAssayResultRequiredColumns.includes(colName);
      return {
        name: colName,
        required: isRequired,
      };
    });

    assayValues['assay_result_properties'] = assayResultValues;

    try {
      const response = await registerAssays(JSON.stringify([assayValues]));
      const {status} = JSON.parse(response);
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
      props = props.map((p: any) => ({...p, value_type: options.overrideValueType}));

    if (options?.append)
      props = [...options.append, ...props];

    return props;
  }

  private async fetchAssayProperties(): Promise<MolTrackProp[]> {
    return this.fetchProperties('ASSAY', {append: [{name: 'name', value_type: 'string', nullable: false}]});
  }

  private async fetchAssayResultProperties(): Promise<MolTrackProp[]> {
    return this.fetchProperties('ASSAY_RESULT', {overrideValueType: 'bool'});
  }

  private names(columns: DG.Column[]): string[] {
    return columns.map((col) => col.name);
  }

  private async buildUI() {
    this.registerButton = ui.bigButton('REGISTER', () => this.handleRegisterClick());
    const clearAllButton = ui.iconFA('eraser', () => this.clearInputs(this.assayInputs), 'Clear all');
    this.view.setRibbonPanels([[this.registerButton, clearAllButton]]);

    const [assayProps, assayResultProps] = await Promise.all([
      this.fetchAssayProperties(),
      this.fetchAssayResultProperties(),
    ]);

    const {section: assaySection, inputs, formBackingObject} = await this.createPropertySection(
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
    this.formBackingObject = {...formBackingObject};

    const df = DG.DataFrame.fromColumns(
      assayResultProps.map((p) => DG.Column.fromStrings(p.name, [p.name])),
    );

    const columnsEditor = ui.input.columns('Columns', {
      value: df.columns.toList(),
      table: df,
      onValueChanged: (value) => {
        this.selectedAssayResultColumns = this.names(value);
      },
      additionalColumns: {
        'required': df.columns.toList(),
      },
      onAdditionalColumnsChanged: (values: { [key: string]: DG.Column[] }) => {
        this.selectedAssayResultRequiredColumns = this.names(values['required']);
      }
    });

    const assayResultAccordion = ui.accordion();
    assayResultAccordion.addPane(
      'Assay result properties',
      () => columnsEditor.root,
      true,
    );

    const container = ui.divV([
      this.messageContainer,
      ui.divV([assaySection, assayResultAccordion.root], 'moltrack-top-sections'),
    ], 'moltrack-register-container');

    this.view.root.append(container);
  }
}
