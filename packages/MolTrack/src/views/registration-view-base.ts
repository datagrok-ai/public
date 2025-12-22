/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MolTrackProp} from '../utils/constants';
import {buildPropertyOptions} from '../utils/utils';

let openedView: DG.ViewBase | null = null;

export abstract class RegistrationViewBase {
  view: DG.View;
  messageContainer: HTMLDivElement = ui.div([], 'moltrack-info-container');
  registerButton?: HTMLButtonElement;
  inputs: DG.InputBase[] = [];
  formBackingObject: Record<string, any> = {};
  invalidForm: boolean = false;
  title: string = 'Registration';
  path: string = '';

  constructor(title?: string) {
    this.view = DG.View.create();
    if (title) this.title = title;

    const titleText = ui.divText(this.title, 'moltrack-title');
    this.messageContainer.appendChild(titleText);
  }

  protected showMessage(isSuccess: boolean, title: string, message?: string) {
    const infoDiv = ui.info(message ?? '', title, true);
    const bar = infoDiv.querySelector('.grok-info-bar') as HTMLElement;
    if (bar) {
      bar.classList.toggle('moltrack-bar-success', isSuccess);
      bar.classList.toggle('moltrack-bar-error', !isSuccess);
    }
    ui.empty(this.messageContainer);
    this.messageContainer.appendChild(infoDiv);
  }

  protected collectNonEmptyInputValues(inputs: DG.InputBase[]): Record<string, any> {
    return inputs
      .filter((input) => input.value !== null && input.value !== undefined && input.value !== '')
      .reduce((acc, input) => {
        acc[input.property.name] = input.stringValue;
        return acc;
      }, {} as Record<string, any>);
  }

  protected clearInputs(inputs: DG.InputBase[]) {
    inputs.forEach((inp) => inp.value = null);
  }

  protected convertToDGProperty(prop: MolTrackProp, options?: any): DG.Property {
    return DG.Property.fromOptions(buildPropertyOptions(prop, options));
  }

  protected async createPropertySection(
    title: string,
    fetchPropsFn: () => Promise<any>,
    convertToDGProperty: (prop: MolTrackProp, options?: any) => DG.Property,
    options?: {
    disableNames?: string[];
    initiallyOpen?: boolean;
    onValidationChange?: (invalid: boolean) => void;
    reservedProperties?: string[];
    generateExample?: (pattern: string) => string;
  },
  ): Promise<{ section: HTMLElement; inputs: DG.InputBase[]; formBackingObject: Record<string, any> }> {
    const {
      disableNames = [],
      initiallyOpen = false,
      onValidationChange,
      reservedProperties = [],
      generateExample = () => '',
    } = options ?? {};

    const disableAll = disableNames.includes('*');
    let props: DG.Property[] = [];
    let propArray: any[] = [];
    const formBackingObject: Record<string, any> = {};

    try {
      const rawProps = await fetchPropsFn();
      const parsed = typeof(rawProps) === 'string' ? JSON.parse(rawProps) : rawProps;
      propArray = parsed.properties ?? parsed;

      props = propArray.map(convertToDGProperty);
      for (const p of props)
        formBackingObject[p.name] = null;
    } catch (err) {
      grok.shell.error(`Failed to fetch properties for "${title}": ${err}`);
    }

    const inputs = props.map((p) => {
      const input = DG.InputBase.forProperty(p, formBackingObject);
      input.onChanged.subscribe(() => {
        const invalid = !input.validate();
        onValidationChange?.(invalid);
      });

      const rawProp = propArray.find((rp) => rp.name === p.name);
      if (rawProp?.pattern) {
        const message = reservedProperties.includes(rawProp.name) ?
          'will be assigned at registration' :
          `e.g., ${generateExample(rawProp.pattern)}`;
        (input.input as HTMLInputElement).placeholder = message;
      }

      if (disableAll || disableNames.includes(p.name))
        input.readOnly = true;

      return input;
    });

    const form = ui.wideForm(inputs);
    form.classList.add('moltrack-compound-form', 'moltrack-form');

    const acc = ui.accordion();
    const accPane = acc.addPane(title, () => form);
    accPane.expanded = initiallyOpen;
    const section = accPane.root;

    return {section, inputs, formBackingObject};
  }

  show() {
    openedView?.close();
    openedView = this.view;
    grok.shell.addPreview(this.view);
  }

  protected abstract handleRegisterClick(): Promise<void>;
}
