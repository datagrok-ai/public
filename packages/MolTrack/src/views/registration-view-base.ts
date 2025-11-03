/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import { MolTrackProp } from '../utils/constants';
import { buildPropertyOptions } from '../utils/utils';

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

  show() {
    openedView?.close();
    openedView = this.view;
    grok.shell.addPreview(this.view);
  }

  protected abstract handleRegisterClick(): Promise<void>;
}
