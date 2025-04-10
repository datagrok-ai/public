/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Subject} from 'rxjs';
import type {ConsistencyInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import type {ValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';
import $ from 'cash-dom';

export type ValidationIconInput = {
  validation?: ValidationResult,
  consistency?: ConsistencyInfo,
};

export class ValidationIcon extends HTMLElement {
  private destroyed$ = new Subject<boolean>();

  private status?: ValidationIconInput;
  private scalar: boolean = true;
  private currentIcon?: HTMLElement;
  private currentPopover?: HTMLElement;

  constructor() {
    super();
    $(this).addClass('rfv2-validation-icon');
  }

  public destroy() {
    this.destroyed$.next(true);
    ui.empty(this);
  }

  connectedCallback() {
  }

  disconnectedCallback() {
  }

  set validationStatus(s: ValidationIconInput | undefined) {
    this.status = s;
    this.update();
  }

  get validationStatus() {
    return this.status;
  }

  set isScalar(v: boolean) {
    this.scalar = v;
  }

  get isScalar() {
    return !!this.scalar;
  }

  consistencyReset() {
    this.dispatchEvent(new CustomEvent('consistency-reset'));
  }

  requestAction(id: string) {
    this.dispatchEvent(new CustomEvent('action-request', {detail: id}));
  }

  update() {
    if (!this.status)
      return;
    const {validation, consistency} = this.status;

    while (this.firstChild && this.removeChild(this.firstChild));
    this.rerenderValidations();

    if (this.currentIcon && this.currentPopover) {
      this.appendChild(this.currentIcon);
      this.appendChild(this.currentPopover);
    }

    const isAnythingToShow =
      (validation?.errors && validation.errors.length) ||
      (validation?.warnings && validation.warnings.length) ||
      (validation?.notifications && validation.notifications.length) ||
      consistency?.inconsistent;

    $(this).css('display', isAnythingToShow ? 'flex': 'none');
  }

  rerenderValidations() {
    if (!this.status)
      return;
    const {validation, consistency} = this.status;

    const iconOptions = (() => {
      if (validation?.errors && validation.errors.length) return {name: 'exclamation-circle', color: 'var(--red-3)'};
      if (validation?.warnings && validation.warnings.length)
        return {name: 'exclamation-circle', color: 'var(--orange-2)'};
      if (validation?.notifications && validation.notifications.length)
        return {name: 'info-circle', color: 'var(--blue-1)'};
      if (consistency?.inconsistent)
        return {name: 'handshake-slash', color: 'var(--blue-1)'};

      return null;
    })();

    if (!iconOptions) return null;

    const icon = ui.iconFA(iconOptions.name, () => this.displayValidation());
    $(icon).css({'pointer-events': 'all'});
    if (iconOptions.color) $(icon).css('color', `${iconOptions.color}!important`);

    $(icon).toggleClass('fal far');
    const popover = this.addPopover(icon);

    this.currentIcon = icon;
    this.currentPopover = popover;
  }

  private renderValidationResults() {
    const root = ui.divV([], {style: {gap: '10px'}});

    const status = this?.status;
    if (status?.validation) {
      (['errors', 'warnings', 'notifications'] as const)
        .filter((category) => !!status.validation?.[category]?.length)
        .flatMap((category) => status.validation![category]!.map((advice) =>({category, advice})))
        .forEach(({category, advice}) => {
          const sectionIconOptions = this.getIconOptions(category)!;
          const icon = ui.iconFA(sectionIconOptions.name);
          $(icon).css({'color': sectionIconOptions.color, 'margin-right': '4px'});

          root.appendChild(ui.divV([
            ui.divH([
              icon,
              ui.divText(advice.description),
            ]),
            ...(advice.actions ?? []).map(
              (action) => ui.link(
                action.actionName,
                () => this.requestAction(action.action),
                undefined, {style: {paddingLeft: '20px'}})),
          ]));
        });
    }

    if (status?.consistency?.inconsistent) {
      const sectionIconOptions = this.getIconOptions('inconsistent')!;
      const icon = ui.iconFA(sectionIconOptions.name);
      $(icon).css({'color': sectionIconOptions.color, 'margin-right': '4px'});

      const consistentValue = status.consistency.assignedValue;
      root.appendChild(ui.divV([
        ui.divH([
          icon,
          ui.divText(`Current value is incosistent. Computed value was ${
            this.scalar ? consistentValue : 'different'
          }`),
        ]),
        ui.link(
          'Reset to consistent value',
          () => this.consistencyReset(),
          undefined, {style: {paddingLeft: '20px'}},
        ),
      ]));
    }

    return root;
  }

  private getIconOptions(category: 'errors' | 'warnings' | 'notifications' | 'inconsistent') {
    if (category === 'errors') return {name: 'exclamation-circle', color: 'var(--red-3)!important'};
    if (category === 'warnings')
      return {name: 'exclamation-circle', color: 'var(--orange-2)!important'};
    if (category === 'notifications')
      return {name: 'info-circle', color: 'var(--blue-1)!important'};
    if (category === 'inconsistent')
      return {name: 'handshake-slash', color: 'var(--blue-1)!important'};
  }

  private addPopover(icon: HTMLElement) {
    const popover = ui.div([], 'd4-tooltip');
    this.stylePopover(popover);
    icon.appendChild(popover);
    return popover;
  }

  private displayValidation() {
    const {currentIcon: icon, currentPopover: popover} = this;
    if (!icon || !popover)
      return;
    this.alignPopover(icon, popover);
    ui.empty(popover);
    const content = this.renderValidationResults();
    popover.appendChild(content);
    popover.showPopover();
  }

  private alignPopover(target: HTMLElement, popover: HTMLElement): void {
    const bounds = target.getBoundingClientRect().toJSON();
    popover.style.inset = 'unset';
    popover.style.top = bounds.y + 'px';
    popover.style.left = (bounds.x + 20) + 'px';
  }

  private stylePopover(popover: HTMLElement): void {
    popover.popover = 'auto';
    $(popover).css({
      'font-style': 'normal',
      'pointer-events': 'all',
      'max-width': '300px',
      'text-wrap': 'pretty',
    });
  }
}

export interface ValidationIconT extends ValidationIcon {};
