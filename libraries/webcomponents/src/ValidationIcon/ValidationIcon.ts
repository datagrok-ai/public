/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Subject, BehaviorSubject, combineLatest} from 'rxjs';
import {debounceTime, distinctUntilChanged, takeUntil} from 'rxjs/operators';
import {isLinkedRestriction} from '@datagrok-libraries/utils/src/validation-utils';

export type RestrictionType = 'disabled' | 'restricted' | 'info' | 'none';

export type ConsistencyInfo = {
  restriction: RestrictionType,
  inconsistent: boolean,
  assignedValue: any,
}

export interface ActionItem {
  actionName: string;
  action: string;
  additionalParams?: Record<string, any>;
}

export interface Advice {
  description: string;
  actions?: ActionItem[];
}

export type ValidationItem = string | Advice;

export interface ValidationResult {
  errors?: ValidationItem[];
  warnings?: ValidationItem[];
  notifications?: ValidationItem[];
}

// detail stays the action id, so hosts that read only detail keep working
export class ActionRequestEvent extends CustomEvent<string> {
  constructor(action: string, public additionalParams?: Record<string, any>) {
    super('action-request', {detail: action});
  }
}

export type ValidationIconInput = {
  validation?: ValidationResult,
  consistency?: ConsistencyInfo,
};

const POPOVER_CLASS = 'dg-validation-icon-popover';
const POPOVER_MAX_HEIGHT = 400;
const POPOVER_VIEWPORT_MARGIN = 8;

let popoverStylesInjected = false;
function injectPopoverStyles(): void {
  if (popoverStylesInjected)
    return;
  popoverStylesInjected = true;
  const style = document.createElement('style');
  style.textContent = `
    .${POPOVER_CLASS} { scrollbar-width: thin; scrollbar-color: var(--grey-3) transparent; }
    .${POPOVER_CLASS}::-webkit-scrollbar { width: 8px; }
    .${POPOVER_CLASS}::-webkit-scrollbar-track { background: transparent; }
    .${POPOVER_CLASS}::-webkit-scrollbar-thumb { background: var(--grey-3); border-radius: 4px; }
    .${POPOVER_CLASS}::-webkit-scrollbar-thumb:hover { background: var(--grey-4); }
  `;
  document.head.appendChild(style);
}

type IconSpec = {name: string, color: string};

const CATEGORIES = ['errors', 'warnings', 'notifications'] as const;

const CATEGORY_ICONS: Record<typeof CATEGORIES[number], IconSpec> = {
  errors: {name: 'exclamation-circle', color: 'var(--red-3)'},
  warnings: {name: 'exclamation-circle', color: 'var(--orange-2)'},
  notifications: {name: 'info-circle', color: 'var(--blue-1)'},
};

const inconsistentIcon = (restriction?: string): IconSpec =>
  ({name: 'unlink', color: isLinkedRestriction(restriction) ? 'var(--red-3)' : 'var(--blue-1)'});

// an inconsistent workflow-set value outranks validation, an inconsistent free value ranks last
function statusIcon({validation, consistency}: ValidationIconInput): IconSpec | undefined {
  if (consistency?.inconsistent && isLinkedRestriction(consistency.restriction))
    return inconsistentIcon(consistency.restriction);
  const category = CATEGORIES.find((c) => validation?.[c]?.length);
  if (category)
    return CATEGORY_ICONS[category];
  if (consistency?.inconsistent)
    return inconsistentIcon(consistency.restriction);
}

function coloredIcon(spec: IconSpec, handler?: () => void, tooltip?: string) {
  const icon = ui.iconFA(spec.name, handler, tooltip);
  icon.style.setProperty('color', spec.color, 'important');
  return icon;
}

function sectionIcon(spec: IconSpec) {
  const icon = coloredIcon(spec);
  icon.style.marginRight = '4px';
  return icon;
}

const indentedLink = (text: string, handler: () => void) =>
  ui.link(text, handler, undefined, {style: {paddingLeft: '20px'}});

export class ValidationIcon extends HTMLElement {
  isDataFrame = false;

  private destroyed$ = new Subject<boolean>();

  private status$ = new BehaviorSubject<ValidationIconInput | undefined>(undefined);
  private isScalar$ = new BehaviorSubject(true);
  private hover$ = new BehaviorSubject(false);
  private currentIcon?: HTMLElement;
  private currentPopover?: HTMLElement;

  constructor() {
    super();

    combineLatest([
      this.status$.pipe(distinctUntilChanged()),
      this.isScalar$.pipe(distinctUntilChanged()),
      this.hover$.pipe(distinctUntilChanged()),
    ]).pipe(
      debounceTime(0),
      takeUntil(this.destroyed$),
    ).subscribe(() => this.update());
  }

  public destroy() {
    this.destroyed$.next(true);
    ui.empty(this);
  }

  set validationStatus(s: ValidationIconInput | undefined) {
    this.status$.next(s);
  }

  get validationStatus() {
    return this.status$.value;
  }

  set useHover(v: boolean) {
    this.hover$.next(v);
  }

  get useHover() {
    return this.hover$.value;
  }

  set isScalar(v: boolean) {
    this.isScalar$.next(v);
  }

  get isScalar() {
    return this.isScalar$.value;
  }

  consistencyReset() {
    this.dispatchEvent(new CustomEvent('consistency-reset'));
  }

  requestAction(id: string, additionalParams?: Record<string, any>) {
    this.dispatchEvent(new ActionRequestEvent(id, additionalParams));
  }

  showDataFrameDiff(cval: any) {
    this.dispatchEvent(new CustomEvent('show-dataframe-diff', {detail: cval}));
  }

  private update() {
    const status = this.status$.value;
    const spec = status && statusIcon(status);
    this.currentIcon = undefined;
    this.currentPopover = undefined;
    if (spec && this.hover$.value) {
      this.currentIcon = coloredIcon(spec);
      ui.tooltip.bind(this.currentIcon, () => this.renderValidationResults());
    } else if (spec) {
      this.currentIcon = coloredIcon(spec, () => this.displayValidation(), 'Show validation');
      this.currentIcon.style.pointerEvents = 'all';
      this.currentIcon.classList.toggle('fal');
      this.currentIcon.classList.toggle('far');
      this.currentPopover = this.createPopover();
    } else if (status && this.isLinked(status.consistency)) {
      this.currentIcon = ui.iconFA('link', null, status.consistency!.restriction === 'disabled' ?
        'Set by the workflow' : 'Set by the workflow. Editing it makes the value inconsistent');
      this.currentIcon.style.setProperty('color', 'var(--grey-4)', 'important');
      this.currentIcon.style.pointerEvents = 'all';
    }
    // the popover is the icon's sibling, so a click inside it does not reach the icon's handler
    this.replaceChildren(...[this.currentIcon, this.currentPopover].filter((el): el is HTMLElement => !!el));
    this.style.display = this.currentIcon ? 'flex' : 'none';
  }

  private renderValidationResults() {
    const root = ui.divV([], {style: {gap: '10px'}});
    const {validation, consistency} = this.status$.value ?? {};

    for (const category of CATEGORIES) {
      for (const advice of validation?.[category] ?? []) {
        const description = typeof advice === 'string' ? advice : advice.description;
        const actions = typeof advice === 'string' ? [] : (advice.actions ?? []);
        root.appendChild(ui.divV([
          ui.divH([sectionIcon(CATEGORY_ICONS[category]), ui.divText(description)], {style: {alignItems: 'center'}}),
          ...actions.map((action) =>
            indentedLink(action.actionName, () => this.requestAction(action.action, action.additionalParams))),
        ]));
      }
    }

    if (consistency?.inconsistent) {
      const consistentValue = consistency.assignedValue;
      root.appendChild(ui.divV([
        ui.divH([
          sectionIcon(inconsistentIcon(consistency.restriction)),
          ui.divText(`Current value is inconsistent. Computed value was ${
            this.isScalar$.value ? consistentValue : 'different'
          }`),
        ]),
        ...(this.isDataFrame ?
          [indentedLink('Show dataframes diff', () => this.showDataFrameDiff(consistentValue))] : []),
        indentedLink('Reset to consistent value', () => this.consistencyReset()),
      ]));
    }

    return root;
  }

  private isLinked(consistency?: ConsistencyInfo) {
    return !!consistency && !consistency.inconsistent && isLinkedRestriction(consistency.restriction);
  }

  private createPopover() {
    injectPopoverStyles();
    const popover = ui.div([], `d4-tooltip ${POPOVER_CLASS}`);
    popover.popover = 'auto';
    Object.assign(popover.style, {
      fontStyle: 'normal',
      pointerEvents: 'all',
      maxWidth: '300px',
      overflowY: 'auto',
      cursor: 'default',
    });
    popover.style.setProperty('text-wrap', 'pretty');
    return popover;
  }

  private displayValidation() {
    const {currentIcon: icon, currentPopover: popover} = this;
    if (!icon || !popover)
      return;
    this.alignPopover(icon, popover);
    ui.empty(popover);
    popover.appendChild(this.renderValidationResults());
    popover.showPopover();
  }

  private alignPopover(target: HTMLElement, popover: HTMLElement): void {
    const bounds = target.getBoundingClientRect();
    const margin = POPOVER_VIEWPORT_MARGIN;
    const viewportHeight = window.innerHeight;
    const spaceBelow = viewportHeight - bounds.top - margin;
    const spaceAbove = bounds.bottom - margin;

    popover.style.inset = 'unset';
    popover.style.left = (bounds.x + 20) + 'px';

    if (spaceBelow >= POPOVER_MAX_HEIGHT) {
      // Enough room below — open downward from the icon.
      popover.style.top = bounds.top + 'px';
      popover.style.bottom = 'unset';
      popover.style.maxHeight = POPOVER_MAX_HEIGHT + 'px';
    } else if (spaceAbove >= POPOVER_MAX_HEIGHT) {
      // Not enough below, enough above — open upward from the icon.
      popover.style.top = 'unset';
      popover.style.bottom = (viewportHeight - bounds.bottom) + 'px';
      popover.style.maxHeight = POPOVER_MAX_HEIGHT + 'px';
    } else {
      // Neither side fits fully — anchor to the viewport top and take up to POPOVER_MAX_HEIGHT.
      popover.style.top = margin + 'px';
      popover.style.bottom = 'unset';
      popover.style.maxHeight = Math.min(POPOVER_MAX_HEIGHT, viewportHeight - 2 * margin) + 'px';
    }
  }
}

export interface ValidationIconT extends ValidationIcon {};
