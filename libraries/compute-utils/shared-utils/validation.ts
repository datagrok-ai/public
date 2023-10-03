import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

// validation/advisory system
export interface ActionItems {
    actionName: string;
    action: Function;
  }

export interface Advice {
    description: string;
    actions?: ActionItems[];
  }

export interface ValidationResultBase {
    // awaiting for validation results
    pending?: boolean;
    errors?: Advice[];
    warnings?: Advice[];
    notifications?: Advice[];
  }

export interface ValidationResult extends ValidationResultBase {
    // revalidation request
    revalidate?: string[];
    // revalidations context
    context?: any;
  }

export function isValidationPassed(result?: ValidationResult) {
  return !result?.errors?.length && !result?.pending;
}

export function makeAdvice(description: string, actions?: ActionItems[]) {
  return {description, actions};
}

export function getErrorMessage(result?: ValidationResult) {
  if (result?.errors)
    return result.errors.map((err) => err.description).join('; ');
}

export interface ValidationPayload {
    errors?: (string | Advice)[],
    warnings?: (string | Advice)[],
    notifications?: (string | Advice)[],
  }

export function makeValidationResult(payload?: ValidationPayload): ValidationResultBase {
  const wrapper = (item: string | Advice) => typeof item === 'string' ? makeAdvice(item) : item;
  return {
    errors: payload?.errors?.map((err) => wrapper(err)),
    warnings: payload?.warnings?.map((warn) => wrapper(warn)),
    notifications: payload?.notifications?.map((note) => wrapper(note)),
  };
}

export function makePendingValidationResult(): ValidationResult {
  return {pending: true};
}

export function makeRevalidation(revalidate: string[], context?: any, result?: ValidationResultBase): ValidationResult {
  return {revalidate, context, ...result};
}

export interface ValidationInfo {
    param: string,
    funcCall: DG.FuncCall,
    lastCall?:DG.FuncCall,
    isRevalidation: boolean,
    isNewOutput: boolean,
    signal: AbortSignal,
    context?: any
  }

export type Validator = (val: any, info: ValidationInfo)
    => Promise<ValidationResult | undefined>;

export type ValidatorFactory = (params: any) => { validator: Validator };

export const nonNullValidator: Validator = async (value: any) => {
  if (value == null)
    return makeValidationResult({errors: ['Missing value']});
};

export function getValidationIcon(messages: ValidationResultBase | undefined) {
  let popover: any;
  let icon: any;
  if (messages?.pending)
    icon = ui.iconFA('spinner', () => {displayValidation(messages, icon, popover);});

  if (messages?.errors) {
    icon = ui.iconFA('exclamation', () => {displayValidation(messages, icon, popover);});
    icon.style.color = 'var(--red-3)';
  } else if (messages?.warnings) {
    icon = ui.iconFA('exclamation', () => {displayValidation(messages, icon, popover);});
    icon.style.color = 'var(--orange-2)';
  } else if (messages?.notifications) {
    icon = ui.iconFA('info', () => {displayValidation(messages, icon, popover);} );
    icon.style.color = 'var(--blue-1)';
  }
  if (icon)
    popover = addPopover(icon);

  return icon;
}

function addPopover(icon: HTMLElement) {
  const popover = ui.div();
  stylePopover(popover);
  icon.appendChild(popover);
  return popover;
}

function displayValidation(messages: ValidationResultBase, icon: HTMLElement, popover: HTMLElement) {
  if (popover && icon) {
    alignPopover(icon, popover);
    while (popover.firstChild && popover.removeChild(popover.firstChild));
    const content = renderDynamicHelp(messages);
    popover.appendChild(content);
    popover.showPopover();
  }
}

function alignPopover(target: HTMLElement, popover: HTMLElement): void {
  const bounds = target.getBoundingClientRect().toJSON();
  popover.style.inset = 'unset';
  popover.style.top = bounds.y + 'px';
  popover.style.left = (bounds.x + 20) + 'px';
}

function stylePopover(popover: HTMLElement): void {
  popover.popover = 'auto';
  popover.style.cursor = 'default';
  popover.style.padding = '10px';
  popover.style.background = '#fdffe5';
  popover.style.border = '1px solid #E4E6CE';
  popover.style.borderRadius = '2px';
  popover.style.boxShadow = '0 0 5px #E4E6CE';
  popover.style.maxWidth = '500px';
}

function renderDynamicHelp(messages: ValidationResultBase) {
  const root = ui.div('', {style: {display: 'flex', flexDirection: 'column'}});
  for (const [category, advices] of Object.entries(messages)) {
    const icon = getAdviceIcon(category);
    if (!icon || !advices)
      continue;
    for (const advice of advices as Advice[]) {
      root.appendChild(ui.span([
        icon,
        advice.description,
        ...(advice.actions ?? []).map(
          (action) => ui.link(action.actionName, action.action, undefined, {style: {paddingLeft: '10px'}})),
      ], {style: {lineHeight: '1.2', marginBottom: '10px'}}));
    }
  }
  return root;
}

function getAdviceIcon(category: string) {
  let icon: HTMLElement | undefined;
  if (category === 'errors') {
    icon = ui.iconFA('exclamation');
    icon.style.color = 'var(--red-3)';
  } else if (category === 'warnings') {
    icon = ui.iconFA('exclamation');
    icon.style.color = 'var(--orange-2)';
  } else if (category === 'notifications') {
    icon = ui.iconFA('info');
    icon.style.color = 'var(--blue-1)';
  }
  if (icon) {
    icon.style.cursor = 'default';
    icon.style.display = 'inline-block';
    icon.style.fontFamily = '"Font Awesome 5 Pro"';
  }
  return icon;
}

