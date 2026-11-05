import * as ui from 'datagrok-api/ui';

/** Shows an engine's description in a dialog: its first sentence as text, the whole
 *  description (implementation notes, credits) in the tooltip, so long descriptions
 *  neither widen nor crowd the dialog. */
export function setEngineDescription(element: HTMLElement, description: string): void {
  const text = description.trim();
  const first = /^(.+?[.!?])\s+(?=[A-Z])/.exec(text)?.[1] ?? text;
  element.textContent = first;
  element.dataset.description = text;
  element.style.cursor = first !== text ? 'help' : '';
}

/** Binds the description tooltip once; it follows {@link setEngineDescription}. */
export function bindEngineDescription(element: HTMLElement): void {
  ui.tooltip.bind(element, () => element.dataset.description && element.dataset.description !== element.textContent ?
    ui.divText(element.dataset.description, {style: {maxWidth: '420px'}}) : null);
}
