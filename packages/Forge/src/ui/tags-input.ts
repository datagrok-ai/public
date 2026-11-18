import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {stringsOf} from '../preparation/preparation-options';

const TOOLTIP = 'Free tags: type one and press Enter. A comma separates two tags.';
const PLACEHOLDER = 'Type a tag and press Enter';

/** The **Tags** input: the platform's chips input, starting with [tags]. */
export function tagsInput(tags: string[]): DG.InputBase {
  const input = ui.input.tags<string>('Tags', {value: tags, allowNew: true, multiValue: true, tooltipText: TOOLTIP});
  // The chips input has no placeholder option; its text box shows the platform's "Type to search..." otherwise. Nor
  // does it turn the browser's autofill off, as the platform's text inputs do: a value picked there makes no chip.
  const box = input.root.querySelector('input.d4-tags-selector-input');
  if (box instanceof HTMLInputElement) {
    box.placeholder = PLACEHOLDER;
    box.autocomplete = 'off';
  }
  return input;
}

export function tagsOfInput(input: DG.InputBase): string[] {
  return stringsOf(input.value);
}
