import * as ui from 'datagrok-api/ui';

import typeahead from 'typeahead-standalone'; // imports library (js)
import '../../css/typeahead-input.css';

export namespace u2 {
  export function typeAhead(label: string = '', source: string[] = [], minLength: number = 1, limit: number = 5):
    HTMLDivElement {
    const inputElement = ui.stringInput(label, '');
    (inputElement.input as HTMLInputElement).placeholder = 'Search';

    typeahead({
      input: <HTMLInputElement>inputElement.input,
      source: {
        local: source,
      },
      minLength: minLength,
      limit: limit,
      hint: false,
    });

    // set classes
    inputElement.root.getElementsByClassName('tt-list')[0].className = 'ui-input-list';
    inputElement.root.getElementsByClassName('tt-input')[0].className = 'ui-input-editor';

    return inputElement.root as HTMLDivElement;
  }
}
