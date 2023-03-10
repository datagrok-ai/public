import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import typeahead from 'typeahead-standalone'; // imports library (js)
import '../../css/typeahead-input.css';

export namespace u2 {
  export class TypeAhead extends DG.InputBase {
    constructor(label: string = '', source: string[] = [],
      options?: {minLength?: number, limit?: number}) {
      const inputElement = ui.stringInput(label, '');
      super(inputElement.dart);

      (this.input as HTMLInputElement).placeholder = 'Search';

      typeahead({
        input: <HTMLInputElement> this.input,
        source: {
          local: source,
        },
        minLength: options?.minLength ?? 1,
        limit: options?.limit ?? 5,
        hint: false,
      });

      this.root.getElementsByClassName('tt-list')[0].className = 'ui-input-list';
      this.root.getElementsByClassName('tt-input')[0].className = 'ui-input-editor';
    }
  }
}
