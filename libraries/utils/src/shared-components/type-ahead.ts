import * as ui from 'datagrok-api/ui';

import typeahead from 'typeahead-standalone'; // imports library (js)
import 'typeahead-standalone/dist/basic.css'; // imports basic styles (css)


export namespace u2 {
    export function typeAhead(source: string[], minLength: number = 1, limit: number = 5, highlight: boolean = false,
      autoSelect: boolean = false, hint: boolean = false, diacritics: boolean = false, debounceRemote: number = 100):
      HTMLDivElement {
      const inputElement = ui.searchInput('', '');

      typeahead({
        input: <HTMLInputElement>inputElement.input,
        source: {
          local: source,
          // prefetch: {...}
          // remote: {...}
        },
        minLength: minLength,
        limit: limit,
        highlight: highlight,
        autoSelect: autoSelect,
        hint: hint,
        diacritics: diacritics,
        debounceRemote: debounceRemote,
      });

      const element = inputElement.root.getElementsByClassName('typeahead-standalone')[0];
      return ui.div([element]);
    }
}
