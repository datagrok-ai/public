import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import typeahead from 'typeahead-standalone';
import {Dictionary, typeaheadConfig} from 'typeahead-standalone/dist/types';
import '../../css/typeahead-input.css';


type TypeAheadConfig = Omit<typeaheadConfig<Dictionary>, 'input' | 'className'>

export namespace u2 {
  export class TypeAhead extends DG.InputBase {
    constructor(config: TypeAheadConfig) {
      const inputElement = ui.stringInput('', '');
      super(inputElement.dart);

      const typeAheadConfig: typeaheadConfig<Dictionary> = Object.assign(
        {input: <HTMLInputElement> this.input}, config);

      typeahead(typeAheadConfig);

      this.root.getElementsByClassName('tt-list')[0].className = 'ui-input-list';
      this.root.getElementsByClassName('tt-input')[0].className = 'ui-input-editor';
    }
  }
}
