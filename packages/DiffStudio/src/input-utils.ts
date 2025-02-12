import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {ANNOT_SEPAR, BRACE_CLOSE, BRACE_OPEN, BRACKET_CLOSE, BRACKET_OPEN, CONTROL_SEP, Input,
  ModelError, IVP} from '@datagrok/diff-studio-tools';

import {ERROR_MSG, INPUT_TYPE, LINK} from './ui-constants';

/** */
type InputsStore = {
  const: DG.InputBase,
  min: DG.InputBase,
  max: DG.InputBase,
  isFixed: boolean,
};

/** */
function getOptions(name: string, modelInput: Input, modelBlock: string): DG.PropertyOptions {
  const options: DG.PropertyOptions = {
    name: name,
    defaultValue: modelInput.value,
    type: DG.TYPE.FLOAT,
    inputType: INPUT_TYPE.FLOAT,
  };


  if (modelInput.annot !== null) {
    let annot = modelInput.annot;
    let descr: string | undefined = undefined;

    let posOpen = annot.indexOf(BRACKET_OPEN);
    let posClose = annot.indexOf(BRACKET_CLOSE);

    if (posOpen !== -1) {
      if (posClose === -1) {
        throw new ModelError(
          `${ERROR_MSG.MISSING_CLOSING_BRACKET}. Correct annotation in the **${modelBlock}** block.`,
          LINK.INTERFACE,
          annot,
        );
      }

      descr = annot.slice(posOpen + 1, posClose);

      annot = annot.slice(0, posOpen);
    }

    posOpen = annot.indexOf(BRACE_OPEN);
    posClose = annot.indexOf(BRACE_CLOSE);

    if (posOpen >= posClose) {
      throw new ModelError(
        `${ERROR_MSG.INCORRECT_BRACES_USE}. Correct annotation in the ***${modelBlock}** block.`,
        LINK.INTERFACE,
        annot,
      );
    }

    let pos: number;
    let key: string;
    let val;

    annot.slice(posOpen + 1, posClose).split(ANNOT_SEPAR).forEach((str) => {
      pos = str.indexOf(CONTROL_SEP);

      if (pos === -1) {
        throw new ModelError(
          `${ERROR_MSG.MISSING_COLON}. Correct annotation in the **${modelBlock}** block.`,
          LINK.INTERFACE,
          annot,
        );
      }

      key = str.slice(0, pos).trim();
      val = str.slice(pos + 1).trim();

      // @ts-ignore
      options[key !== 'caption' ? key : 'friendlyName'] = strToVal(val);
    });

    options.description = descr ?? '';
    options.name = options.friendlyName ?? options.name;
    options.friendlyName = options.name;
  }

  return options;
}; // getOptions

/** */
function getAnalysisInputs(ivp: IVP): Map<string, InputsStore[]> {
  const analysisInputs = new Map<string, InputsStore[]>();

  return analysisInputs;
} // getAnalysisInputs
