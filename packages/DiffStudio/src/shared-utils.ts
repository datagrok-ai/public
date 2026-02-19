import {MISC, TITLE, ERROR_MSG, LINK} from './ui-constants';
import {CONTROL_EXPR} from './constants';
import {CONTROL_SEP, BRACE_OPEN, BRACE_CLOSE, BRACKET_OPEN, BRACKET_CLOSE, ANNOT_SEPAR, Input, ARG_INPUT_KEYS_MAPPING} from './scripting-tools';

/** Diff Studio model error */
export class ModelError extends Error {
  private helpUrl: string;
  private toHighlight: string = undefined;

  constructor(message: string, helpUrl: string, toHighlight?: string) {
    super(message);
    this.helpUrl = helpUrl;
    this.toHighlight = toHighlight;
  }

  public getHelpUrl() {
    return this.helpUrl;
  }

  public getToHighlight() {
    return this.toHighlight;
  }
}; // ModelError


/**  String-to-value */
export const strToVal = (s: string) => {
  const num = Number(s);
  return !isNaN(num) ? num : s === 'true' ? true : s === 'false' ? false : s;
};


/** Copy of getLookupsInfo without UI warnings */
export function getLookupsInfoData(inputsLookup: string) {
  const info = new Map<string, string>();

  const braceOpenIdx = inputsLookup.indexOf(BRACE_OPEN);
  const braceCloseIdx = inputsLookup.indexOf(BRACE_CLOSE);

  if (braceOpenIdx < 0) {
    return null;
  }

  if (braceCloseIdx < 0) {
    return null;
  }

  // extract name
  info.set('name', inputsLookup.slice(0, braceOpenIdx).replaceAll(' ', ''));

  // extract features
  const options = inputsLookup.slice(braceOpenIdx + 1, braceCloseIdx).split(ANNOT_SEPAR);
  let sepIdx: number;

  for (const opt of options) {
    sepIdx = opt.indexOf(CONTROL_SEP);

    if (sepIdx < 0) {
      return null;
    }

    info.set(opt.slice(0, sepIdx).trim(), opt.slice(sepIdx + 1).trim());
  }

  // extract tooltip
  const bracketOpenIdx = inputsLookup.indexOf(BRACKET_OPEN);
  if (bracketOpenIdx > 0) {
    const bracketCloseIdx = inputsLookup.indexOf(BRACKET_CLOSE);

    if (bracketCloseIdx < 0) {
      return null;
    }

    info.set(MISC.TOOLTIP, inputsLookup.slice(bracketOpenIdx + 1, bracketCloseIdx));
  }

  if (info.get(MISC.CHOICES) === undefined) {
    return null;
  }

  return {
    name: info.get(MISC.NAME) ?? '',
    caption: info.get(MISC.CAPTION) ?? (info.get(MISC.NAME) ?? ''),
    category: info.get(MISC.CATEGORY) ?? TITLE.MISC,
    tooltip: info.get(MISC.TOOLTIP) ?? '',
    choices: info.get(MISC.CHOICES),
  };
} // getLookupsInfoData


/** Return options with respect to the model input specification */
export function getOptions(name: string, modelInput: Input, modelBlock: string, startingInputs?: Map<string, number>) {
  const options: Record<string,any> = {
    name: name,
    defaultValue: modelInput.value,
    type: 'double',
    inputType: 'Float',
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
    let val: string;

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

      options[key !== 'caption' ? key : 'friendlyName'] = strToVal(val);
    });

    options.description = descr ?? '';
    options.friendlyName = options.friendlyName ?? options.name;
    options.caption = options.friendlyName;
  }

  if (modelBlock === CONTROL_EXPR.ARG) {
    options.friendlyName = options.name;
    options.caption = options.name;
    options.name = ARG_INPUT_KEYS_MAPPING[options.name];
  }

  if (startingInputs) {
    options.defaultValue = startingInputs
      .get(options.name!.replace(' ', '').toLowerCase()) ?? options.defaultValue;
    modelInput.value = options.defaultValue;
  }

  return options;
}; // getOptions
