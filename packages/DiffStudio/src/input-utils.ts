import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {ANNOT_SEPAR, BRACE_CLOSE, BRACE_OPEN, BRACKET_CLOSE, BRACKET_OPEN, CONTROL_SEP, Input,
  ModelError, IVP, ARG_INPUT_KEYS, CONTROL_EXPR, SCRIPTING} from '@datagrok/diff-studio-tools';

import {ERROR_MSG, INPUT_TYPE, LINK, TITLE} from './ui-constants';
import {strToVal} from './utils';

/** Store of analysis inputs */
type InputsStore = {
  const: DG.InputBase,
  min: DG.InputBase,
  max: DG.InputBase,
  isFixed: boolean,
};

/** Set tooltip of inputs */
function setTooltip(inps: InputsStore, hint: string): void {
  inps.const.setTooltip(hint);
  const caption = inps.const.caption;
  inps.min.setTooltip(`Min value of ${caption}`);
  inps.max.setTooltip(`Max value of ${caption}`);
}

/** Returns options of the Diff Studio model input */
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

      options[key !== 'caption' ? key : 'friendlyName'] = strToVal(val);
    });

    options.description = descr ?? '';
    options.name = options.friendlyName ?? options.name;
    options.friendlyName = options.name;
  }

  return options;
}; // getOptions

/** Pull input to appropriate category & add tooltip */
function categorizeInput(options: DG.PropertyOptions, inps: InputsStore,
  inputsByCategories: Map<string, InputsStore[]>): void {
  const category = options.category;

  if (category === undefined)
    inputsByCategories.get(TITLE.MISC)?.push(inps);
  else if (inputsByCategories.has(category))
        inputsByCategories.get(category)!.push(inps);
  else
    inputsByCategories.set(category, [inps]);

  setTooltip(inps, options.description!);
};

/** Return inputs store created using options */
function inpStoreFromOpts(options: DG.PropertyOptions): InputsStore {
  const constInp = ui.input.forProperty(DG.Property.fromOptions(options));

  const minInp = ui.input.forProperty(DG.Property.fromOptions(options));
  minInp.caption = `${minInp.caption} (min)`;

  const maxInp = ui.input.forProperty(DG.Property.fromOptions(options));
  maxInp.caption = `${maxInp.caption} (max)`;

  return {
    const: constInp,
    min: minInp,
    max: maxInp,
    isFixed: true,
  };
}

/** Return inputs for analysis DS model */
function getAnalysisInputs(ivp: IVP): Map<string, InputsStore[]> {
  const inputsByCategories = new Map<string, InputsStore[]>();
  inputsByCategories.set(TITLE.MISC, []);
  let options: DG.PropertyOptions;

  // Inputs for argument
  for (const key of ARG_INPUT_KEYS) {
    //@ts-ignore
    options = getOptions(key, ivp.arg[key], CONTROL_EXPR.ARG);
    categorizeInput(options, inpStoreFromOpts(options), inputsByCategories);
  } // arg

  // Inputs for initial values
  ivp.inits.forEach((val, key) => {
    options = getOptions(key, val, CONTROL_EXPR.INITS);
    categorizeInput(options, inpStoreFromOpts(options), inputsByCategories);
  }); // inits

  // Inputs for parameters
  if (ivp.params !== null) {
    ivp.params.forEach((val, key) => {
      options = getOptions(key, val, CONTROL_EXPR.PARAMS);
      categorizeInput(options, inpStoreFromOpts(options), inputsByCategories);
    });
  } // params

  // Inputs for loop
  if (ivp.loop !== null) {
    options = getOptions(SCRIPTING.COUNT, ivp.loop.count, CONTROL_EXPR.LOOP);
    options.inputType = INPUT_TYPE.INT; // since it's an integer
    options.type = DG.TYPE.INT; // since it's an integer
    categorizeInput(options, inpStoreFromOpts(options), inputsByCategories);
  }

  return inputsByCategories;
} // getAnalysisInputs

/** */
function getAnalysisInpGroup(inps: InputsStore): HTMLElement[] {
  const roots: HTMLElement[] = [];

  const conInp = inps.const;
  const minInp = inps.min;
  const maxInp = inps.max;

  const showInps = (isFixed: boolean) => {
    conInp.root.hidden = !isFixed;
    minInp.root.hidden = isFixed;
    maxInp.root.hidden = isFixed;
  };

  const turnOnInp = ui.input.toggle(' ', {
    value: false,
    onValueChanged: (val) => {
      if (val) {
        inps.isFixed = !val;
        showInps(!val);
        turnOffInp.value = true;
      }
    },
  });
  turnOnInp.captionLabel.hidden = true;

  const turnOffInp = ui.input.toggle(' ', {
    value: false,
    onValueChanged: (val) => {
      if (!val) {
        inps.isFixed = !val;
        showInps(!val);
        turnOnInp.value = false;
      }
    },
  });
  turnOffInp.captionLabel.hidden = true;

  const getSwitchMock = ui.div([], 'sa-switch-input');

  conInp.root.insertBefore(turnOnInp.root, conInp.captionLabel);
  minInp.root.insertBefore(turnOffInp.root, minInp.captionLabel);
  maxInp.root.insertBefore(getSwitchMock, maxInp.captionLabel);

  roots.push(conInp.root);
  roots.push(minInp.root);
  roots.push(maxInp.root);

  showInps(true);

  return roots;
}

/** Return inputs form */
function getAnalysisForm(header: string, inputsByCategories: Map<string, InputsStore[]>,
  topCategory: string | null = null): HTMLElement {
  //const form = ui.div([ui.h1(header)], {style: {overflowY: 'scroll', width: '100%'}});

  const form = ui.form([]);
  form.append(ui.h1(header));

  const miscInputs = inputsByCategories.get(TITLE.MISC);

  if (inputsByCategories.size === 1) {
    miscInputs.forEach((inps) => {
      getAnalysisInpGroup(inps).forEach((root) => form.append(root));
    });
  } else {
    if (topCategory !== null) {
      form.append(ui.h2(topCategory));

      inputsByCategories.get(topCategory).forEach((inps) => {
        getAnalysisInpGroup(inps).forEach((root) => form.append(root));
      });
    }

    inputsByCategories.forEach((inps, category) => {
      if ((category !== TITLE.MISC) && (category !== topCategory)) {
        form.append(ui.h2(category));
        inps.forEach((inp) => {
          getAnalysisInpGroup(inp).forEach((root) => form.append(root));
        });
      }
    });

    if ((miscInputs.length > 0) && (topCategory !== TITLE.MISC)) {
      form.append(ui.h2(TITLE.MISC));
      miscInputs.forEach((inp) => {
        getAnalysisInpGroup(inp).forEach((root) => form.append(root));
      });
    }
  }

  form.style.overflowY = 'show';

  //form.classList.add('ui-form');

  return form;
} // getAnalysisForm

export function getForm(ivp: IVP): HTMLElement {
  const inps = getAnalysisInputs(ivp);

  return getAnalysisForm('Fit', inps);
}
