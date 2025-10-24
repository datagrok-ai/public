/* eslint-disable valid-jsdoc */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

/** Lookup table string consts */
enum LOOKUP {
  LOAD = 'Failed to load lookup table: ',
  PLATFORM = 'the platform issue',
  FUNCTION = 'incorrect function',
  NO_DF = 'no dataframe',
  INCORRECT = 'incorrect dataframe, ',
  ROWS = 'at least one row is needed.',
  NULLS = 'missing values are not allowed.',
  NUMS = 'no numerical columns.',
  CHOICES = 'first column must contain strings.',
  DEFAULT = 'Default',
};

/** Inputs table constants */
enum INPUTS_DF {
  MIN_ROWS_COUNT = 1,
  INP_NAMES_IDX = 0,
  INPUT_SETS_COL_IDX = 0,
};

/** Annotating constants */
enum ANNOT {
  CHOICES = 'choices',
  NAME = 'name',
  CATEGORY = 'category',
  CAPTION = 'caption',
  TOOLTIP = 'tooltip',
  MISC = 'Misc'
};

type Format = string | null | undefined;

/** Lookup data specification */
type LookupData = {
  arr: Int32Array | Uint32Array | Float32Array | Float64Array,
  format: Format,
};

/** Return value formatted with respect to the given format */
export function getFormatted(value: number, format: Format): number {
  if ((format === null) || (format === undefined))
    return value;

  if (format.includes('E') || format.includes('e'))
    return value;

  const precision = format.length - 2;

  if ((precision < 1) || (precision > 20))
    return value;

  return Number(value.toPrecision(precision));
}

/** Load dataframe using the command */
async function loadTable(command: string): Promise<DG.DataFrame | null> {
  const funcCall = grok.functions.parse(command);

  if (!(funcCall instanceof DG.FuncCall)) {
    grok.shell.warning(`${LOOKUP.LOAD}, ${LOOKUP.FUNCTION}`);
    return null;
  }

  const calledFuncCall = await funcCall.call(undefined, undefined, {processed:true, report:false});
  const output = calledFuncCall.getOutputParamValue();

  if (!(output instanceof DG.DataFrame)) {
    grok.shell.warning(`${LOOKUP.LOAD}, ${LOOKUP.NO_DF}`);
    return null;
  }

  return output;
}

/** Check correctness of lookup table */
function isLookupTableCorrect(table: DG.DataFrame | null): boolean {
  if (table === null)
    return false;

  const cols = table.columns;

  // check rows count
  if (table.rowCount < INPUTS_DF.MIN_ROWS_COUNT) {
    grok.shell.warning(`${LOOKUP.LOAD}${LOOKUP.INCORRECT}${LOOKUP.ROWS}`);
    return false;
  }

  // check nulls & numerical cols
  let numColsCount = 0;

  for (const col of cols) {
    if (col.stats.missingValueCount > 0) {
      grok.shell.warning(`${LOOKUP.LOAD}${LOOKUP.INCORRECT}${LOOKUP.NULLS}`);
      return false;
    }

    if (col.isNumerical)
      ++numColsCount;
  }

  if (numColsCount === 0) {
    grok.shell.warning(`${LOOKUP.LOAD}${LOOKUP.INCORRECT}${LOOKUP.NUMS}`);
    return false;
  }

  // check column with names of inputs
  if (cols.byIndex(INPUTS_DF.INP_NAMES_IDX).type !== DG.COLUMN_TYPE.STRING) {
    grok.shell.warning(`${LOOKUP.LOAD}${LOOKUP.INCORRECT}${LOOKUP.CHOICES}`);
    return false;
  }

  return true;
} // isLookupTableCorrect

/** Return table with inputs */
async function getInputsTable(command: string): Promise<DG.DataFrame | null> {
  try {
    const table = await loadTable(command);

    if (isLookupTableCorrect(table))
      return table;
  } catch (err) {
    const msg = (err instanceof Error) ? err.message : `check ${command}`;
    grok.shell.warning(`${LOOKUP.LOAD} ${msg}`);
  }

  return null;
}

/** Return specification of lookup table input */
function getLookupsInfo(inputsLookup: string) {
  const info = new Map<string, string>();

  const braceOpenIdx = inputsLookup.indexOf('{');
  const braceCloseIdx = inputsLookup.indexOf('}');

  if (braceOpenIdx < 0) {
    grok.shell.warning('Missing "{"');
    return null;
  }

  if (braceCloseIdx < 0) {
    grok.shell.warning('Missing "}"');
    return null;
  }

  // extract name
  info.set(ANNOT.NAME, inputsLookup.slice(0, braceOpenIdx).replaceAll(' ', ''));

  // extract features
  const options = inputsLookup.slice(braceOpenIdx + 1, braceCloseIdx).split(';');
  let sepIdx: number;

  for (const opt of options) {
    sepIdx = opt.indexOf(':');

    if (sepIdx < 0) {
      grok.shell.warning('Missing ":"');
      return null;
    }

    info.set(opt.slice(0, sepIdx).trim(), opt.slice(sepIdx + 1).trim());
  }

  // extract tooltip
  const bracketOpenIdx = inputsLookup.indexOf('[');
  if (bracketOpenIdx > 0) {
    const bracketCloseIdx = inputsLookup.indexOf(']');

    if (bracketCloseIdx < 0) {
      grok.shell.warning('Missing "]"');
      return null;
    }

    info.set(ANNOT.TOOLTIP, inputsLookup.slice(bracketOpenIdx + 1, bracketCloseIdx));
  }

  if (info.get(ANNOT.CHOICES) === undefined) {
    grok.shell.warning(`Missing "${ANNOT.CHOICES}"-expression`);
    return null;
  }

  return {
    name: info.get(ANNOT.NAME) ?? '',
    caption: info.get(ANNOT.CAPTION) ?? (info.get(ANNOT.NAME) ?? ''),
    category: info.get(ANNOT.CATEGORY) ?? ANNOT.MISC,
    tooltip: info.get(ANNOT.TOOLTIP) ?? '',
    choices: info.get(ANNOT.CHOICES),
  };
} // getLookupsInfo

/** Return values lookup choice input */
export async function getLookupChoiceInput(inputsLookup: string, constIputs: Map<string, DG.InputBase>) {
  if (inputsLookup === undefined)
    return null;

  const info = getLookupsInfo(inputsLookup);

  if ((info === null) || (info.choices === undefined))
    return null;

  const inputsDf = await getInputsTable(info.choices);

  if (inputsDf === null)
    return null;

  const cols = inputsDf.columns;
  const rowCount = inputsDf.rowCount;

  const inpSetsNames = cols.byIndex(INPUTS_DF.INPUT_SETS_COL_IDX).toList();
  const choices = [LOOKUP.DEFAULT as string].concat(inpSetsNames);

  const defaultInputs = new Map<string, any>();

  constIputs.forEach((input, name) => defaultInputs.set(name, input.value));

  const tableInputs = new Map<string, Map<string, number>>(); // set <-> {(input <-> value)}
  const colsRaw = new Map<string, LookupData>();

  for (const col of cols) {
    if (col.isNumerical) {
      colsRaw.set(col.name, {
        arr: col.getRawData(),
        format: col.meta.format,
      });
    }
  }

  for (let row = 0; row < rowCount; ++row) {
    const inputs = new Map<string, number>();
    colsRaw.forEach((info, name) => inputs.set(name, getFormatted(info.arr[row], info.format)));
    tableInputs.set(inpSetsNames[row], inputs);
  }

  // create input for lookup table use
  const lookupChoiceInput = ui.input.choice<string>(info.caption, {
    items: choices,
    nullable: false,
    value: choices[0],
    tooltipText: info.tooltip,
    onValueChanged: (value) => {
      if (value === LOOKUP.DEFAULT)
        constIputs.forEach((input, name) => input.value = defaultInputs.get(name));
      else {
        const colInputs = tableInputs.get(value);
        constIputs.forEach((input, name) => input.value = colInputs!.get(name) ?? input.value);
      }
    },
  });

  return {
    input: lookupChoiceInput,
    category: info.category ?? 'Misc',
  };
}
