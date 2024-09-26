import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

/** Lookup table string consts */
export enum LOOKUP {
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
export enum INPUTS_DF {
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

/** Load dataframe using the command */
async function loadTable(command: string): Promise<DG.DataFrame | null> {
  const funcCall = grok.functions.parse(command);
  
  if (!(funcCall instanceof DG.FuncCall)) {
    grok.shell.warning(`${LOOKUP.LOAD}, ${LOOKUP.FUNCTION}`);
    return null;
  }
  
  const calledFuncCall = await funcCall.call();
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
export async function getInputsTable(command: string): Promise<DG.DataFrame | null> {
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
export function getLookupsInfo(inputsLookup: string) {
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