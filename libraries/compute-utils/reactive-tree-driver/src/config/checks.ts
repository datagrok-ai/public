import * as DG from 'datagrok-api/dg';
import {RuleExpr, RuleValidatorEffect} from './PipelineConfiguration';

/** The annotation options the driver validates. Keys and values match the function
 *  annotation syntax; `check` links use the same object. */
export type CheckOptions = {
  nullable?: boolean;
  min?: number;
  max?: number;
  /** Regex literal, `/pattern/flags`. */
  validator?: string;
  choices?: any[];
  /** Column kind, column type or semantic type. */
  type?: string;
  semType?: string;
  /** Id of the dataframe io the column belongs to. */
  table?: string;
  allowNulls?: boolean;
};

export type CheckSeverity = 'error' | 'warning' | 'notification';

export type CheckKey = 'required' | 'min' | 'max' | 'validator' | 'choices' | 'type' | 'semType' | 'table' | 'allowNulls';

export type ExpandedCheck = {
  key: CheckKey;
  needsTable: boolean;
  params: {when: RuleExpr, effects: RuleValidatorEffect[]};
};

export type CheckExtras = {
  when?: RuleExpr;
  message?: RuleExpr;
  severity?: CheckSeverity;
};

export const checkOptionKeys: (keyof CheckOptions)[] =
  ['nullable', 'min', 'max', 'validator', 'choices', 'type', 'semType', 'table', 'allowNulls'];

// aliases shared by annotation-derived and config checks
export const VALUE = 'value';
export const TABLE = 'table';
export const TARGET = 'target';

const present = {'!': {missing: [VALUE]}};
const value = {var: VALUE};

const regexLiteral = /^\/(.*)\/([a-z]*)$/s;

export function parseRegexLiteral(literal: string): {pattern: string, flags: string} | undefined {
  const match = regexLiteral.exec(literal);
  return match ? {pattern: match[1], flags: match[2]} : undefined;
}

function conditions(options: CheckOptions): {key: CheckKey, needsTable: boolean, when: RuleExpr, message: string}[] {
  const out: {key: CheckKey, needsTable: boolean, when: RuleExpr, message: string}[] = [];
  const add = (key: CheckKey, when: RuleExpr, message: string, needsTable = false) =>
    out.push({key, needsTable, when: {and: [present, when]}, message});
  if (options.nullable === false)
    out.push({key: 'required', needsTable: false, when: {missing: [VALUE]}, message: 'Missing value'});
  if (options.min != null)
    add('min', {'<': [value, options.min]}, `Must be at least ${options.min}`);
  if (options.max != null)
    add('max', {'>': [value, options.max]}, `Must be at most ${options.max}`);
  if (options.validator != null) {
    const regex = parseRegexLiteral(options.validator)!;
    add('validator', {'!': {regex: [value, regex.pattern, regex.flags]}}, `Must match ${options.validator}`);
  }
  if (options.choices != null)
    add('choices', {'!': {in: [value, options.choices]}}, `Must be one of: ${options.choices.join(', ')}`);
  if (options.type != null)
    add('type', {'!': {columnIs: [value, options.type]}}, `Column must be ${options.type}`);
  if (options.semType != null)
    add('semType', {'!': {columnIs: [value, options.semType]}}, `Column must have semantic type ${options.semType}`);
  if (options.table != null) {
    add('table', {and: [{'!': {missing: [TABLE]}}, {'!': {in: [{var: `${VALUE}.name`}, {columns: [{var: TABLE}]}]}}]},
      'Column does not belong to the table', true);
  }
  if (options.allowNulls === false)
    add('allowNulls', {'>': [{nulls: value}, 0]}, 'Column has missing values');
  return out;
}

export function validateCheckOptions(id: string, options: CheckOptions) {
  for (const key of Object.keys(options)) {
    if (!checkOptionKeys.includes(key as keyof CheckOptions))
      throw new Error(`Check ${id}: unknown option ${key}`);
  }
  if (options.validator != null && !parseRegexLiteral(options.validator))
    throw new Error(`Check ${id}: validator must be a regex literal /pattern/flags`);
  if (options.choices != null && !Array.isArray(options.choices))
    throw new Error(`Check ${id}: choices must be an array`);
  for (const key of ['min', 'max'] as const) {
    if (options[key] != null && typeof options[key] !== 'number')
      throw new Error(`Check ${id}: ${key} must be a number`);
  }
}

/** Expands annotation-style options into one validator per option, in the rule params shape.
 *  Values are read from the `value` alias, the table from `table`, results go to `target`. */
export function expandChecks(options: CheckOptions, extras: CheckExtras = {}): ExpandedCheck[] {
  return conditions(options).map(({key, needsTable, when, message}) => ({
    key,
    needsTable,
    params: {
      when: extras.when == null ? when : {and: [extras.when, when]},
      effects: [{effect: extras.severity ?? 'error', targets: [TARGET], message: extras.message ?? message}],
    },
  }));
}

function parseBool(val: string | undefined): boolean | undefined {
  if (val === 'true') return true;
  if (val === 'false') return false;
  return undefined;
}

function parseNumber(val: string | undefined): number | undefined {
  if (val == null || val === '')
    return undefined;
  const num = Number(val);
  return Number.isNaN(num) ? undefined : num;
}

function parseChoices(val: string | undefined): any[] | undefined {
  if (val == null)
    return undefined;
  try {
    const parsed = JSON.parse(val);
    return Array.isArray(parsed) ? parsed : undefined;
  } catch {
    return undefined;
  }
}

/** Reads the supported options off a function parameter. Function-backed choices, GrokScript
 *  validators and other unsupported forms are left out. */
export function parseAnnotationChecks(prop: DG.Property): CheckOptions {
  const options = prop.options ?? {};
  const checks: CheckOptions = {};
  const isNumeric = prop.propertyType === DG.TYPE.INT || prop.propertyType === DG.TYPE.FLOAT ||
    prop.propertyType === DG.TYPE.BIG_INT;
  const isColumn = prop.propertyType === DG.TYPE.COLUMN;
  if (isNumeric) {
    const min = parseNumber(options.min);
    const max = parseNumber(options.max);
    if (min != null) checks.min = min;
    if (max != null) checks.max = max;
  }
  if (options.validator != null && parseRegexLiteral(options.validator))
    checks.validator = options.validator;
  const choices = parseChoices(options.choices);
  if (choices && DG.TYPES_SCALAR.has(prop.propertyType))
    checks.choices = choices;
  if (isColumn) {
    // the platform stores a column `type:` annotation as the `columns` option and the type filter
    const type = options.type ?? options.columns ?? prop.columnTypeFilter;
    if (type) checks.type = type;
    const semType = options.semType || prop.semType;
    if (semType) checks.semType = semType;
    if (options.table) checks.table = options.table;
    if (parseBool(options.allowNulls) === false) checks.allowNulls = false;
  }
  return checks;
}

export function isOptionalAnnotation(prop: DG.Property): boolean {
  return prop.options?.optional === 'true' || prop.options?.nullable === 'true';
}
