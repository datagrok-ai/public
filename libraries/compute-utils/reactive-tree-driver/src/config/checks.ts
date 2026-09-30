import * as DG from 'datagrok-api/dg';
import {RuleEffect, RuleExpr, RuleLogic, RuleSource, RuleValidatorEffect} from './PipelineConfiguration';

/** The annotation options the driver validates. Keys and values match the function
 *  annotation syntax; `check` links use the same object. */
export type CheckOptions = {
  nullable?: boolean;
  /** Synonym of `nullable`, as in annotations. */
  optional?: boolean;
  min?: number;
  max?: number;
  /** A regex literal `/pattern/flags`, or a GrokScript expression over the step's inputs and `value`. */
  validator?: string;
  /** GrokScript expression; the input is hidden while it is `false`. */
  visible?: string;
  /** Named validator functions, evaluated by the platform (1.28+). */
  validators?: string[];
  choices?: any[];
  /** Column kind, column type or semantic type. */
  type?: string;
  semType?: string;
  /** Id of the dataframe io the column belongs to. */
  table?: string;
  allowNulls?: boolean;
};

export type CheckSeverity = 'error' | 'warning' | 'notification';

export type CheckKey =
  'required' | 'min' | 'max' | 'validator' | 'validators' | 'visible' | 'choices' | 'type' | 'semType' | 'table' | 'allowNulls';

export type CheckFamily = 'validator' | 'meta';

export type ExpandedCheck = {
  key: CheckKey;
  family: CheckFamily;
  needsTable: boolean;
  /** The check reads the node's FuncCall through the `call` alias. */
  needsCall: boolean;
  /** The check's expression reads the step's other inputs by name. */
  needsInputs: boolean;
  params: {when: RuleLogic, effects: RuleEffect[], sources?: Record<string, RuleSource>};
};

export type CheckExtras = {
  when?: RuleLogic;
  message?: RuleExpr;
  severity?: CheckSeverity;
};

export const checkOptionKeys: (keyof CheckOptions)[] =
  ['nullable', 'optional', 'min', 'max', 'validator', 'validators', 'visible', 'choices', 'type', 'semType', 'table', 'allowNulls'];

// aliases shared by annotation-derived and config checks
export const VALUE = 'value';
export const TABLE = 'table';
export const TARGET = 'target';
export const CALL = 'call';
const VERDICTS = 'verdicts';

const present = {'!': {missing: [VALUE]}};
const value = {var: VALUE};

const regexLiteral = /^\/(.*)\/([a-z]*)$/s;

export function parseRegexLiteral(literal: string): {pattern: string, flags: string} | undefined {
  const match = regexLiteral.exec(literal);
  return match ? {pattern: match[1], flags: match[2]} : undefined;
}

type Condition = {
  key: CheckKey;
  family?: CheckFamily;
  needsTable: boolean;
  needsCall: boolean;
  needsInputs?: boolean;
  when: RuleLogic;
  /** Ready-made effects for non-validator families. */
  effects?: RuleEffect[];
  /** A fixed message, or the alias of a source whose verdicts carry the messages. */
  message?: string;
  verdicts?: string;
  /** An expression producing the message, overriding `message`. */
  verdictMessage?: RuleExpr;
  sources?: Record<string, RuleSource>;
};

function conditions(options: CheckOptions): Condition[] {
  const out: Condition[] = [];
  const add = (key: CheckKey, when: RuleLogic, message: string, needsTable = false) =>
    out.push({key, needsTable, needsCall: false, when: {and: [present, when]}, message});
  if (options.nullable === false)
    out.push({key: 'required', needsTable: false, needsCall: false, when: {missing: [VALUE]}, message: 'Missing value'});
  if (options.min != null)
    add('min', {'<': [value, options.min]}, `Must be at least ${options.min}`);
  if (options.max != null)
    add('max', {'>': [value, options.max]}, `Must be at most ${options.max}`);
  if (options.validator != null) {
    const regex = parseRegexLiteral(options.validator);
    if (regex)
      add('validator', {'!': {regex: [value, regex.pattern, regex.flags]}}, `Must match ${options.validator}`);
    else {
      out.push({
        key: 'validator', needsTable: false, needsCall: false, needsInputs: true,
        when: {and: [present, {'!!': {scriptVerdict: options.validator}}]},
        message: options.validator, verdictMessage: {scriptVerdict: options.validator},
      });
    }
  }
  if (options.visible != null) {
    out.push({
      key: 'visible', family: 'meta', needsTable: false, needsCall: false, needsInputs: true,
      when: {'==': [{script: options.visible}, false]},
      effects: [{effect: 'hide', targets: [TARGET]}],
    });
  }
  if (options.validators?.length) {
    out.push({
      key: 'validators', needsTable: false, needsCall: false,
      sources: {[VERDICTS]: {validators: {input: VALUE, names: options.validators}}},
      when: present, verdicts: VERDICTS,
    });
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
  if (options.optional != null && options.nullable != null && options.optional !== options.nullable)
    throw new Error(`Check ${id}: nullable and optional disagree`);
  for (const key of ['validator', 'visible'] as const) {
    if (options[key] != null && typeof options[key] !== 'string')
      throw new Error(`Check ${id}: ${key} must be a regex literal or a GrokScript expression`);
  }
  if (options.choices != null && !Array.isArray(options.choices))
    throw new Error(`Check ${id}: choices must be an array`);
  if (options.validators != null &&
      (!Array.isArray(options.validators) || options.validators.some((name) => typeof name !== 'string')))
    throw new Error(`Check ${id}: validators must be an array of function names`);
  for (const key of ['min', 'max'] as const) {
    if (options[key] != null && typeof options[key] !== 'number')
      throw new Error(`Check ${id}: ${key} must be a number`);
  }
}

/** Expands annotation-style options into one validator per option, in the rule params shape.
 *  Values are read from the `value` alias, the table from `table`, results go to `target`. */
export function expandChecks(options: CheckOptions, extras: CheckExtras = {}): ExpandedCheck[] {
  if (options.optional != null)
    options = {...options, nullable: options.optional};
  return conditions(options).map((condition) => {
    const {key, needsTable, needsCall, when, message, verdicts, verdictMessage, sources} = condition;
    const family = condition.family ?? 'validator';
    let effects: RuleEffect[] = [];
    if (family !== 'validator')
      effects = condition.effects!;
    else if (verdicts != null && extras.message == null && extras.severity == null)
      effects.push({effect: 'verdicts', targets: [TARGET], source: verdicts});
    else {
      const text = extras.message ?? verdictMessage ?? message ?? {map: [{var: verdicts!}, {var: 'message'}]};
      effects.push({effect: extras.severity ?? 'error', targets: [TARGET], message: text});
    }
    return {
      key,
      family,
      needsTable,
      needsCall,
      needsInputs: condition.needsInputs ?? false,
      params: {
        when: extras.when == null ? when : {and: [extras.when, when]},
        effects,
        ...(sources ? {sources} : {}),
      },
    };
  });
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
  if (options.validator)
    checks.validator = options.validator;
  if (options.visible)
    checks.visible = options.visible;
  // the platform keeps `validators` as an array; other options arrive as strings
  const validators = Array.isArray(options.validators) ? options.validators : parseChoices(options.validators);
  if (validators?.length && validators.every((name: unknown) => typeof name === 'string'))
    checks.validators = validators;
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
