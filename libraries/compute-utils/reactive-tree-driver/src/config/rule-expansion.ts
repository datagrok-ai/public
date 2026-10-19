import * as DG from 'datagrok-api/dg';
import {AnnotationLinkKind, LinkSpecString} from '../data/common-types';
import {
  PipelineCheckConfiguration, PipelineLinkConfiguration, PipelineLinkConfigurationBase,
  PipelineLinkConfigurationInput, PipelineMetaConfiguration, PipelineRuleConfiguration, PipelineValidatorConfiguration,
  RuleEffect, RuleExpr, RuleMessageAction, RuleSource,
} from './PipelineConfiguration';
import {CALL, CheckOptions, expandChecks, TABLE, TARGET, VALUE, validateCheckOptions} from './checks';
import {LinkIOParsed, parseLinkIO} from './LinkSpec';
import {compileCheckFormulas, compileRuleFormulas} from './rule-formula';
import {FuncCallIODescription, IOType, expandDeferredIOs, normalizeLinkSpec} from './config-processing-utils';
import {
  ruleDataHandler, ruleMetaHandler, ruleNodeMetaHandler, ruleValidatorHandler,
} from '../runtime/rule-handlers';
import {descriptionOutputs} from '../runtime/StateTreeNodes';
import {indexFromEnd} from '../utils';
import {ruleTargets, usedAliases} from '../runtime/rule-expressions';
import {sourceReads} from '../runtime/rule-sources';
import {DriverLogger, reportError} from '../data/Logger';

const metaEffects = new Set(['hide', 'show', 'items', 'meta']);
const validatorEffects = new Set(['error', 'warning', 'notification', 'verdicts']);
const dataEffects = new Set(['set', 'clear', 'assign']);
const messageEffects = new Set(['error', 'warning', 'notification']);
const actionEffects = new Set(['hideAction', 'showAction']);
const driverValues = new Set(['$position', '$steps']);

type Family = 'meta' | 'validator' | 'data' | 'nodemeta' | 'pipelineValidator';

export function isRuleLink(
  link: PipelineLinkConfigurationInput<LinkSpecString>,
): link is PipelineRuleConfiguration<LinkSpecString> {
  return link.type === 'rule';
}

export function isCheckLink(
  link: PipelineLinkConfigurationInput<LinkSpecString>,
): link is PipelineCheckConfiguration<LinkSpecString> {
  return link.type === 'check';
}

/** The rules a step's annotations stand for: for each input whose choices the platform evaluates,
 *  the item list, and for the first `propagateChoice: all` input the lookup writing the picked row
 *  into the other scalar inputs, with a warning for cells that do not fit their inputs. */
export function annotationRules(
  nqName: string, io: FuncCallIODescription[], logger?: DriverLogger,
): PipelineRuleConfiguration<LinkSpecString>[] {
  // without the platform evaluation (before 1.28) the rules would only ever be off
  if (typeof (DG.FuncCall.prototype as any).evalParamChoices !== 'function')
    return [];
  const inputs = io.filter((item) => item.direction === 'input');
  // aliases next to the io names the template queries produce
  const free = (name: string): string => inputs.some((other) => other.id === name) ? free(`${name}_`) : name;
  // the first lookup already writes every other scalar input, so a second one would have no targets of its own
  const lookupKey = inputs.find((item) => item.dynamicChoices?.propagate)?.id;
  for (const item of inputs) {
    if (item.dynamicChoices?.propagate && item.id !== lookupKey) {
      reportError('warning', 'configProcessing',
        `Step ${nqName}: propagateChoice on '${item.id}' is ignored, '${lookupKey}' already fills the step's inputs`,
        logger);
    }
  }
  return inputs.filter((item) => item.dynamicChoices).flatMap(({id}) => {
    const lookup = id === lookupKey;
    const choices = free(`${id}_choices`);
    const target = free(`${id}_target`);
    const source = {[choices]: {choices: {input: id}}};
    const rules: PipelineRuleConfiguration<LinkSpecString>[] = [{
      id: `::${id}:choices`,
      type: 'rule',
      debounce: 0,
      from: `_(template):inputs(${nqName})`,
      to: `${target}:${id}`,
      sources: source,
      effects: [
        {effect: 'items', targets: target, items: {var: `${choices}.items`}, when: {'!!': {var: choices}}},
        ...(lookup ? [{
          effect: 'warning' as const, targets: target, message: {var: `${choices}.rowErrors`},
          when: {'!!': {var: choices}},
        }] : []),
      ],
    }];
    if (lookup) {
      rules.push({
        id: `::${id}:lookup`,
        type: 'rule',
        runOnInit: true,
        from: `${id}:${id}`,
        to: `_(template):inputs(${nqName}, ${id}|$nonscalar|$linked)`,
        sources: source,
        effects: [{
          effect: 'assign', values: {var: `${choices}.row`}, ignoreCase: true, restriction: 'restricted',
          when: {'!!': {var: `${choices}.row`}},
        }],
      });
    }
    return rules;
  });
}

/** The checks a step's annotations stand for, one link per input and check; a GrokScript expression
 *  sees every input of the step under its own name. */
export function annotationCheckLinks(
  io: FuncCallIODescription[],
): {kind: AnnotationLinkKind, link: PipelineLinkConfiguration<LinkSpecString>}[] {
  const inputs = io.filter((item) => item.direction === 'input');
  const inputQueries = inputs.filter((other) => other.id !== VALUE).map((other) => `${other.id}:${other.id}`);
  return inputs.flatMap((item) =>
    expandChecks({...item.checks, nullable: item.nullable}, {validatorsViaCall: true})
      .map(({key, family, needsCall, needsInputs, params}) => {
        const kind: AnnotationLinkKind = key === 'required' ? 'required' : 'check';
        const common = {
          id: `::${item.id}:${key}`,
          from: [`${VALUE}:${item.id}`, ...(needsCall ? [`${CALL}(call,optional):.`] : []),
            ...(needsInputs ? inputQueries : [])],
          to: `${TARGET}:${item.id}`,
          params,
        };
        const link: PipelineMetaConfiguration<LinkSpecString> | PipelineValidatorConfiguration<LinkSpecString> =
          family === 'meta' ? {...common, type: 'meta', handler: ruleMetaHandler} :
            {...common, type: 'validator', handler: ruleValidatorHandler, debounce: 0};
        return {kind, link};
      }));
}

export function expandLink(
  link: PipelineLinkConfigurationInput<LinkSpecString>,
): PipelineLinkConfiguration<LinkSpecString>[] {
  return isRuleLink(link) ? expandRule(compileRuleFormulas(link)) :
    isCheckLink(link) ? expandCheck(compileCheckFormulas(link)) : [link];
}

function singleQuery(id: string, field: string, query: LinkSpecString | undefined): string | undefined {
  if (query == null)
    return undefined;
  if (typeof query !== 'string' && !Array.isArray(query))
    throw new Error(`Check ${id}: ${field} must be a query`);
  if (Array.isArray(query)) {
    if (query.length !== 1)
      throw new Error(`Check ${id}: ${field} must be a single query`);
    return query[0];
  }
  return query;
}

function expandCheck(check: PipelineCheckConfiguration<LinkSpecString>): PipelineLinkConfiguration<LinkSpecString>[] {
  const {id} = check;
  const options: CheckOptions = {...check.check};
  validateCheckOptions(id, options);
  if (options.nullable === true || options.optional === true)
    throw new Error(`Check ${id}: nullable: true only relaxes the default check as a script annotation`);
  const io = singleQuery(id, 'io', check.io)!;
  const table = singleQuery(id, 'table', options.table);
  for (const alias of usedAliases(check.when)) {
    if (alias !== VALUE && alias !== TABLE)
      throw new Error(`Check ${id}: when references unknown alias ${alias}, use ${VALUE} or ${TABLE}`);
  }
  const expanded = expandChecks(options, {when: check.when, message: check.message, severity: check.severity});
  if (!expanded.length)
    throw new Error(`Check ${id}: no options to check`);
  const vars = Object.entries(check.vars ?? {}).map(([alias, query]) => {
    if (alias === VALUE)
      throw new Error(`Check ${id}: vars alias ${VALUE} is the checked io`);
    if (alias.startsWith('$'))
      throw new Error(`Check ${id}: vars alias ${alias} is reserved for the driver`);
    return `${alias}:${singleQuery(id, `vars.${alias}`, query)}`;
  });
  return expanded.map(({key, family, needsTable, needsInputs, params}) => {
    const from = [`${VALUE}:${io}`];
    if (needsTable)
      from.push(`${TABLE}:${table}`);
    if (needsInputs)
      from.push(...vars);
    const common = {
      id: `${id}::${key}`,
      from,
      to: [`${TARGET}:${io}`],
      not: check.not,
      base: check.base,
      nodePriority: check.nodePriority,
      params,
    };
    if (family === 'meta') {
      const link: PipelineMetaConfiguration<LinkSpecString> = {...common, type: 'meta', handler: ruleMetaHandler};
      return link;
    }
    const link: PipelineValidatorConfiguration<LinkSpecString> = {
      ...common, type: 'validator', handler: ruleValidatorHandler, debounce: check.debounce ?? 0,
    };
    return link;
  });
}

// `nodes` collects the (node) aliases; parseLinkIO allows them only in `to`
function aliasesOf(ruleId: string, ios: LinkSpecString | undefined, ioType: IOType, nodes?: Set<string>) {
  const aliases = new Map<string, string>();
  // a node meta link does not expand its actions queries
  const allowed = ioType === 'actions' ? ['optional'] : ['optional', 'template', 'node'];
  for (const raw of normalizeLinkSpec(ios)) {
    for (const parsed of expandDeferredIOs(parseLinkIO(raw, ioType), ruleId)) {
      const badFlag = parsed.flags?.find((flag) => !allowed.includes(flag));
      if (badFlag)
        throw new Error(`Rule ${ruleId}: (${badFlag}) flag is not allowed in rule queries (${raw})`);
      if (parsed.flags?.includes('node')) {
        if (parsed.flags.includes('template'))
          throw new Error(`Rule ${ruleId}: (node) and (template) flags cannot be combined (${raw})`);
        nodes!.add(parsed.name);
      }
      if (parsed.name.startsWith('$'))
        throw new Error(`Rule ${ruleId}: alias ${parsed.name} is reserved for the driver`);
      aliases.set(parsed.name, raw);
    }
  }
  return aliases;
}

function effectExpressions(effect: RuleEffect): RuleExpr[] {
  switch (effect.effect) {
  case 'describe': return [effect.value];
  case 'items': return [effect.items];
  case 'meta': return Object.values(effect.meta);
  case 'error':
  case 'warning':
  case 'notification': return [effect.message, ...(effect.actions ?? []).flatMap((item) =>
    item.additionalParams === undefined ? [item.actionName] : [item.actionName, item.additionalParams])];
  case 'verdicts': return [];
  case 'set': return [effect.value];
  case 'assign': return [effect.values];
  default: return [];
  }
}

const offeredActions = (effect: RuleEffect): RuleMessageAction[] =>
  (effect as {actions?: RuleMessageAction[]}).actions ?? [];

function ruleLink(
  family: Family, link: Omit<PipelineLinkConfigurationBase<LinkSpecString>, 'actions'>,
  rule: PipelineRuleConfiguration<LinkSpecString>, actions: string[],
): PipelineLinkConfiguration<LinkSpecString> {
  switch (family) {
  case 'meta': return {...link, type: 'meta', handler: ruleMetaHandler};
  case 'validator':
    return {...link, type: 'validator', handler: ruleValidatorHandler, debounce: rule.debounce,
      ...(actions.length ? {actions} : {})};
  case 'data': return {...link, type: 'data', handler: ruleDataHandler, runOnInit: rule.runOnInit};
  case 'nodemeta': return {...link, type: 'nodemeta', handler: ruleNodeMetaHandler, actions};
  case 'pipelineValidator':
    return {...link, type: 'pipelineValidator', handler: ruleValidatorHandler, debounce: rule.debounce};
  }
}

// a template query yields every io it matches, so a rule link keeps only the aliases its own effects target,
// and a separate effect's link only the inputs it reads
const linkAliasesOf = new WeakMap<object, {to: Set<string>, from?: Set<string>}>();

export function narrowRuleLink<L extends {from: LinkIOParsed[], to: LinkIOParsed[]}>(
  expanded: object, processed: L,
): L {
  const aliases = linkAliasesOf.get(expanded);
  if (aliases) {
    processed.to = processed.to.filter((io) => aliases.to.has(io.name));
    if (aliases.from)
      processed.from = processed.from.filter((io) => aliases.from!.has(io.name));
  }
  return processed;
}

// expects compiled formulas: every effect and source is an object
function expandRule(rule: PipelineRuleConfiguration<LinkSpecString>): PipelineLinkConfiguration<LinkSpecString>[] {
  const {id, when} = rule;
  const effects = rule.effects as RuleEffect[];
  const sources = rule.sources as Record<string, RuleSource> | undefined;
  if (rule.handler)
    throw new Error(`Rule ${id}: handler is not allowed, rules use built-in handlers`);
  if (!effects?.length)
    throw new Error(`Rule ${id}: effects list is empty`);

  const fromAliases = aliasesOf(id, rule.from, 'input');
  const nodeAliases = new Set<string>();
  const toAliases = aliasesOf(id, rule.to, 'output', nodeAliases);
  const actionAliases = aliasesOf(id, rule.actions, 'actions');
  for (const alias of actionAliases.keys()) {
    if (fromAliases.has(alias) || toAliases.has(alias))
      throw new Error(`Rule ${id}: actions alias ${alias} collides with an io alias`);
  }
  const from = [...normalizeLinkSpec(rule.from)];
  const isRuleAlias = (name: string) => fromAliases.has(name) || Object.hasOwn(sources ?? {}, name);
  const checkExpr = (expr: RuleExpr | undefined, isAlias = isRuleAlias) => {
    const fields = new Set<string>();
    for (const alias of usedAliases(expr, fields)) {
      if (driverValues.has(alias)) {
        if (alias === '$position' && !rule.base)
          throw new Error(`Rule ${id}: $position needs a base`);
        continue;
      }
      if (!isAlias(alias))
        throw new Error(`Rule ${id}: expression references unknown input alias ${alias}`);
    }
    // reduce's own names are not aliases there
    for (const name of fields) {
      if (isAlias(name) && name !== 'current' && name !== 'accumulator') {
        throw new Error(`Rule ${id}: ${name} inside map, filter, all, some, none or reduce is a field of the ` +
          `element, not the alias ${name}; read the field with var("${name}")`);
      }
    }
  };
  // a source reads inputs and the sources declared before it, which resolve first
  const sourceAliases = Object.keys(sources ?? {});
  const declared = new Set<string>();
  const isReadable = (name: string) => fromAliases.has(name) || declared.has(name);
  const checkOrder = (alias: string, source: RuleSource) => {
    for (const name of sourceReads(source)) {
      if (name === alias)
        throw new Error(`Rule ${id}: source ${alias} reads itself`);
      if (sourceAliases.indexOf(name) > sourceAliases.indexOf(alias))
        throw new Error(`Rule ${id}: source ${alias} reads ${name}, which is declared after it`);
    }
  };
  const checkArgs = (alias: string, source: RuleSource, args: unknown) => {
    if (args == null)
      return;
    if (typeof args !== 'object' || Array.isArray(args))
      throw new Error(`Rule ${id}: source ${alias} args must map parameters to expressions`);
    checkOrder(alias, source);
    Object.values(args).forEach((arg) => checkExpr(arg, isReadable));
  };
  const expandedSources: Record<string, RuleSource> = {};
  const callQueries = new Map<string, string>();
  for (const [alias, source] of Object.entries(sources ?? {})) {
    declared.add(alias);
    if (alias.startsWith('$'))
      throw new Error(`Rule ${id}: source alias ${alias} is reserved for the driver`);
    // JS orders integer keys first, which would break the declaration order
    if (/^(0|[1-9]\d*)$/.test(alias))
      throw new Error(`Rule ${id}: source alias ${alias} must not be a number`);
    if (fromAliases.has(alias))
      throw new Error(`Rule ${id}: source alias ${alias} collides with an input alias`);
    if ('js' in source) {
      const {args, fn} = source.js;
      if (Array.isArray(args))
        checkOrder(alias, source);
      if (!Array.isArray(args) || args.some((arg) => !isReadable(arg)))
        throw new Error(`Rule ${id}: source ${alias} args must be input or earlier source aliases`);
      if (typeof fn !== 'function')
        throw new Error(`Rule ${id}: source ${alias} fn must be a function`);
      expandedSources[alias] = source;
      continue;
    }
    if ('func' in source) {
      const {name, args} = source.func;
      if (typeof name !== 'string' || !name)
        throw new Error(`Rule ${id}: source ${alias} name must be a function name`);
      checkArgs(alias, source, args);
      expandedSources[alias] = source;
      continue;
    }
    if ('file' in source) {
      if (typeof source.file !== 'string' || !source.file)
        throw new Error(`Rule ${id}: source ${alias} file must be a path`);
      expandedSources[alias] = source;
      continue;
    }
    if ('table' in source) {
      const {table} = source;
      const csv = typeof table === 'string' ? table : (table as {csv?: unknown} | null)?.csv;
      if (!(table instanceof DG.DataFrame) && (typeof csv !== 'string' || !csv))
        throw new Error(`Rule ${id}: source ${alias} table must be a dataframe or CSV text`);
      expandedSources[alias] = source;
      continue;
    }
    if ('query' in source) {
      const {connection, sql, args} = source.query;
      if (typeof connection !== 'string' || !connection || typeof sql !== 'string' || !sql)
        throw new Error(`Rule ${id}: source ${alias} needs a connection and sql`);
      checkArgs(alias, source, args);
      expandedSources[alias] = source;
      continue;
    }
    if (!('validators' in source) && !('choices' in source))
      throw new Error(`Rule ${id}: unknown source kind for alias ${alias}`);
    const {input} = 'choices' in source ? source.choices : source.validators;
    if (!fromAliases.has(input))
      throw new Error(`Rule ${id}: source ${alias} references unknown input alias ${input}`);
    if ('validators' in source) {
      const {names} = source.validators;
      if (names != null && (!Array.isArray(names) || names.some((name) => typeof name !== 'string')))
        throw new Error(`Rule ${id}: source ${alias} names must be an array of function names`);
      if (names) {
        expandedSources[alias] = {validators: {input, names}};
        continue;
      }
    }
    // annotation validators and choices need the step's FuncCall: derive it from the input's query
    const raw = fromAliases.get(input)!;
    const tail = raw.slice(raw.indexOf(':') + 1);
    const cut = tail.lastIndexOf('/');
    // a query without a node path addresses an io of the node the rule is defined on
    const callQuery = `${CALL}(call,optional):${cut < 0 ? '.' : tail.slice(0, cut)}`;
    parseLinkIO(callQuery, 'input');
    if (!from.includes(callQuery))
      from.push(callQuery);
    callQueries.set(alias, callQuery);
    expandedSources[alias] = 'choices' in source ? {choices: {input, call: CALL}} : {validators: {input, call: CALL}};
  }

  checkExpr(when);

  // describe targets a description field, so it is not one of the io aliases assign defaults to
  const describeTargets = new Set(effects.flatMap((effect) =>
    effect.effect === 'describe' && effect.targets != null ? ruleTargets(effect.targets) : []));
  const effectTargets = (effect: RuleEffect) => {
    if (effect.targets == null) {
      if (effect.effect !== 'assign')
        throw new Error(`Rule ${id}: effect ${effect.effect} needs targets`);
      return [...toAliases.keys()].filter((alias) => !nodeAliases.has(alias) && !describeTargets.has(alias));
    }
    return ruleTargets(effect.targets);
  };

  const targetFamily = (effect: RuleEffect, target: string): Family => {
    if (actionAliases.has(target)) {
      if (!actionEffects.has(effect.effect))
        throw new Error(`Rule ${id}: actions alias ${target} takes only hideAction or showAction`);
      return 'nodemeta';
    }
    if (actionEffects.has(effect.effect))
      throw new Error(`Rule ${id}: effect ${effect.effect} targets unknown actions alias ${target}`);
    if (!toAliases.has(target))
      throw new Error(`Rule ${id}: effect ${effect.effect} targets unknown output alias ${target}`);
    if (nodeAliases.has(target)) {
      if (!messageEffects.has(effect.effect))
        throw new Error(`Rule ${id}: (node) alias ${target} takes only error, warning or notification`);
      return 'pipelineValidator';
    }
    if (effect.effect === 'describe') {
      for (const parsed of parseLinkIO(toAliases.get(target)!, 'output')) {
        const field = indexFromEnd(parsed.segments);
        if (parsed.flags?.includes('template') || field?.type !== 'selector' ||
          !(descriptionOutputs as readonly string[]).includes(field.ids[0])) {
          throw new Error(`Rule ${id}: describe target ${target} must end with one of ` +
            `${descriptionOutputs.join(', ')}`);
        }
      }
      return 'nodemeta';
    }
    return metaEffects.has(effect.effect) ? 'meta' : validatorEffects.has(effect.effect) ? 'validator' : 'data';
  };

  const families = new Map<RuleEffect, Family>();
  const aliasKinds = new Map<string, string>();
  const targeted = new Set<string>();
  const visibilityTargets = new Set<string>();
  for (const effect of effects) {
    if (!metaEffects.has(effect.effect) && !validatorEffects.has(effect.effect) && !dataEffects.has(effect.effect) &&
      !actionEffects.has(effect.effect) && effect.effect !== 'describe')
      throw new Error(`Rule ${id}: unknown effect ${(effect as any).effect}`);
    if (effect.effect === 'verdicts' && !('validators' in (sources?.[effect.source] ?? {})))
      throw new Error(`Rule ${id}: verdicts effect references unknown validators source ${effect.source}`);
    if (effect.separate != null && typeof effect.separate !== 'boolean')
      throw new Error(`Rule ${id}: separate of effect ${effect.effect} must be true or false`);
    if (effect.additionalTriggers != null) {
      if (effect.separate === false)
        throw new Error(`Rule ${id}: additionalTriggers of effect ${effect.effect} needs a separate link`);
      for (const alias of ruleTargets(effect.additionalTriggers)) {
        if (!fromAliases.has(alias))
          throw new Error(`Rule ${id}: additionalTriggers of ${effect.effect} names unknown input alias ${alias}`);
      }
    }
    const offered = offeredActions(effect);
    if (!Array.isArray(offered) || (offered.length && !messageEffects.has(effect.effect)))
      throw new Error(`Rule ${id}: actions are a list of message buttons, on error, warning or notification only`);
    for (const item of offered) {
      if (!actionAliases.has(item.actions))
        throw new Error(`Rule ${id}: message action references unknown actions alias ${item.actions}`);
      if (typeof item.action !== 'string' || !item.action)
        throw new Error(`Rule ${id}: message action on ${item.actions} needs an action id`);
      targeted.add(item.actions);
    }
    if ((effect.effect === 'hideAction' || effect.effect === 'showAction') &&
      (typeof effect.action !== 'string' || !effect.action))
      throw new Error(`Rule ${id}: ${effect.effect} needs an action id`);
    for (const target of effectTargets(effect)) {
      const family = targetFamily(effect, target);
      if (family === 'pipelineValidator' && offered.length)
        throw new Error(`Rule ${id}: (node) alias ${target} messages take no actions`);
      if (families.has(effect) && families.get(effect) !== family)
        throw new Error(`Rule ${id}: effect ${effect.effect} mixes targets of different kinds`);
      families.set(effect, family);
      const kind = family === 'nodemeta' ? (actionAliases.has(target) ? 'action' : 'field') :
        family === 'pipelineValidator' ? 'node' : 'io';
      if (aliasKinds.has(target) && aliasKinds.get(target) !== kind)
        throw new Error(`Rule ${id}: alias ${target} is targeted as both ${aliasKinds.get(target)} and ${kind}`);
      aliasKinds.set(target, kind);
      targeted.add(target);
      const key = effect.effect === 'hide' || effect.effect === 'show' ? `output alias ${target}` :
        effect.effect === 'hideAction' || effect.effect === 'showAction' ? `action ${target}/${effect.action}` :
          undefined;
      if (key) {
        if (visibilityTargets.has(key))
          throw new Error(`Rule ${id}: hide/show applied twice to ${key}`);
        visibilityTargets.add(key);
      }
    }
    effectExpressions(effect).forEach((expr) => checkExpr(expr));
    checkExpr(effect.when);
  }
  for (const [kind, aliases] of [['output', toAliases], ['actions', actionAliases]] as const) {
    for (const alias of aliases.keys()) {
      if (!targeted.has(alias))
        throw new Error(`Rule ${id}: ${kind} alias ${alias} is not targeted by any effect`);
    }
  }

  // a link re-runs on every rewire when its own expressions read $steps; `when` and sources are shared
  const readsStepsIn = (expr: RuleExpr | undefined) => usedAliases(expr).includes('$steps');
  const sharedSteps = readsStepsIn(when) || Object.values(expandedSources).some((source) =>
    sourceReads(source).includes('$steps'));

  // a separate effect's link reads what the effect, the rule's `when` and the sources these use read
  const separateInputs = (effect: RuleEffect) => {
    const read = new Set([...effectExpressions(effect), effect.when, when].flatMap((expr) => usedAliases(expr)));
    if (effect.effect === 'verdicts')
      read.add(effect.source);
    const aliases = Object.keys(expandedSources);
    // a source reads only earlier ones, so one pass from the last collects the sources read through others
    for (const alias of [...aliases].reverse()) {
      if (read.has(alias))
        sourceReads(expandedSources[alias]).forEach((name) => read.add(name));
    }
    const used = aliases.filter((alias) => read.has(alias));
    ruleTargets(effect.additionalTriggers ?? []).forEach((alias) => read.add(alias));
    return {
      inputs: [...fromAliases.keys()].filter((alias) => read.has(alias)),
      calls: [...new Set(used.flatMap((alias) => callQueries.get(alias) ?? []))],
      sources: used,
      readsSteps: read.has('$steps'),
    };
  };

  const raws = (aliases: Map<string, string>, names: Set<string>) =>
    [...new Set([...aliases].filter(([alias]) => names.has(alias)).map(([, raw]) => raw))];

  const makeLink = (
    family: Family, linkEffects: RuleEffect[], suffix: string, separate?: ReturnType<typeof separateInputs>,
  ): PipelineLinkConfiguration<LinkSpecString> => {
    const targets = new Set(linkEffects.flatMap(effectTargets));
    const readsSteps = separate ? separate.readsSteps :
      sharedSteps || linkEffects.some((effect) => [...effectExpressions(effect), effect.when].some(readsStepsIn));
    const linkSources = separate ?
      (separate.sources.length ? Object.fromEntries(separate.sources.map((alias) => [alias, expandedSources[alias]])) :
        undefined) :
      sources && expandedSources;
    const inputs = separate && new Set([...separate.inputs, ...(separate.calls.length ? [CALL] : [])]);
    const link = {
      id: `${id}::${suffix}`,
      from: separate ? [...raws(fromAliases, inputs!), ...separate.calls] : from,
      to: raws(toAliases, targets),
      not: rule.not,
      base: rule.base,
      nodePriority: rule.nodePriority,
      dataFrameMutations: Array.isArray(rule.dataFrameMutations) && inputs ?
        rule.dataFrameMutations.filter((alias) => inputs.has(alias)) : rule.dataFrameMutations,
      params: {
        when, effects: linkEffects, ...(linkSources ? {sources: linkSources} : {}), ...(readsSteps ? {readsSteps} : {}),
      },
    };
    const actionNames = new Set([...targets, ...linkEffects.flatMap(offeredActions).map((item) => item.actions)]);
    const typed = ruleLink(family, link, rule, raws(actionAliases, actionNames));
    linkAliasesOf.set(typed, {to: targets, from: inputs});
    return typed;
  };

  const isSeparate = (effect: RuleEffect) => !!effect.separate || effect.additionalTriggers != null;
  const result: PipelineLinkConfiguration<LinkSpecString>[] = [];
  for (const family of ['meta', 'validator', 'data', 'nodemeta', 'pipelineValidator'] as const) {
    const familyEffects = effects.filter((effect) => families.get(effect) === family && !isSeparate(effect));
    if (familyEffects.length)
      result.push(makeLink(family, familyEffects, family));
  }
  effects.forEach((effect, idx) => {
    if (!isSeparate(effect))
      return;
    const family = families.get(effect)!;
    const separate = separateInputs(effect);
    if (!separate.inputs.length && !separate.calls.length && !separate.readsSteps &&
      !(family === 'data' && rule.runOnInit)) {
      throw new Error(`Rule ${id}: separate effect ${effect.effect} reads no input, ` +
        `add additionalTriggers${family === 'data' ? ' or runOnInit' : ''}`);
    }
    result.push(makeLink(family, [effect], `${family}:${idx + 1}`, separate));
  });
  return result;
}
