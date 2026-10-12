import * as DG from 'datagrok-api/dg';
import {AnnotationLinkKind, LinkSpecString} from '../data/common-types';
import {
  PipelineCheckConfiguration, PipelineHandlerConfiguration, PipelineLinkConfiguration, PipelineLinkConfigurationInput,
  PipelineMetaConfiguration, PipelineRuleConfiguration, PipelineValidatorConfiguration, RuleEffect, RuleExpr, RuleSource,
} from './PipelineConfiguration';
import {CALL, CheckOptions, expandChecks, TABLE, TARGET, VALUE, validateCheckOptions} from './checks';
import {parseLinkIO} from './LinkSpec';
import {compileCheckFormulas, compileRuleFormulas} from './rule-formula';
import {FuncCallIODescription, IOType, expandDeferredIOs, normalizeLinkSpec} from './config-processing-utils';
import {ruleDataHandler, ruleMetaHandler, ruleValidatorHandler} from '../runtime/rule-handlers';
import {ruleTargets, usedAliases} from '../runtime/rule-expressions';
import {DriverLogger, reportError} from '../data/Logger';

const metaEffects = new Set(['hide', 'show', 'items', 'meta']);
const validatorEffects = new Set(['error', 'warning', 'notification', 'verdicts']);
const dataEffects = new Set(['set', 'clear', 'assign']);

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

export function expandLinks(
  links: PipelineLinkConfigurationInput<LinkSpecString>[],
): PipelineLinkConfiguration<LinkSpecString>[] {
  return links.flatMap((link) => isRuleLink(link) ? expandRule(compileRuleFormulas(link)) :
    isCheckLink(link) ? expandCheck(compileCheckFormulas(link)) : [link]);
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

function aliasesOf(ruleId: string, ios: LinkSpecString | undefined, ioType: IOType) {
  const aliases = new Map<string, string>();
  for (const raw of normalizeLinkSpec(ios)) {
    for (const parsed of expandDeferredIOs(parseLinkIO(raw, ioType), ruleId)) {
      const badFlag = parsed.flags?.find((flag) => flag !== 'optional' && flag !== 'template');
      if (badFlag)
        throw new Error(`Rule ${ruleId}: (${badFlag}) flag is not allowed in rule queries (${raw})`);
      if (parsed.name.startsWith('$'))
        throw new Error(`Rule ${ruleId}: alias ${parsed.name} is reserved for the driver`);
      aliases.set(parsed.name, raw);
    }
  }
  return aliases;
}

function effectExpressions(effect: RuleEffect): RuleExpr[] {
  switch (effect.effect) {
  case 'items': return [effect.items];
  case 'meta': return Object.values(effect.meta);
  case 'error':
  case 'warning':
  case 'notification': return [effect.message];
  case 'verdicts': return [];
  case 'set': return [effect.value];
  case 'assign': return [effect.values];
  default: return [];
  }
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
  const toAliases = aliasesOf(id, rule.to, 'output');
  const from = [...normalizeLinkSpec(rule.from)];
  const isAlias = (name: string) => fromAliases.has(name) || Object.hasOwn(sources ?? {}, name);
  const checkExpr = (expr: RuleExpr | undefined) => {
    const fields = new Set<string>();
    for (const alias of usedAliases(expr, fields)) {
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
  const checkArgs = (alias: string, args: unknown) => {
    if (args == null)
      return;
    if (typeof args !== 'object' || Array.isArray(args))
      throw new Error(`Rule ${id}: source ${alias} args must map parameters to expressions`);
    Object.values(args).forEach(checkExpr);
  };
  const expandedSources: Record<string, RuleSource> = {};
  for (const [alias, source] of Object.entries(sources ?? {})) {
    if (alias.startsWith('$'))
      throw new Error(`Rule ${id}: source alias ${alias} is reserved for the driver`);
    if (fromAliases.has(alias))
      throw new Error(`Rule ${id}: source alias ${alias} collides with an input alias`);
    if ('js' in source) {
      const {args, fn} = source.js;
      if (!Array.isArray(args) || args.some((arg) => !fromAliases.has(arg)))
        throw new Error(`Rule ${id}: source ${alias} args must be input aliases`);
      if (typeof fn !== 'function')
        throw new Error(`Rule ${id}: source ${alias} fn must be a function`);
      expandedSources[alias] = source;
      continue;
    }
    if ('func' in source) {
      const {name, args} = source.func;
      if (typeof name !== 'string' || !name)
        throw new Error(`Rule ${id}: source ${alias} name must be a function name`);
      checkArgs(alias, args);
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
      checkArgs(alias, args);
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
    expandedSources[alias] = 'choices' in source ? {choices: {input, call: CALL}} : {validators: {input, call: CALL}};
  }

  checkExpr(when);

  const effectTargets = (effect: RuleEffect) => {
    if (effect.targets == null) {
      if (effect.effect !== 'assign')
        throw new Error(`Rule ${id}: effect ${effect.effect} needs targets`);
      return [...toAliases.keys()];
    }
    return ruleTargets(effect.targets);
  };

  const targeted = new Set<string>();
  const visibilityTargets = new Set<string>();
  for (const effect of effects) {
    if (!metaEffects.has(effect.effect) && !validatorEffects.has(effect.effect) && !dataEffects.has(effect.effect))
      throw new Error(`Rule ${id}: unknown effect ${(effect as any).effect}`);
    if (effect.effect === 'verdicts' && !('validators' in (sources?.[effect.source] ?? {})))
      throw new Error(`Rule ${id}: verdicts effect references unknown validators source ${effect.source}`);
    for (const target of effectTargets(effect)) {
      if (!toAliases.has(target))
        throw new Error(`Rule ${id}: effect ${effect.effect} targets unknown output alias ${target}`);
      targeted.add(target);
      if (effect.effect === 'hide' || effect.effect === 'show') {
        if (visibilityTargets.has(target))
          throw new Error(`Rule ${id}: hide/show applied twice to output alias ${target}`);
        visibilityTargets.add(target);
      }
    }
    effectExpressions(effect).forEach(checkExpr);
    checkExpr(effect.when);
  }
  for (const alias of toAliases.keys()) {
    if (!targeted.has(alias))
      throw new Error(`Rule ${id}: output alias ${alias} is not targeted by any effect`);
  }

  const common = {
    from,
    not: rule.not,
    base: rule.base,
    nodePriority: rule.nodePriority,
    dataFrameMutations: rule.dataFrameMutations,
  };

  const familyLink = (family: Set<string>) => {
    const familyEffects = effects.filter((effect) => family.has(effect.effect));
    if (!familyEffects.length)
      return undefined;
    const familyTargets = new Set(familyEffects.flatMap(effectTargets));
    const to = [...new Set([...toAliases].filter(([alias]) => familyTargets.has(alias)).map(([, raw]) => raw))];
    return {to, params: {when, effects: familyEffects, ...(sources ? {sources: expandedSources} : {})}};
  };

  const result: PipelineLinkConfiguration<LinkSpecString>[] = [];
  const meta = familyLink(metaEffects);
  if (meta) {
    const link: PipelineMetaConfiguration<LinkSpecString> = {
      ...common, ...meta, id: `${id}::meta`, type: 'meta', handler: ruleMetaHandler,
    };
    result.push(link);
  }
  const validator = familyLink(validatorEffects);
  if (validator) {
    const link: PipelineValidatorConfiguration<LinkSpecString> = {
      ...common, ...validator, id: `${id}::validator`, type: 'validator', handler: ruleValidatorHandler,
      debounce: rule.debounce,
    };
    result.push(link);
  }
  const data = familyLink(dataEffects);
  if (data) {
    const link: PipelineHandlerConfiguration<LinkSpecString> = {
      ...common, ...data, id: `${id}::data`, type: 'data', handler: ruleDataHandler, runOnInit: rule.runOnInit,
    };
    result.push(link);
  }
  return result;
}
