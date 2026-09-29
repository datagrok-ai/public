import {LinkSpecString} from '../data/common-types';
import {
  PipelineCheckConfiguration, PipelineHandlerConfiguration, PipelineLinkConfiguration, PipelineLinkConfigurationInput,
  PipelineMetaConfiguration, PipelineRuleConfiguration, PipelineValidatorConfiguration, RuleEffect, RuleExpr, RuleSource,
} from './PipelineConfiguration';
import {CALL, CheckOptions, expandChecks, TABLE, TARGET, VALUE, validateCheckOptions} from './checks';
import {parseLinkIO} from './LinkSpec';
import {IOType, normalizeLinkSpec} from './config-processing-utils';
import {ruleDataHandler, ruleMetaHandler, ruleValidatorHandler} from '../runtime/rule-handlers';
import {ruleTargets, usedAliases} from '../runtime/rule-expressions';

const metaEffects = new Set(['hide', 'show', 'items', 'meta']);
const validatorEffects = new Set(['error', 'warning', 'notification', 'verdicts']);
const dataEffects = new Set(['set', 'clear']);

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

export function expandLinks(
  links: PipelineLinkConfigurationInput<LinkSpecString>[],
): PipelineLinkConfiguration<LinkSpecString>[] {
  return links.flatMap((link) => isRuleLink(link) ? expandRule(link) : isCheckLink(link) ? expandCheck(link) : [link]);
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
    if ([VALUE, TABLE, TARGET, CALL].includes(alias))
      throw new Error(`Check ${id}: vars alias ${alias} is reserved`);
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
    for (const parsed of parseLinkIO(raw, ioType)) {
      const badFlag = parsed.flags?.find((flag) => flag !== 'optional');
      if (badFlag)
        throw new Error(`Rule ${ruleId}: (${badFlag}) flag is not allowed in rule queries (${raw})`);
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
  default: return [];
  }
}

function expandRule(rule: PipelineRuleConfiguration<LinkSpecString>): PipelineLinkConfiguration<LinkSpecString>[] {
  const {id, when, effects, sources} = rule;
  if (rule.handler)
    throw new Error(`Rule ${id}: handler is not allowed, rules use built-in handlers`);
  if (!effects?.length)
    throw new Error(`Rule ${id}: effects list is empty`);

  const fromAliases = aliasesOf(id, rule.from, 'input');
  const toAliases = aliasesOf(id, rule.to, 'output');
  const from = [...normalizeLinkSpec(rule.from)];
  const expandedSources: Record<string, RuleSource> = {};
  for (const [alias, source] of Object.entries(sources ?? {})) {
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
    if (!('validators' in source))
      throw new Error(`Rule ${id}: unknown source kind for alias ${alias}`);
    const {input, names} = source.validators;
    if (!fromAliases.has(input))
      throw new Error(`Rule ${id}: source ${alias} references unknown input alias ${input}`);
    if (names != null && (!Array.isArray(names) || names.some((name) => typeof name !== 'string')))
      throw new Error(`Rule ${id}: source ${alias} names must be an array of function names`);
    if (names) {
      expandedSources[alias] = {validators: {input, names}};
      continue;
    }
    // the annotation's validators need the step's FuncCall: derive it from the input's query
    if (fromAliases.has(CALL))
      throw new Error(`Rule ${id}: input alias ${CALL} is reserved for sources`);
    const raw = fromAliases.get(input)!;
    const tail = raw.slice(raw.indexOf(':') + 1);
    const cut = tail.lastIndexOf('/');
    if (cut < 0)
      throw new Error(`Rule ${id}: source ${alias} needs an io query for ${input}`);
    const callQuery = `${CALL}(call,optional):${tail.slice(0, cut)}`;
    parseLinkIO(callQuery, 'input');
    if (!from.includes(callQuery))
      from.push(callQuery);
    expandedSources[alias] = {validators: {input, call: CALL}};
  }

  const checkExpr = (expr: RuleExpr) => {
    for (const alias of usedAliases(expr)) {
      if (!fromAliases.has(alias) && !(alias in (sources ?? {})))
        throw new Error(`Rule ${id}: expression references unknown input alias ${alias}`);
    }
  };
  checkExpr(when);

  const targeted = new Set<string>();
  const visibilityTargets = new Set<string>();
  for (const effect of effects) {
    if (!metaEffects.has(effect.effect) && !validatorEffects.has(effect.effect) && !dataEffects.has(effect.effect))
      throw new Error(`Rule ${id}: unknown effect ${(effect as any).effect}`);
    if (effect.effect === 'verdicts' && !('validators' in (sources?.[effect.source] ?? {})))
      throw new Error(`Rule ${id}: verdicts effect references unknown validators source ${effect.source}`);
    for (const target of ruleTargets(effect.targets)) {
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
    const familyTargets = new Set(familyEffects.flatMap((effect) => ruleTargets(effect.targets)));
    const to = [...toAliases].filter(([alias]) => familyTargets.has(alias)).map(([, raw]) => raw);
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
