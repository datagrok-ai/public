import {LinkSpecString} from '../data/common-types';
import {
  PipelineCheckConfiguration, PipelineHandlerConfiguration, PipelineLinkConfiguration, PipelineLinkConfigurationInput,
  PipelineMetaConfiguration, PipelineRuleConfiguration, PipelineValidatorConfiguration, RuleEffect, RuleExpr,
} from './PipelineConfiguration';
import {CheckOptions, expandChecks, TABLE, TARGET, VALUE, validateCheckOptions} from './checks';
import {parseLinkIO} from './LinkSpec';
import {IOType, normalizeLinkSpec} from './config-processing-utils';
import {ruleDataHandler, ruleMetaHandler, ruleValidatorHandler} from '../runtime/rule-handlers';
import {ruleTargets, usedAliases} from '../runtime/rule-expressions';

const metaEffects = new Set(['hide', 'show', 'items', 'meta']);
const validatorEffects = new Set(['error', 'warning', 'notification']);
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
  if (options.nullable === true)
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
  return expanded.map(({key, needsTable, params}) => {
    const link: PipelineValidatorConfiguration<LinkSpecString> = {
      id: `${id}::${key}`,
      type: 'validator',
      from: needsTable ? [`${VALUE}:${io}`, `${TABLE}:${table}`] : [`${VALUE}:${io}`],
      to: [`${TARGET}:${io}`],
      not: check.not,
      base: check.base,
      nodePriority: check.nodePriority,
      debounce: check.debounce ?? 0,
      handler: ruleValidatorHandler,
      params,
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
  case 'set': return [effect.value];
  default: return [];
  }
}

function expandRule(rule: PipelineRuleConfiguration<LinkSpecString>): PipelineLinkConfiguration<LinkSpecString>[] {
  const {id, when, effects} = rule;
  if (rule.handler)
    throw new Error(`Rule ${id}: handler is not allowed, rules use built-in handlers`);
  if (!effects?.length)
    throw new Error(`Rule ${id}: effects list is empty`);

  const fromAliases = aliasesOf(id, rule.from, 'input');
  const toAliases = aliasesOf(id, rule.to, 'output');

  const checkExpr = (expr: RuleExpr) => {
    for (const alias of usedAliases(expr)) {
      if (!fromAliases.has(alias))
        throw new Error(`Rule ${id}: expression references unknown input alias ${alias}`);
    }
  };
  checkExpr(when);

  const targeted = new Set<string>();
  const visibilityTargets = new Set<string>();
  for (const effect of effects) {
    if (!metaEffects.has(effect.effect) && !validatorEffects.has(effect.effect) && !dataEffects.has(effect.effect))
      throw new Error(`Rule ${id}: unknown effect ${(effect as any).effect}`);
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
    from: rule.from,
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
    return {to, params: {when, effects: familyEffects}};
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
