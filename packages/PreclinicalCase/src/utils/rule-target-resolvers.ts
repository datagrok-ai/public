import {IssueDetail} from '../types/validation-result';

export type TargetResolverResult = {
  targets: Set<string>;
  context: Set<string>;
};

type RuleTargetResolver = (issue: IssueDetail) => TargetResolverResult;

const ISO_8601_REGEX = /^\d{4}-\d{2}-\d{2}(T\d{2}:\d{2}:\d{2}(\.\d{1,3})?(Z|[+-]\d{2}:\d{2})?)?$/;

function ensureTargetsNotEmpty(result: TargetResolverResult): TargetResolverResult {
  if (result.targets.size === 0) {
    for (const v of result.context)
      result.targets.add(v);
    result.context.clear();
  }
  return result;
}

function resolveISO8601Targets(issue: IssueDetail): TargetResolverResult {
  const targets = new Set<string>();
  const context = new Set<string>();

  for (let i = 0; i < issue.variables.length; i++) {
    const variable = issue.variables[i];
    const value = i < issue.values.length ? issue.values[i] : '';

    if (!value || value === '' || value.toLowerCase() === 'not in dataset' ||
      value === 'null' || !ISO_8601_REGEX.test(String(value).trim()))
      targets.add(variable);
    else
      context.add(variable);
  }

  return ensureTargetsNotEmpty({targets, context});
}

function resolveSTRESNTargets(issue: IssueDetail): TargetResolverResult {
  const targets = new Set<string>();
  const context = new Set<string>();

  for (let i = 0; i < issue.variables.length; i++) {
    const variable = issue.variables[i];
    if (variable.toUpperCase().endsWith('STRESN'))
      targets.add(variable);
    else
      context.add(variable);
  }

  return ensureTargetsNotEmpty({targets, context});
}

const ruleTargetResolvers: {[coreId: string]: RuleTargetResolver} = {
  'CORE-000547': resolveISO8601Targets,
  'CORE-000353': resolveISO8601Targets,
  'CORE-000542': resolveSTRESNTargets,
  'CORE-000863': resolveSTRESNTargets,
};

export function resolveVariableRoles(issue: IssueDetail): TargetResolverResult {
  const resolver = ruleTargetResolvers[issue.core_id];
  if (resolver)
    return resolver(issue);

  return {targets: new Set<string>(issue.variables), context: new Set<string>()};
}
