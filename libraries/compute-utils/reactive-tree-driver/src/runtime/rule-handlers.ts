import {
  Handler, MetaHandler, RuleDataEffect, RuleExpr, RuleMetaEffect, RuleSource, RuleValidatorEffect, Validator,
} from '../config/PipelineConfiguration';
import {ValidationResult} from '../data/common-types';
import {IControllerBase} from '../RuntimeControllers';
import {buildRuleContext, evaluate, isOn, ruleTargets, RuleContext} from './rule-expressions';
import {resolveSources, ValidatorVerdict} from './rule-sources';

function ruleParams<E>(controller: IControllerBase) {
  return {
    when: controller.getParam('when') as RuleExpr | undefined,
    effects: (controller.getParam('effects') ?? []) as E[],
    sources: controller.getParam('sources') as Record<string, RuleSource> | undefined,
  };
}

// rules without sources stay synchronous so they remain batchable
function withContext(
  controller: IControllerBase, sources: Record<string, RuleSource> | undefined, run: (ctx: RuleContext) => void,
) {
  const base = buildRuleContext(controller);
  if (!sources)
    return run(base);
  const resolved = resolveSources(controller, sources);
  if (resolved instanceof Promise)
    return resolved.then((values) => run({...base, ...values}));
  return run({...base, ...resolved});
}

function matchedTargets<E extends {targets: string | string[]}>(controller: IControllerBase, effect: E) {
  const outputs = controller.getMatchedOutputs();
  return ruleTargets(effect.targets).filter((target) => outputs.has(target));
}

export const ruleMetaHandler: MetaHandler = ({controller}) => {
  const {when, effects, sources} = ruleParams<RuleMetaEffect>(controller);
  return withContext(controller, sources, (ctx) => {
    const on = isOn(when, ctx);
    const metas: Record<string, Record<string, any>> = {};
    for (const effect of effects) {
      for (const target of matchedTargets(controller, effect)) {
        const meta = metas[target] ??= {};
        if (effect.effect === 'hide')
          meta.hidden = on;
        else if (effect.effect === 'show')
          meta.hidden = !on;
        else if (effect.effect === 'items' && on)
          meta.items = evaluate(effect.items, ctx);
        else if (effect.effect === 'meta' && on) {
          for (const [key, expr] of Object.entries(effect.meta))
            meta[key] = evaluate(expr, ctx);
        }
      }
    }
    for (const [target, meta] of Object.entries(metas))
      controller.setViewMeta(target, meta);
  });
};

const severityKey = {error: 'errors', warning: 'warnings', notification: 'notifications'} as const;

export const ruleValidatorHandler: Validator = ({controller}) => {
  const {when, effects, sources} = ruleParams<RuleValidatorEffect>(controller);
  return withContext(controller, sources, (ctx) => {
    const on = isOn(when, ctx);
    const results: Record<string, ValidationResult | undefined> = {};
    for (const effect of effects) {
      for (const target of matchedTargets(controller, effect)) {
        results[target] ??= undefined;
        if (!on)
          continue;
        if (effect.effect === 'verdicts') {
          const verdicts: ValidatorVerdict[] = ctx[effect.source] ?? [];
          if (!verdicts.length)
            continue;
          const result = results[target] ??= {};
          for (const {message, isError} of verdicts)
            (result[isError ? 'errors' : 'warnings'] ??= []).push({description: message});
          continue;
        }
        // an array message yields one item per element, an empty array none
        const message = evaluate(effect.message, ctx);
        const descriptions = (Array.isArray(message) ? message : [message]).map((item) => String(item));
        if (!descriptions.length)
          continue;
        const result = results[target] ??= {};
        (result[severityKey[effect.effect]] ??= []).push(...descriptions.map((description) => ({description})));
      }
    }
    for (const [target, result] of Object.entries(results))
      controller.setValidation(target, result);
  });
};

export const ruleDataHandler: Handler = ({controller}) => {
  const {when, effects, sources} = ruleParams<RuleDataEffect>(controller);
  return withContext(controller, sources, (ctx) => {
    const on = isOn(when, ctx);
    for (const effect of effects) {
      for (const target of matchedTargets(controller, effect)) {
        if (effect.effect === 'set') {
          if (on)
            controller.setAll(target, evaluate(effect.value, ctx), effect.restriction ?? 'none');
          else
            controller.clearRestriction(target);
        } else if (effect.effect === 'clear' && on)
          controller.setAll(target, null, effect.restriction ?? 'none');
      }
    }
  });
};
