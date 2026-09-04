/* Step text → step definition, through cucumber expressions. Every definition is tried; when
   several match, the one with the fewest parameters (the most literal text) wins, then the longer
   expression; an exact tie is reported as ambiguous rather than picked silently. */
import {CucumberExpression, ParameterType, ParameterTypeRegistry} from '@cucumber/cucumber-expressions';
import {parameterTypes, ParameterTypeDef, StepDef, steps} from './registry.js';

export interface MatchedArg {
  /** Parameter type name: `string`, `int`, `element`, … */
  type: string;
  /** The text as written (quotes included for `{string}`). */
  text: string;
  /** The transformed value. */
  value: unknown;
}

export interface StepMatch {
  def: StepDef;
  args: MatchedArg[];
}

export interface MatchResult {
  match?: StepMatch;
  ambiguous?: StepMatch[];
}

interface Compiled {
  def: StepDef;
  expression: CucumberExpression;
  params: number;
}

export class StepMatcher {
  private readonly _compiled: Compiled[] = [];

  constructor(defs: StepDef[] = steps, types: ParameterTypeDef[] = parameterTypes) {
    const registry = new ParameterTypeRegistry();
    for (const t of types) {
      registry.defineParameterType(new ParameterType(t.name, t.regexp, null,
        t.transformer ?? ((s: string) => s), true, false));
    }
    for (const def of defs) {
      const expression = new CucumberExpression(def.expression, registry);
      const params = (def.expression.match(/\{[^}]+\}/g) ?? []).length;
      this._compiled.push({def, expression, params});
    }
  }

  get definitions(): StepDef[] {
    return this._compiled.map((c) => c.def);
  }

  /** Fewer parameters first, then more literal text matched (the step text minus what the
   * parameters captured); a definition that spells the words out beats one that captures them. */
  match(text: string): MatchResult {
    const hits: (StepMatch & {params: number; literal: number})[] = [];
    for (const c of this._compiled) {
      const args = c.expression.match(text);
      if (!args)
        continue;
      const matched = args.map((a) => ({
        type: a.parameterType.name ?? 'string',
        text: a.group.value ?? '',
        value: a.getValue(null),
      }));
      const literal = text.length - matched.reduce((n, a) => n + a.text.length, 0);
      hits.push({def: c.def, params: c.params, literal, args: matched});
    }
    if (hits.length === 0)
      return {};
    hits.sort((a, b) => a.params - b.params || b.literal - a.literal);
    const best = hits[0];
    const ties = hits.filter((h) => h.params === best.params && h.literal === best.literal);
    if (ties.length > 1)
      return {ambiguous: ties.map(({def, args}) => ({def, args}))};
    return {match: {def: best.def, args: best.args}};
  }

  /** Definitions sharing the longest leading word run with `text` — the "did you mean" list. */
  suggest(text: string, limit = 3): StepDef[] {
    const words = text.toLowerCase().split(/\s+/);
    const score = (def: StepDef) => {
      const own = def.expression.toLowerCase().split(/\s+/);
      let n = 0;
      while (n < words.length && n < own.length && (own[n] === words[n] || own[n].startsWith('{')))
        n++;
      return n;
    };
    return this._compiled.map((c) => c.def).map((def) => ({def, n: score(def)}))
      .filter((x) => x.n > 0).sort((a, b) => b.n - a.n).slice(0, limit).map((x) => x.def);
  }
}
