/* Gherkin → a flat model the compiler consumes: keyword types resolved (And/But inherit), outlines
   expanded into one scenario per example row, rule backgrounds folded into their scenarios. */
import {AstBuilder, GherkinClassicTokenMatcher, Parser} from '@cucumber/gherkin';
import {IdGenerator} from '@cucumber/messages';
import type * as messages from '@cucumber/messages';

export type StepType = 'Context' | 'Action' | 'Outcome';

export interface StepModel {
  keyword: string;
  type: StepType;
  text: string;
  line: number;
  table?: string[][];
  docString?: string;
}

export interface ScenarioModel {
  name: string;
  tags: string[];
  steps: StepModel[];
  line: number;
}

export interface FeatureModel {
  path: string;
  name: string;
  description: string;
  tags: string[];
  background: StepModel[];
  scenarios: ScenarioModel[];
}

export class GherkinParseError extends Error {}

export function parseFeature(path: string, source: string): FeatureModel {
  const parser = new Parser(new AstBuilder(IdGenerator.uuid()), new GherkinClassicTokenMatcher());
  let doc: messages.GherkinDocument;
  try {
    doc = parser.parse(source);
  } catch (e) {
    throw new GherkinParseError(`${path}: ${(e as Error).message}`);
  }
  const feature = doc.feature;
  const model: FeatureModel = {path, name: feature?.name ?? '', description: feature?.description?.trim() ?? '',
    tags: feature?.tags.map((t) => t.name) ?? [], background: [], scenarios: []};
  if (!feature)
    return model;
  for (const child of feature.children) {
    if (child.background)
      model.background.push(...toSteps(child.background.steps));
    if (child.scenario)
      model.scenarios.push(...expand(child.scenario, []));
    if (child.rule)
      walkRule(child.rule, model);
  }
  return model;
}

function walkRule(rule: messages.Rule, model: FeatureModel): void {
  const tags = rule.tags.map((t) => t.name);
  const background: StepModel[] = [];
  for (const child of rule.children) {
    if (child.background)
      background.push(...toSteps(child.background.steps));
    if (child.scenario) {
      for (const scenario of expand(child.scenario, tags))
        model.scenarios.push({...scenario, steps: [...background, ...scenario.steps]});
    }
  }
}

function toSteps(steps: readonly messages.Step[]): StepModel[] {
  let type: StepType = 'Context';
  return steps.map((s) => {
    const kt = s.keywordType as string | undefined;
    if (kt === 'Context' || kt === 'Action' || kt === 'Outcome')
      type = kt;
    return {keyword: s.keyword.trim(), type, text: s.text, line: s.location.line,
      table: s.dataTable ? s.dataTable.rows.map((r) => r.cells.map((c) => c.value)) : undefined,
      docString: s.docString?.content};
  });
}

function expand(scenario: messages.Scenario, inherited: string[]): ScenarioModel[] {
  const tags = [...inherited, ...scenario.tags.map((t) => t.name)];
  const base = toSteps(scenario.steps);
  if (scenario.examples.length === 0)
    return [{name: scenario.name, tags, steps: base, line: scenario.location.line}];
  const out: ScenarioModel[] = [];
  for (const examples of scenario.examples) {
    const header = examples.tableHeader?.cells.map((c) => c.value) ?? [];
    const exampleTags = examples.tags.map((t) => t.name);
    for (const row of examples.tableBody) {
      const values = row.cells.map((c) => c.value);
      const subst = (s: string) => header.reduce((acc, h, i) => acc.split(`<${h}>`).join(values[i] ?? ''), s);
      const steps = base.map((st) => ({...st, text: subst(st.text),
        table: st.table?.map((r) => r.map(subst)),
        docString: st.docString === undefined ? undefined : subst(st.docString)}));
      const label = header.map((h, i) => `${h}=${values[i] ?? ''}`).join(', ');
      out.push({name: `${subst(scenario.name)} [${label}]`, tags: [...tags, ...exampleTags], steps, line: row.location.line});
    }
  }
  return out;
}
