/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {PromptEngine} from './prompt-engine';
import {ChatGptFuncParams, ExecutePlanResult, ExecutionContext, FunctionMeta, PackageInfo, Plan} from './interfaces';

export class ChatGptAssistant {
  private delimiter = '####';

  constructor(private promptEngine: PromptEngine) {}

  private async chat(prompt: string, system?: string): Promise<string> {
    return await this.promptEngine.generate(
      `${this.delimiter}${prompt}${this.delimiter}`,
      system ??
      `You are an assistant that outputs only valid JSON.
      Do NOT include explanations or extra text.
      User input is delimited by ${this.delimiter}.
      Always respond strictly in JSON format.`
    );
  }

  private parseJSON<T>(text: string, fallback: T): T {
    try {
      return JSON.parse(text) as T;
    } catch {
      console.warn('Failed to parse model output as JSON:', text);
      return fallback;
    }
  }

  private getProperties(f: DG.Func): ChatGptFuncParams {
    const props: ChatGptFuncParams = {};
    for (const p of f.inputs) {
      props[p.name] = {
        type: p.propertyType,
        description: p.description ?? '',
        default: p.options['default'],
      };
    }
    return props;
  }

  // TODO: Substitute with the list of available roles along with descriptions
  private packagesInfo: PackageInfo[] = [
    {
      name: 'Admetica',
      description: 'Assessment of ADMET Properties',
    },
    {
      name: 'Bio',
      description: 'Bioinformatics support (import/export of sequences, HELM and other format conversion, visualization, analysis). [See more](https://github.com/datagrok-ai/public/blob/master/packages/Bio/README.md) for details.',
    },
    {
      name: 'Chem',
      description:
        'Cheminformatics support: import, rendering, sketching, calculation of properties, predictive models, augmentations, multiple analyses',
    },
    {
      name: 'Chembl',
      description: 'Chembl integration, commonly used queries and browser',
    },
    {
      name: 'Retrosynthesis',
      description: 'Creates retrosynthesis paths for the selected molecule',
    }
  ];

  private async selectPackages(userGoal: string): Promise<string[]> {
    const prompt = `User goal: "${userGoal}"
    Available packages:
${JSON.stringify(this.packagesInfo, null, 2)}

Determine which packages are relevant for achieving this goal.
Respond ONLY in JSON with the folllowing structure, no Markdown or extra text.
{
  "selected_packages": [
    { "name": "<package_name>", "reason": "<why this helps>" }
  ]
}

Rules:
- Include only packages that can meaningfully contribute.
- Rank them from most to least relevant.
- Do not include packages that clearly don't apply.`;

    const res = await this.chat(prompt);
    const parsed = this.parseJSON<{ selected_packages: { name: string, reason: string }[] }>(
      res,
      {selected_packages: [{name: 'Chem', reason: 'Default fallback'}]}
    );
    return parsed.selected_packages.map((p) => p.name);
  }

  private readonly reasoningSystemPrompt = `
You are a reasoning assistant that plans how to achieve user goals using available functions.

Rules:
- Output ONLY valid JSON.
- Produce the smallest possible sequence of steps.
- Each step should be:
  { "action": "call_function", "function": "<name>", "inputs": {...}, "outputs": [...] }
- ALL inputs must always appear in the "inputs" object, even if they can use a default value.
- When an input value comes from a previous step’s output, prefix it with a "$".
- If one function can directly achieve the goal, output only that.
- Do not use helper or lookup functions unless required by the goal itself.
- Include a short "analysis" explaining your reasoning.`;

  private async planFunctions(userGoal: string, functions: FunctionMeta[]): Promise<Plan> {
    const prompt = `
User request: "${userGoal}"
Available functions:
${JSON.stringify(functions, null, 2)}

Analyze what needs to be done and produce the minimal sequence of steps.
Respond strictly as JSON with this structure:
{
  "goal": "<restated_goal>",
  "analysis": ["...reasoning..."],
  "steps": [ { "action": "call_function", "function": "...", "inputs": {...}, "outputs": [...] } ]
}`;
    const result = await this.chat(prompt, this.reasoningSystemPrompt);
    return this.parseJSON(result, {goal: '', analysis: [], steps: []});
  }

  private async getFunctionsByPackage(pkg: string): Promise<FunctionMeta[]> {
    const funcs = DG.Func.find({package: pkg});
    return funcs.map((f) => {
      return {
        name: f.name,
        description: (f.nqName || f.package.name) + ': ' + (f.description || f.options['description'] || f.friendlyName || f.name),
        inputs: {type: 'object', properties: this.getProperties(f)},
        outputs: f.outputs.map((o) => ({
          name: o.name,
          type: o.propertyType,
        }))
      };
    });
  }

  public async executePlan(plan: Plan): Promise<ExecutePlanResult | null> {
    if (!plan?.steps?.length) {
      console.warn('No steps to execute in plan.');
      return null;
    }

    const context: ExecutionContext = {};
    for (const [index, step] of plan.steps.entries()) {
      if (step.action !== 'call_function') continue;

      const func = DG.Func.find({name: step.function})[0];
      if (!func) {
        console.warn(`Function not found: ${step.function}`);
        continue;
      }

      const resolvedInputs: Record<string, any> = {};
      for (const [key, val] of Object.entries(step.inputs ?? {})) {
        if (typeof val === 'string' && val.startsWith('$')) {
          const varName = val.slice(1);
          resolvedInputs[key] = context[varName]['value'];
        } else
          resolvedInputs[key] = val;
      }

      let result;
      try {
        result = await func.apply(resolvedInputs);
      } catch (err) {
        console.error(`Error executing ${func.name}:`, err);
        result = null;
      }

      if (Array.isArray(step.outputs)) {
        const outputName = step.outputs[0] ?? func.name;
        const funcOutputMeta = func.outputs.find((o) => o.name === outputName);
        context[outputName] = {value: result, meta: funcOutputMeta};
      }
    }

    const finalKey = plan.steps.at(-1)?.outputs?.[0];
    const finalResult = finalKey ? context[finalKey] : null;

    return {context, finalResult};
  }

  public async plan(userGoal: string): Promise<Plan> {
    const packages = await this.selectPackages(userGoal);
    await setTimeout(() => {}, 25000);

    const allFunctions: FunctionMeta[] = [];
    for (const pkg of packages) {
      const funcs = await this.getFunctionsByPackage(pkg);
      allFunctions.push(...funcs);
      await setTimeout(() => {}, 25000);
    }

    const plan = await this.planFunctions(userGoal, allFunctions);
    await setTimeout(() => {}, 25000);

    return plan;
  }

  public async execute(plan: Plan): Promise<ExecutePlanResult | null> {
    return await this.executePlan(plan);
  }
}
