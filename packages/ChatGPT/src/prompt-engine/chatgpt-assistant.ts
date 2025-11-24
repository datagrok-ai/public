/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {PromptEngine} from './prompt-engine';
import {ChatGptFuncParams, ExecutePlanResult, ExecutionContext, FunctionMeta, JsonSchema, PackageInfo, PackageSelection, PackageSelectionSchema, Plan, PlanSchema} from './interfaces';

export class ChatGptAssistant {
  private delimiter = '####';

  constructor(private promptEngine: PromptEngine) {}

  private async chat<T>(
    prompt: string,
    options: { system?: string; schema: JsonSchema }
  ): Promise<T> {
    const systemPrompt = options.system ?? `
You are an assistant that outputs only valid JSON.
Do NOT include explanations or extra text.
User input is delimited by ${this.delimiter}.
Always respond strictly in JSON format.
`;

    const raw = await this.promptEngine.generate(
      `${this.delimiter}${prompt}${this.delimiter}`,
      systemPrompt,
      options.schema
    );

    try {
      return JSON.parse(raw) as T;
    } catch (err) {
      throw new Error(`Failed to parse JSON from LLM output: ${raw}\nError: ${err}`);
    }
  }

  private getProperties(f: DG.Func): ChatGptFuncParams {
    const props: ChatGptFuncParams = {};
    for (const p of f.inputs) {
      props[p.name] = {
        type: p.propertyType,
        description: p.description ?? '',
        default: p.options['default'],
        optional: p.options['optional'],
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
        'Cheminformatics support: import, rendering, sketching, calculation of properties, predictive models, augmentations, multiple analyses. Includes multiple different calculators(descriptors, properties, solubility(LogS through properties) and many more). Can be used in conjunction to peptides, conversion tools to optimize/maximize/minimize different properties.',
    },
    {
      name: 'Chembl',
      description: 'Chembl integration, commonly used queries and browser',
    },
    {
      name: 'ClinicalCase',
      description: 'Analysis of the clinical data represented in the [SDTM](https://www.cdisc.org/standards/foundational/sdtm) and [SEND](https://www.cdisc.org/standards/foundational/send) formats.',
    },
    {
      name: 'Retrosynthesis',
      description: 'Creates retrosynthesis paths for the selected molecule',
    },
    {
      name: 'Biologics',
      description: 'Database with ADC, protein sequences, peptides, drugs and their assay data. contains tools and queries into the database. compounds identifiers for peptides, sequences, drugs, ADCs etc. like GROKPEP-######, GROKADC-###### etc.',
    },
    {
      name: 'SequenceTranslator',
      description: 'Contains tools for enumerating sequences at different positions and monomers. very useful when trying to maximize/minimize/optimize certain properties of peptides/proteins by sequence enumeration.',
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

    const res = await this.chat<PackageSelection>(prompt, {schema: PackageSelectionSchema});
    console.log('Selected packages:', res);
    return res.selected_packages.map((p: any) => p.name);
  }

  private readonly reasoningSystemPrompt = `
You are a reasoning assistant that plans how to achieve user goals using available functions.

General output rules:
- Output ONLY valid JSON.
- Produce the smallest possible sequence of steps.
- Each step should be:
  { "action": "call_function", "function": "<name>", "inputs": {...}, "outputs": [...] }

STRICT INPUT RULES (applies to ALL models):
- You MUST include EVERY REQUIRED input parameter defined for the function.
- Optional inputs:
  - If a value is provided, include it in the "inputs" object.
  - If no value is provided, the field may be omitted.
- Inputs with default values are treated as required:
  - You MUST include them in the "inputs", even if the value equals the default.
- Leaving out a required input or a defaulted input is considered an ERROR.
- Inputs should appear exactly as defined in the function metadata.

Referencing previous steps:
- If an input value comes from a previous step's output, prefix it with "$".

Planning rules:
- Produce the minimal number of steps needed to achieve the user’s goal.
- If one function can directly satisfy the goal, output only that step.
- Do not use helper or lookup functions unless required by the goal.
- Ensure all steps are valid according to the function metadata.
- Include a short "analysis" explaining your reasoning.
- When asked to optimize or maximize/minimize something, make sure that the calculation of the said property is also done.
`;

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
    return await this.chat<Plan>(prompt, {
      system: this.reasoningSystemPrompt,
      schema: PlanSchema,
    });
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
        const funcOutputMeta = func.outputs.length === 1 ? func.outputs[0] : func.outputs.find((o) => o.name === outputName);
        context[outputName] = {value: result, meta: funcOutputMeta};
      }
    }

    const finalKey = plan.steps.at(-1)?.outputs?.[0];
    const finalResult = finalKey ? context[finalKey] : null;
    if (!finalKey) {
      // some functions along the steps can be such that they just alter the dataFrame and do nothing, so search for the last dataFrame in the context
      for (let i = plan.steps.length - 2; i >= 0; i--) {
        const step = plan.steps[i];
        const outputMeta = step.outputs?.length ? context[step.outputs[0]]?.meta : null;
        if (outputMeta?.propertyType === 'dataframe')
          return {context, finalResult: context[step.outputs[0]]};
        else if (step.outputs?.length)
          break;
      }
    }

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
