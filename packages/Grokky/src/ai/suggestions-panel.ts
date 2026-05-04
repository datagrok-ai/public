import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ClaudeRuntimeClient} from '../claude/runtime-client';
import {buildViewContext} from '../claude/exec-blocks';

export const SIMPLE_REFS = ['value', 'column', 'columnName', 'table'] as const;
export type SimpleRef = typeof SIMPLE_REFS[number];
export type RefKey = SimpleRef | 'columnByName';

export type ContextRefValue = {$ref: SimpleRef} | {$ref: 'columnByName', name: string};
export type ArgValue = string | number | boolean | null | ContextRefValue;

export interface SuggestedAction {
  packageName: string;
  functionName: string;
  label: string;
  why: string;
  args: Record<string, ArgValue>;
}

export type SuggestionContext = DG.SemanticValue | DG.Column | DG.DataFrame;

const REF_DESCRIPTIONS: Record<RefKey, string> = {
  value: 'user-selected cell value',
  column: 'user-selected DG.Column',
  columnName: 'user-selected column name as a string',
  table: 'current DG.DataFrame',
  columnByName: 'DG.Column with the given name from the current table',
};

interface RefBindings {
  value?: unknown;
  column?: DG.Column;
  columnName?: string;
  table?: DG.DataFrame;
}

function buildRefs(ctx: SuggestionContext): RefBindings {
  const column = ctx instanceof DG.Column ? ctx :
    ctx instanceof DG.SemanticValue ? ctx.cell?.column :
      undefined;
  const table = column?.dataFrame ?? (ctx instanceof DG.DataFrame ? ctx : undefined);
  return {
    value: ctx instanceof DG.SemanticValue ? ctx.value : undefined,
    column,
    columnName: column?.name,
    table,
  };
}

function availableRefs(refs: RefBindings): RefKey[] {
  const simple = SIMPLE_REFS.filter((k) => refs[k] !== undefined);
  return refs.table != null ? [...simple, 'columnByName'] : simple;
}

const LITERAL_BRANCHES = [{type: 'string'}, {type: 'number'}, {type: 'boolean'}, {type: 'null'}];

const COLUMN_BY_NAME_BRANCH = {
  type: 'object',
  additionalProperties: false,
  properties: {
    $ref: {type: 'string', const: 'columnByName'},
    name: {type: 'string', minLength: 1},
  },
  required: ['$ref', 'name'],
};

function buildArgsSchema(available: RefKey[]): object {
  const simple = available.filter((r): r is SimpleRef => r !== 'columnByName');
  const description = [
    'Args for grok.functions.call. Each value is a JSON literal or a context $ref.',
    'Available refs:',
    ...available.map((r) => `  - ${r}: ${REF_DESCRIPTIONS[r]}`),
  ].join('\n');
  return {
    type: 'object',
    description,
    additionalProperties: {
      oneOf: [
        ...LITERAL_BRANCHES,
        ...(simple.length ? [{
          type: 'object',
          additionalProperties: false,
          properties: {$ref: {type: 'string', enum: simple}},
          required: ['$ref'],
        }] : []),
        ...(available.includes('columnByName') ? [COLUMN_BY_NAME_BRANCH] : []),
      ],
    },
  };
}

export function buildSuggestionsSchema(ctx: SuggestionContext): object {
  return {
    type: 'object',
    properties: {
      actions: {
        type: 'array',
        minItems: 1,
        maxItems: 4,
        items: {
          type: 'object',
          properties: {
            packageName: {type: 'string', description: 'PascalCase Datagrok package name (Chem, Admetica, Bio, …)'},
            functionName: {type: 'string', description: 'Function name as registered, e.g. drugLikeness'},
            label: {type: 'string', description: 'Imperative 3-5 word label'},
            why: {
              type: 'string',
              description: 'One short sentence grounded in the current value AND the table/columns present in ' +
                'the view (e.g. "you have a Molecule column with no ADMET annotations yet"). Avoid generic ' +
                'descriptions of what the function does.',
            },
            args: buildArgsSchema(availableRefs(buildRefs(ctx))),
          },
          required: ['packageName', 'functionName', 'label', 'why', 'args'],
        },
      },
    },
    required: ['actions'],
  };
}

function isContextRef(v: unknown): v is ContextRefValue {
  return !!v && typeof v === 'object' && '$ref' in v;
}

export function resolveArgs(args: Record<string, ArgValue>, ctx: SuggestionContext): Record<string, unknown> {
  const refs = buildRefs(ctx);
  const out: Record<string, unknown> = {};
  for (const [name, raw] of Object.entries(args ?? {})) {
    if (!isContextRef(raw)) {
      out[name] = raw;
      continue;
    }
    const v = raw.$ref === 'columnByName' ?
      refs.table?.col(raw.name) :
      refs[raw.$ref];
    if (v == null) {
      const detail = raw.$ref === 'columnByName' ? `${raw.$ref} name="${raw.name}"` : raw.$ref;
      console.warn(`AI Suggestions: unresolved $ref "${detail}" for arg "${name}"`);
      continue;
    }
    out[name] = v;
  }
  return out;
}

export function buildSuggestionsPrompt(ctx: SuggestionContext): string {
  const view = buildViewContext(grok.shell.v);
  let header: string;
  if (ctx instanceof DG.SemanticValue) {
    header = [
      `The user selected a value with semantic type "${ctx.semType ?? 'unknown'}".`,
      `Value: ${ctx.value}`,
    ].join('\n');
  } else if (ctx instanceof DG.Column) {
    header = `The user selected column "${ctx.name}" ` +
      `(type: ${ctx.type}${ctx.semType ? `, semType: ${ctx.semType}` : ''}) with ${ctx.length} rows.`;
  } else
    header = `The user selected table "${ctx.name}" with ${ctx.rowCount} rows.`;

  return [
    header,
    view ? `View: ${view}` : '',
    '',
    'Suggest the best 3-4 next actions for this selection, ranked by relevance, drawn from registered',
    'functions across available Datagrok packages (use the packages knowledge already in your',
    'context). Each action must be callable as packageName:functionName via grok.functions.call.',
    'Prefer actions that produce something immediately useful (a panel, dialog, or viewer).',
  ].filter(Boolean).join('\n');
}

export class AISuggestionsWidget extends DG.Widget {
  private readonly sessionId: string;
  private readonly ctx: SuggestionContext;

  constructor(ctx: SuggestionContext) {
    super(ui.div());
    this.ctx = ctx;
    this.sessionId = `ai-suggestions-${Date.now()}-${Math.random().toString(36).slice(2, 8)}`;
    this.root = ui.wait(() => this.render());
  }

  detach(): void {
    try {
      ClaudeRuntimeClient.getInstance().abort(this.sessionId);
    } catch { /* ignore */ }
    super.detach();
  }

  private async render(): Promise<HTMLElement> {
    try {
      const result = await ClaudeRuntimeClient.getInstance().query(
        buildSuggestionsPrompt(this.ctx),
        {sessionId: this.sessionId, outputSchema: buildSuggestionsSchema(this.ctx)});

      const actions: SuggestedAction[] = result?.actions ?? [];
      if (!actions.length)
        return ui.divText('No suggestions for this selection.');

      return ui.table(actions, (a) => [
        ui.link(a.label, () => this.invoke(a), `${a.packageName}:${a.functionName}`),
        ui.divText(a.why),
      ]);
    } catch (e: any) {
      console.warn('AI Suggestions: failed to fetch', e);
      const err = ui.divV([
        ui.divText('Could not get suggestions.'),
        ui.link('Retry', () => err.replaceWith(ui.wait(() => this.render()))),
      ]);
      return err;
    }
  }

  private async invoke(a: SuggestedAction): Promise<void> {
    const func = DG.Func.find({name: a.functionName, package: a.packageName})[0];
    if (!func) {
      grok.shell.warning(`Function not found: ${a.packageName}:${a.functionName}`);
      return;
    }
    try {
      await func.prepare(resolveArgs(a.args ?? {}, this.ctx)).call();
    } catch (e: any) {
      grok.shell.error(`Failed to run ${a.label}: ${e?.message ?? e}`);
    }
  }
}
