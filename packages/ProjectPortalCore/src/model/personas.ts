/** Simple / Advanced and quick actions per persona (F-08, FR-30–32). The mode decides how much of
 * the program page shows at once; the persona decides which quick-action cards head it. Both are
 * the user's own choice, remembered per user. A portfolio may map its groups to personas, to pick
 * the persona of a user who has not chosen one, and may declare its own cards per persona. Pure:
 * the view says which actions are possible for the program it shows. */
import type {Setting} from './config';

export type PortalMode = 'simple' | 'advanced';
export const PORTAL_MODES: PortalMode[] = ['simple', 'advanced'];

export type Persona = 'chemist' | 'biologist' | 'modeler' | 'process_scientist' | 'lead';
export const PERSONAS: Persona[] = ['chemist', 'biologist', 'modeler', 'process_scientist', 'lead'];
export const PERSONA_LABELS: Record<Persona, string> = {
  chemist: 'Chemist', biologist: 'Biologist', modeler: 'Modeler', process_scientist: 'Process scientist', lead: 'Lead',
};

/** The actions a card may open; each is something the program page already does. */
export type QuickActionCode = 'open_data' | 'open_other_data' | 'needs_attention' | 'missing_assays' |
  'suggested_models' | 'save_run' | 'choose_run' | 'follow' | 'open_workspace' | 'gate_decision' | 'portfolio_table';
export const QUICK_ACTIONS: QuickActionCode[] = ['open_data', 'open_other_data', 'needs_attention', 'missing_assays',
  'suggested_models', 'save_run', 'choose_run', 'follow', 'open_workspace', 'gate_decision', 'portfolio_table'];

/** What a card says when the portfolio does not say otherwise. */
export const QUICK_ACTION_TEXT: Record<QuickActionCode, {label: string, description: string}> = {
  open_data: {label: 'Open data', description: 'The program dataset as a table'},
  open_other_data: {label: 'Open further data', description: 'The portfolio\'s further dataset for this program'},
  needs_attention: {label: 'Needs attention', description: 'The compounds or lots the attention rules flag'},
  missing_assays: {label: 'Missing assays', description: 'What the next gate needs and rows lack'},
  suggested_models: {label: 'Suggested models', description: 'Catalogued models that fit, best first'},
  save_run: {label: 'Save a run', description: 'Save a run of this session to the program'},
  choose_run: {label: 'Choose a run', description: 'Select the run a gate criterion rests on'},
  follow: {label: 'Follow', description: 'Hear when new results land'},
  open_workspace: {label: 'Open workspace', description: 'The program\'s folder'},
  gate_decision: {label: 'Gate decision', description: 'Readiness, then advance, hold or terminate'},
  portfolio_table: {label: 'Portfolio table', description: 'Every program of the portfolio, with filters'},
};

/** The cards each persona gets unless the portfolio declares its own. */
export const DEFAULT_CARDS: Record<Persona, QuickActionCode[]> = {
  chemist: ['open_data', 'needs_attention', 'missing_assays', 'suggested_models', 'follow'],
  biologist: ['open_other_data', 'open_data', 'missing_assays', 'follow', 'open_workspace'],
  modeler: ['suggested_models', 'save_run', 'choose_run', 'missing_assays', 'open_data'],
  process_scientist: ['save_run', 'choose_run', 'open_data', 'open_workspace', 'follow'],
  lead: ['gate_decision', 'portfolio_table', 'needs_attention', 'follow', 'open_workspace'],
};

/** A card the portfolio declares that runs a function instead of a built-in action. The function
 * gets `portfolio`, `program` and `stage` (codes) when it has inputs by those names; a view it
 * returns is opened (FR-32). */
export interface FunctionCard {
  name: string;
  /** Package:Function. */
  function: string;
  description?: string;
}

export type CardConfig = QuickActionCode | FunctionCard;

export interface PersonaConfig {
  /** The persona of a user who has not chosen one and is in none of the groups below. */
  default?: Persona;
  /** Groups mapped to personas; the first group the user is in decides. */
  groups?: {group: string, persona: Persona}[];
  /** Cards per persona, in order; a persona left out keeps the defaults. */
  cards?: Partial<Record<Persona, CardConfig[]>>;
}

/** Problems with a portfolio's `personas`, as sentences. */
export function validatePersonas(p: PersonaConfig | undefined): string[] {
  if (p === undefined)
    return [];
  const errors: string[] = [];
  const isPersona = (x: unknown): x is Persona => PERSONAS.includes(x as Persona);
  if (p.default !== undefined && !isPersona(p.default))
    errors.push(`personas: default "${p.default}" is not one of ${PERSONAS.join(', ')}`);
  for (const g of p.groups ?? []) {
    if (!g?.group?.trim?.())
      errors.push('personas: a group mapping needs a group name');
    if (!isPersona(g?.persona))
      errors.push(`personas: group "${g?.group}" maps to "${g?.persona}", which is not one of ${PERSONAS.join(', ')}`);
  }
  for (const [persona, cards] of Object.entries(p.cards ?? {})) {
    if (!isPersona(persona))
      errors.push(`personas: cards for "${persona}", which is not one of ${PERSONAS.join(', ')}`);
    if (!Array.isArray(cards)) {
      errors.push(`personas: cards for ${persona} must be a list`);
      continue;
    }
    for (const c of cards) {
      if (typeof c === 'string') {
        if (!QUICK_ACTIONS.includes(c as QuickActionCode))
          errors.push(`personas: ${persona} card "${c}" is not one of ${QUICK_ACTIONS.join(', ')}`);
      } else if (!c?.name?.trim?.() || !/^[^:\s]+:[^:\s]+$/.test(c?.function ?? ''))
        errors.push(`personas: ${persona} card "${c?.name ?? '?'}" needs a name and a function as Package:Function`);
    }
  }
  return errors;
}

/** The persona of a user who has not chosen one: the first mapped group they are in, else the
 * portfolio's default, else by setting (discovery: chemist; development and rollout: process
 * scientist). */
export function defaultPersona(config: PersonaConfig | null | undefined, myGroups: string[],
  setting: Setting | string | null | undefined): Persona {
  const mine = new Set(myGroups.map((g) => g.toLowerCase()));
  const mapped = (config?.groups ?? []).find((g) => mine.has(g.group.toLowerCase()));
  if (mapped)
    return mapped.persona;
  if (config?.default)
    return config.default;
  return setting === 'discovery' ? 'chemist' : 'process_scientist';
}

/** A remembered value, or null when it is not one we know (a setting written by a newer version). */
export const asPersona = (v: string | null | undefined): Persona | null => PERSONAS.includes(v as Persona) ? v as Persona : null;
export const asMode = (v: string | null | undefined): PortalMode | null => PORTAL_MODES.includes(v as PortalMode) ? v as PortalMode : null;

/** One quick-action card as the page draws it. */
export interface QuickCard {
  /** Unique on the page: the action code, or `fn-<function>` for a function card. */
  key: string;
  action: QuickActionCode | null;
  function: string | null;
  label: string;
  description: string;
  /** Why it cannot be used for this program now; null when it can. */
  unavailable: string | null;
}

/** The cards of [persona] for one program. [unavailable] says, per action, why it cannot be used
 * now (no further dataset, no run criterion on this gate, not an approver); [labels] renames an
 * action for this program (`Open SAR table`). A card twice in the list shows once. */
export function quickCards(persona: Persona, config: PersonaConfig | null | undefined,
  unavailable: Partial<Record<QuickActionCode, string>> = {},
  labels: Partial<Record<QuickActionCode, {label?: string, description?: string}>> = {}): QuickCard[] {
  const list: CardConfig[] = config?.cards?.[persona] ?? DEFAULT_CARDS[persona];
  const seen = new Set<string>();
  const out: QuickCard[] = [];
  for (const c of list) {
    const card: QuickCard = typeof c === 'string' ?
      {key: c, action: c, function: null, label: labels[c]?.label ?? QUICK_ACTION_TEXT[c]?.label ?? c,
        description: labels[c]?.description ?? QUICK_ACTION_TEXT[c]?.description ?? '', unavailable: unavailable[c] ?? null} :
      {key: `fn-${c.function}`, action: null, function: c.function, label: c.name, description: c.description ?? c.function,
        unavailable: null};
    if (seen.has(card.key))
      continue;
    seen.add(card.key);
    out.push(card);
  }
  return out;
}

/** Program-page sections Simple mode keeps behind More: the analytics (missing assays, suggested
 * models, attention), the run list behind comparisons, the program's roles, and the gate log. */
export type PageSection = 'recommendation' | 'attention' | 'saved_runs' | 'roles' | 'snapshots' | 'history';
export const MORE_SECTIONS: PageSection[] = ['recommendation', 'attention', 'saved_runs', 'roles', 'snapshots', 'history'];

/** Board buttons Simple mode leaves out: the ones for the portfolio's keepers. */
export type BoardButton = 'sync_history' | 'export' | 'all_portfolios';
export const ADVANCED_BOARD_BUTTONS: BoardButton[] = ['sync_history', 'export', 'all_portfolios'];

/** Whether a section shows: always in Advanced; in Simple only once More is open. */
export const sectionShown = (mode: PortalMode, section: PageSection, moreOpen: boolean): boolean =>
  mode === 'advanced' || !MORE_SECTIONS.includes(section) || moreOpen;
