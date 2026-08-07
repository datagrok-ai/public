/** Column-mapping editor for `Chem:mpoScoreByProfile` — one column choice per profile property,
 *  stored as the JSON object `computeMpo` takes. */

import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import type {CustomEditorContext, CustomInputEditor} from '../../utils/func-input-overrides';
import type {FlowNode} from '../../rete/scheme';

/** Chosen profile's property list per node — weakly held, deliberately NOT a node property (serializing it would light up Save on a flow nobody edited). */
const resolved = new WeakMap<object, {profile: string; properties: string[]}>();
/** Lookups in flight — a validator called on every render must start at most one fetch. */
const resolving = new WeakMap<object, string>();

export function cacheProfileProperties(node: object, profile: string, properties: string[]): void {
  resolved.set(node, {profile, properties});
}

/** Case-insensitive, punctuation- and space-free comparison key. */
export function normalizeName(s: string): string {
  return String(s ?? '').toLowerCase().replace(/[^a-z0-9]/g, '');
}

/** Best column for a property: exact, then normalized, then containment match — shortest candidate wins (a longer name is more likely a different measurement). */
export function guessColumnFor(property: string, columnNames: readonly string[]): string | null {
  const exact = columnNames.find((c) => c === property);
  if (exact) return exact;

  const target = normalizeName(property);
  if (target.length === 0) return null;

  const normal = columnNames.filter((c) => normalizeName(c) === target);
  if (normal.length > 0) return normal[0];

  const contains = columnNames
    .filter((c) => normalizeName(c).includes(target))
    .sort((a, b) => a.length - b.length || a.localeCompare(b));
  return contains[0] ?? null;
}

/** Fill blanks only — explicit entries (including an explicit empty string) are never overwritten. */
export function autoMap(
  properties: readonly string[],
  columnNames: readonly string[],
  current: Record<string, string>,
): Record<string, string> {
  const next: Record<string, string> = {...current};
  const taken = new Set(Object.values(next).filter((v) => v.length > 0));
  for (const prop of properties) {
    if (prop in next) continue;
    const guess = guessColumnFor(prop, columnNames.filter((c) => !taken.has(c)));
    if (guess) {
      next[prop] = guess;
      taken.add(guess);
    }
  }
  return next;
}

/** Drop entries that aren't properties of the current profile — a profile switch must not leave the old mapping behind. */
export function pruneMapping(
  mapping: Record<string, string>, properties: readonly string[],
): Record<string, string> {
  const keep = new Set(properties);
  return Object.fromEntries(Object.entries(mapping).filter(([k]) => keep.has(k)));
}

/** Parse the stored JSON, tolerating blank / malformed values. */
export function parseMapping(stored: unknown): Record<string, string> {
  if (stored === null || stored === undefined || String(stored).trim().length === 0) return {};
  try {
    const parsed = typeof stored === 'string' ? JSON.parse(stored) : stored;
    if (!parsed || typeof parsed !== 'object' || Array.isArray(parsed)) return {};
    const out: Record<string, string> = {};
    for (const [k, v] of Object.entries(parsed as Record<string, unknown>))
      out[k] = v === null || v === undefined ? '' : String(v);
    return out;
  } catch {
    return {};
  }
}

/** Drop unmapped properties before storing — `computeMpo` treats an empty and a missing entry the same. */
export function serializeMapping(mapping: Record<string, string>): string {
  const filled = Object.entries(mapping).filter(([, v]) => v && v.length > 0);
  return filled.length === 0 ? '' : JSON.stringify(Object.fromEntries(filled));
}

/** MPO desirability curves are defined over numbers — only numeric columns can supply a property. */
export function mappableColumns(columns: readonly DG.Column[] | null): string[] {
  if (!columns) return [];
  return columns.filter((c) => {
    try {
      return c.matches('numerical');
    } catch {
      return false;
    }
  }).map((c) => c.name);
}

/** Profile properties left without a column — a partial mapping makes `computeMpo` silently score fewer properties (a wrong number, not an error). */
export function unmappedProperties(properties: readonly string[], stored: unknown): string[] {
  const mapping = parseMapping(stored);
  return properties.filter((p) => (mapping[p] ?? '').trim().length === 0);
}

/** Property list for the chosen profile, or null while a resolve is in flight (the node stays unready until it lands). */
export function cachedProfileProperties(node: FlowNode, profileName: string): string[] | null {
  const hit = resolved.get(node);
  if (hit && hit.profile === profileName) return hit.properties;

  if (resolving.get(node) === profileName) return null;
  resolving.set(node, profileName);
  void profileProperties(profileName).then((properties) => {
    resolving.delete(node);
    // The profile may have changed again while this was in flight.
    if (String(node.inputValues['profileName'] ?? '').trim() !== profileName) return;
    cacheProfileProperties(node, profileName, properties);
    // Recompute the hint / run gate — nothing else re-renders a node because a background lookup finished.
    node.editorBridge?.notifyParamsChanged(node.id);
  }).catch(() => {
    resolving.delete(node);
  });
  return null;
}

/** Readiness for an MPO node. Fails CLOSED while the property list is unknown — running unverified would compute a wrong score; the lookup clears the state itself. */
export function mpoMappingRequirements(node: FlowNode): string[] {
  const profileName = String(node.inputValues['profileName'] ?? '').trim();
  if (!profileName)
    return []; // the profile itself is already reported missing — don't say it twice

  const properties = cachedProfileProperties(node, profileName);
  if (properties === null) return ['Column mapping (reading the profile…)'];
  if (properties.length === 0) return []; // a profile that scores nothing needs no mapping

  const unmapped = unmappedProperties(properties, node.inputValues['columnMapping']);
  return unmapped.length === 0 ? [] : [`Column mapping — unmapped: ${unmapped.join(', ')}`];
}

/** Profile property names via Chem — Flow doesn't parse the profile format. Empty on any failure. */
async function profileProperties(profileName: string): Promise<string[]> {
  if (!profileName) return [];
  try {
    const f = DG.Func.find({package: 'Chem', name: 'getMpoProfileProperties'})[0];
    if (!f) return [];
    const res = await f.apply({profileName});
    return Array.isArray(res) ? res.map((x) => String(x)) : [];
  } catch {
    return [];
  }
}

export function mpoColumnMappingEditor(param: DG.Property, ctx: CustomEditorContext): CustomInputEditor {
  const ed: CustomInputEditor = {} as CustomInputEditor;
  const host = ui.div([], 'ff-mpo-mapping');
  let mapping: Record<string, string> = {};
  /** Guards an out-of-order async build — the first profile's properties must not land after the second's. */
  let buildToken = 0;

  /** The whole editor collapses to one notice — no half-usable form. */
  const blocked = (text: string): void => {
    ui.empty(host);
    host.setAttribute('data-blocked', 'true');
    host.appendChild(ui.divV([
      ui.iconFA('link'),
      ui.divText(text, 'ff-mpo-mapping-blocked-text'),
    ], 'ff-mpo-mapping-blocked'));
  };

  const build = async (): Promise<void> => {
    const token = ++buildToken;
    const profileName = String(ctx.inputValue('profileName') ?? '');
    if (!profileName)
      return blocked('Choose a profile above — the properties it scores will be listed here.');

    // Resolve properties BEFORE the no-columns check — the blocked state still needs to know whether the profile scores anything.
    const properties = await profileProperties(profileName);
    if (token !== buildToken) return; // superseded by a newer profile
    if (properties.length === 0)
      return blocked(`Profile “${profileName}” scores no properties.`);

    const columnNames = mappableColumns(ctx.columns('table'));
    if (columnNames.length === 0) {
      return blocked('Connect a table and run the flow up to this node — ' +
        'the mapping needs its numeric columns.');
    }

    const filled = autoMap(properties, columnNames, pruneMapping(mapping, properties));
    if (serializeMapping(filled) !== serializeMapping(mapping)) {
      mapping = filled;
      ed.onChanged?.(serializeMapping(mapping));
    } else
      mapping = filled;

    ui.empty(host);
    host.removeAttribute('data-blocked');
    for (const prop of properties) {
      const input = ui.input.choice(prop, {
        value: mapping[prop] ?? '',
        items: ['', ...columnNames],
        nullable: true,
        tooltipText: `Numeric column supplying “${prop}”`,
        onValueChanged: (v) => {
          mapping[prop] = String(v ?? '');
          ed.onChanged?.(serializeMapping(mapping));
        },
      });
      input.root.setAttribute('data-mpo-property', prop);
      host.appendChild(input.root);
    }
  };

  // The profile is a sibling input and the panel won't re-render while the user is in it — react to the edit directly.
  ctx.watch('profileName', () => void build());

  ed.element = host;
  ed.getValue = (): unknown => serializeMapping(mapping);
  ed.setValue = (v): void => {
    mapping = parseMapping(v);
    void build();
  };
  return ed;
}
