/** Column-mapping editor for `Chem:mpoScoreByProfile`.
 *
 *  An MPO profile scores a fixed set of PROPERTIES ("MW", "LogP", …). Which
 *  column of the incoming table supplies each one is a per-flow decision, and
 *  the interactive MPO dialog is the only place that asks — which is why the
 *  function is unusable on a canvas without this. The editor renders one column
 *  choice per profile property and stores the result as a JSON object, the
 *  shape `computeMpo` takes.
 *
 *  Two rules the UI has to obey, both learned the hard way:
 *  - **Reactive.** The profile lives in a sibling input, and the panel does not
 *    re-render while focus is inside it — so the editor subscribes to that
 *    input rather than waiting to be rebuilt.
 *  - **Never a form you can't fill.** With no table wired (or none computed
 *    yet) there are no columns to offer, and a row of empty combos reads as
 *    "broken". The editor renders only a notice in that state.
 *
 *  Auto-mapping fills BLANKS only — an explicit choice (including an explicit
 *  "don't map this") is never overwritten, so a table change tops up the gaps
 *  without undoing the user's work. */

import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import type {CustomEditorContext, CustomInputEditor} from '../../utils/func-input-overrides';
import type {FlowNode} from '../../rete/scheme';

/** The chosen profile's property list, per node, so {@link mpoMappingRequirements}
 *  can check the mapping synchronously.
 *
 *  Runtime-only and weakly held, deliberately NOT a node property: serializing
 *  it would change the saved document the moment a lookup lands, lighting up
 *  Save on a flow nobody edited. After a reload the list is simply unknown
 *  again, the validator fails closed, and the resolve below clears it. Keyed by
 *  node identity, so the several editors that can share a page never collide. */
const resolved = new WeakMap<object, {profile: string; properties: string[]}>();
/** Nodes with a lookup in flight, so a validator called once per render starts
 *  at most one fetch. */
const resolving = new WeakMap<object, string>();

/** Record a resolved list — also called by the editor, so opening the panel and
 *  the readiness check share one answer. */
export function cacheProfileProperties(node: object, profile: string, properties: string[]): void {
  resolved.set(node, {profile, properties});
}

/** Comparison key: case-insensitive, punctuation- and space-free, so
 *  `Molecular Weight`, `molecular_weight` and `MolecularWeight` all collide. */
export function normalizeName(s: string): string {
  return String(s ?? '').toLowerCase().replace(/[^a-z0-9]/g, '');
}

/** Best column for a profile property, or null when nothing is close enough.
 *
 *  In order: an exact name match, a normalized match, then a containment match
 *  (`LogP` → `Calculated LogP`) — the shortest candidate wins, since a longer
 *  name containing the property is more likely to be a different measurement
 *  (`LogP error`). Case-insensitive throughout. Pure — unit-tested. */
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

/** Fill blanks in `current` from `properties` × `columnNames`; explicit entries
 *  (including an explicit empty string) are left alone. Pure. */
export function autoMap(
  properties: readonly string[],
  columnNames: readonly string[],
  current: Record<string, string>,
): Record<string, string> {
  const next: Record<string, string> = {...current};
  const taken = new Set(Object.values(next).filter((v) => v.length > 0));
  for (const prop of properties) {
    if (prop in next) continue; // user has decided about this one
    const guess = guessColumnFor(prop, columnNames.filter((c) => !taken.has(c)));
    if (guess) {
      next[prop] = guess;
      taken.add(guess);
    }
  }
  return next;
}

/** Drop entries that aren't properties of the current profile — switching
 *  profiles must not leave the previous one's mapping behind in the JSON. */
export function pruneMapping(
  mapping: Record<string, string>, properties: readonly string[],
): Record<string, string> {
  const keep = new Set(properties);
  return Object.fromEntries(Object.entries(mapping).filter(([k]) => keep.has(k)));
}

/** Parse the stored JSON, tolerating blank / malformed values (a half-typed
 *  mapping must not throw while the panel renders). */
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

/** Drop unmapped properties before storing — `computeMpo` falls back to the
 *  property name for anything absent, so an empty entry and a missing entry
 *  mean the same thing to it, and the shorter JSON reads better in a script. */
export function serializeMapping(mapping: Record<string, string>): string {
  const filled = Object.entries(mapping).filter(([, v]) => v && v.length > 0);
  return filled.length === 0 ? '' : JSON.stringify(Object.fromEntries(filled));
}

/** MPO desirability curves are defined over numbers, so only numeric columns
 *  can supply a property — offering the rest just invites a runtime failure. */
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

/** Profile properties left without a column, given what the node has stored.
 *
 *  Every property a profile scores must be mapped: `computeMpo` silently skips
 *  the ones it can't resolve, so a partial mapping produces a score computed
 *  over fewer properties than the user asked for — a wrong number, not an
 *  error. Reported as a missing requirement, which both raises the node's
 *  "Needs input" hint and excludes it (and everything downstream) from runs.
 *
 *  Fails OPEN when the property list isn't known yet (no profile chosen, or the
 *  node has never been opened and carries no cached list): the function itself
 *  refuses to run on an incomplete mapping, so an unknown state must not block
 *  a flow that would otherwise work. Pure and synchronous. */
export function unmappedProperties(properties: readonly string[], stored: unknown): string[] {
  const mapping = parseMapping(stored);
  return properties.filter((p) => (mapping[p] ?? '').trim().length === 0);
}

/** The property list for the CURRENTLY CHOSEN profile, or null when it isn't
 *  known — in which case a resolve is started and the node stays unready until
 *  it lands. */
export function cachedProfileProperties(node: FlowNode, profileName: string): string[] | null {
  const hit = resolved.get(node);
  if (hit && hit.profile === profileName) return hit.properties;

  if (resolving.get(node) === profileName) return null; // already on its way
  resolving.set(node, profileName);
  void profileProperties(profileName).then((properties) => {
    resolving.delete(node);
    // The profile may have changed again while this was in flight.
    if (String(node.inputValues['profileName'] ?? '').trim() !== profileName) return;
    cacheProfileProperties(node, profileName, properties);
    // Recompute the hint / run gate now that the answer is known — nothing else
    // re-renders a node because a background lookup finished.
    node.editorBridge?.notifyParamsChanged(node.id);
  }).catch(() => {
    resolving.delete(node);
  });
  return null;
}

/** Readiness for an MPO node: every property the chosen profile scores must be
 *  mapped to a column.
 *
 *  Fails CLOSED while the property list is unknown — the mapping cannot be
 *  verified, and letting the node run would compute a score over whichever
 *  properties happened to resolve. The lookup above clears the state on its
 *  own, so this is a brief "resolving…" rather than a dead end. */
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

/** Profile property names, via Chem (so Flow doesn't parse the profile format).
 *  Empty on any failure — the editor then shows its "pick a profile" state. */
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
  /** Guards against an out-of-order async build: picking two profiles quickly
   *  must not let the first one's properties land after the second's. */
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

    // Resolve the properties BEFORE the no-columns check — the blocked state
    // still has to know whether the profile scores anything, and returning
    // early here is what used to skip the lookup entirely.
    const properties = await profileProperties(profileName);
    if (token !== buildToken) return; // superseded by a newer profile
    if (properties.length === 0)
      return blocked(`Profile “${profileName}” scores no properties.`);

    const columnNames = mappableColumns(ctx.columns('table'));
    if (columnNames.length === 0) {
      return blocked('Connect a table and run the flow up to this node — ' +
        'the mapping needs its numeric columns.');
    }

    // Switching profiles leaves the previous one's entries behind; drop them
    // before auto-filling, or the JSON grows stale keys.
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

  // The profile is a sibling input; the panel won't re-render while the user is
  // interacting with it, so react to the edit directly.
  ctx.watch('profileName', () => void build());

  ed.element = host;
  ed.getValue = (): unknown => serializeMapping(mapping);
  ed.setValue = (v): void => {
    mapping = parseMapping(v);
    void build();
  };
  return ed;
}
