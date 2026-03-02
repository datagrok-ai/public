/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {RDModule, RDReaction} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {NamedReaction} from './types';
import {REACTIONS_STORAGE_PATH, DEFAULT_REACTIONS_FILE, USER_REACTIONS_SUFFIX} from './consts';
import {_package} from '../../package';

let _cachedReactions: NamedReaction[] | null = null;
let _cachedDefaults: NamedReaction[] | null = null;

/** Invalidate the cached reaction list so next load re-reads from storage. */
export function invalidateReactionsCache(): void {
  _cachedReactions = null;
}

/**
 * Load all reactions: defaults merged with user-defined.
 * Results are cached; call invalidateReactionsCache() after mutations.
 */
export async function loadAllReactions(): Promise<NamedReaction[]> {
  if (_cachedReactions)
    return _cachedReactions;

  const defaults = await loadDefaultReactions();
  let userReactions: NamedReaction[] = [];
  try {
    userReactions = await loadUserReactions();
  } catch (e) {
    console.warn('Chem | Reactions: failed to load user reactions', e);
  }

  // user reactions override defaults by id
  const byId = new Map<string, NamedReaction>();
  for (const r of defaults) byId.set(r.id, r);
  for (const r of userReactions) byId.set(r.id, {...r, isUserDefined: true});
  _cachedReactions = Array.from(byId.values());
  return _cachedReactions;
}

/** Load the built-in default reactions from the shipped JSON file (files/reactions/default-reactions.json). */
export async function loadDefaultReactions(): Promise<NamedReaction[]> {
  if (_cachedDefaults)
    return _cachedDefaults;

  try {
    const text = await _package.files.readAsText(DEFAULT_REACTIONS_FILE);
    const parsed = JSON.parse(text);
    if (Array.isArray(parsed)) {
      _cachedDefaults = parsed.map((r: any) => ({...r, isUserDefined: false}));
      return _cachedDefaults;
    }
  } catch (e) {
    console.warn('Chem | Reactions: failed to load default reactions from package file', e);
  }
  return [];
}

/** Load user-defined reactions for the current user from AppData. */
export async function loadUserReactions(): Promise<NamedReaction[]> {
  const userName = getUserFileName();
  const filePath = `${REACTIONS_STORAGE_PATH}/${userName}${USER_REACTIONS_SUFFIX}`;
  try {
    const exists = await grok.dapi.files.exists(filePath);
    if (!exists)
      return [];
    const text = await grok.dapi.files.readAsText(filePath);
    const parsed = JSON.parse(text) as NamedReaction[];
    return Array.isArray(parsed) ? parsed.map((r) => ({...r, isUserDefined: true})) : [];
  } catch (e) {
    console.warn(`Chem | Reactions: error reading user reactions from ${filePath}`, e);
    return [];
  }
}

/** Save a user-defined reaction (add or update). */
export async function saveUserReaction(reaction: NamedReaction): Promise<void> {
  const userReactions = await loadUserReactions();
  const idx = userReactions.findIndex((r) => r.id === reaction.id);
  const reactionToSave: NamedReaction = {...reaction, isUserDefined: true, author: reaction.author ?? DG.User.current().friendlyName};
  if (idx >= 0)
    userReactions[idx] = reactionToSave;
  else
    userReactions.push(reactionToSave);
  await writeUserReactions(userReactions);
  invalidateReactionsCache();
}

/** Delete a user-defined reaction by ID. Returns true if found and deleted. */
export async function deleteUserReaction(reactionId: string): Promise<boolean> {
  const userReactions = await loadUserReactions();
  const idx = userReactions.findIndex((r) => r.id === reactionId);
  if (idx < 0)
    return false;
  userReactions.splice(idx, 1);
  await writeUserReactions(userReactions);
  invalidateReactionsCache();
  return true;
}

/** Validate a reaction SMARTS/SMIRKS string using RDKit. */
export function validateReactionSmarts(
  rdkit: RDModule, smarts: string,
): {valid: boolean; error?: string; numReactants?: number; numProducts?: number} {
  // basic format check first
  const parts = smarts.split('>>');
  if (parts.length !== 2)
    return {valid: false, error: 'Reaction SMARTS must contain exactly one ">>".'};

  let rxn: RDReaction | null = null;
  try {
    rxn = rdkit.get_rxn(smarts);
    if (!rxn)
      return {valid: false, error: 'RDKit could not parse the reaction SMARTS.'};
    const reactantCount = parts[0].split('.').length;
    const productCount = parts[1].split('.').length;
    return {valid: true, numReactants: reactantCount, numProducts: productCount};
  } catch (e: any) {
    return {valid: false, error: e?.message ?? 'Invalid reaction SMARTS.'};
  } finally {
    try {rxn?.delete();} catch (_) {/* noop */}
  }
}

/** Export a list of reactions to a JSON string for sharing. */
export function exportReactions(reactions: NamedReaction[]): string {
  const cleaned = reactions.map((r) => ({...r, isUserDefined: undefined}));
  return JSON.stringify(cleaned, null, 2);
}

/** Import reactions from a JSON string. Validates each one. Returns only valid reactions. */
export function importReactions(rdkit: RDModule, json: string): {valid: NamedReaction[]; errors: string[]} {
  const errors: string[] = [];
  let parsed: any[];
  try {
    parsed = JSON.parse(json);
    if (!Array.isArray(parsed))
      return {valid: [], errors: ['JSON must be an array of reaction objects.']};
  } catch (e: any) {
    return {valid: [], errors: [`Invalid JSON: ${e?.message}`]};
  }

  const valid: NamedReaction[] = [];
  for (let i = 0; i < parsed.length; i++) {
    const r = parsed[i];
    if (!r.id || !r.name || !r.reactionSmarts || !r.mode || !r.category) {
      errors.push(`Reaction at index ${i} is missing required fields (id, name, reactionSmarts, mode, category).`);
      continue;
    }
    const validation = validateReactionSmarts(rdkit, r.reactionSmarts);
    if (!validation.valid) {
      errors.push(`Reaction "${r.name}" (index ${i}): ${validation.error}`);
      continue;
    }
    valid.push({...r, isUserDefined: true});
  }
  return {valid, errors};
}

/** Generate a URL-safe slug ID from a reaction name. */
export function generateReactionId(name: string): string {
  const base = name.toLowerCase()
    .replace(/[^a-z0-9]+/g, '-')
    .replace(/^-|-$/g, '');
  return `${base}-${Date.now().toString(36)}`;
}

// ---- internal helpers ----

function getUserFileName(): string {
  try {
    return DG.User.current().login.replace(/[^a-zA-Z0-9_-]/g, '_');
  } catch {
    return 'default';
  }
}

async function writeUserReactions(reactions: NamedReaction[]): Promise<void> {
  const userName = getUserFileName();
  const filePath = `${REACTIONS_STORAGE_PATH}/${userName}${USER_REACTIONS_SUFFIX}`;
  const content = JSON.stringify(reactions, null, 2);
  await grok.dapi.files.writeAsText(filePath, content);
}
