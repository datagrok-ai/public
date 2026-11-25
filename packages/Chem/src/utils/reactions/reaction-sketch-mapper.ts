/**
 * Auto-map in Chem's reaction sketchers, the Reaction Sketcher's and the reaction cell editor's (a Map tool as in
 * ChemDraw, Auto-map and Clear all in its menu): Crux Sketch's `mapReaction` hook (crux-sketch RXN-020), answered by
 * Chem's RXNMapper on a Crux worker (reaction-mapping.ts), one step at a time, within the Reaction panel's limits: the
 * mapper's one-time load capped at MAP_LOAD_TIMEOUT_MS, each step at MAP_TIMEOUT_MS (5 s). A step the mapper does not
 * map, is late on or fails on maps nothing (null), as the hook asks; Crux Sketch then says which molecules it could not
 * pair. A mapping the mapper is not sure of (MAP-01: under 50 %) is said in a warning balloon, in the Checks' words. A
 * person's Auto-map is an edit like any other: the cell editor writes it only on OK (a mapping is never written into a
 * cell unasked, and this one was asked for).
 *
 * A vendored Crux Sketch from before its Map tool ignores the option (its createSketcher reads only the options it
 * knows), so the hook is passed already and takes effect when that sketcher is vendored.
 */
import * as grok from 'datagrok-api/grok';
import {MappedReaction, mapReaction, mapperExported, mapperReady} from './reaction-mapping';
import {MAP_LOAD_TIMEOUT_MS, MAP_TIMEOUT_MS} from '../../widgets/reaction-panel';
import {mappingFlags} from '../footguns/rules';

/** Crux Sketch's `ReactionMapper`: a step's reaction SMILES in, the same mapped (or null) out. */
export type SketchReactionMapper = (reaction: string) => Promise<string | null>;

type Mapper = (reaction: string, timeoutMs: number) => Promise<MappedReaction>;
let testMapper: Mapper | null = null;
let loading: Promise<void> | null = null;

/** Tests: the mapper the hook calls instead of Chem's, with no load; null: Chem's again. */
export function setSketchMapperForTests(f: Mapper | null): void {
  testMapper = f;
}

/** `p`'s value, or 'late' when it has not come within `ms` (a late value or error is ignored). */
async function within<T>(p: Promise<T>, ms: number): Promise<T | 'late'> {
  p.catch(() => {});
  let timer: ReturnType<typeof setTimeout> | undefined;
  const late = new Promise<'late'>((resolve) => timer = setTimeout(() => resolve('late'), ms));
  try {
    return await Promise.race([p, late]);
  } finally {
    clearTimeout(timer);
  }
}

/** One step mapped by the hook: the mapper loaded (once a page, again after a failed load), the step mapped within its
 * limit; null for anything else, the error to the console. */
async function mapStep(reaction: string): Promise<string | null> {
  try {
    if (!testMapper) {
      loading ??= mapperReady().catch((e) => {
        loading = null;
        throw e;
      });
      if (await within(loading, MAP_LOAD_TIMEOUT_MS) === 'late')
        return null;
    }
    const answer = await within((testMapper ?? mapReaction)(reaction, MAP_TIMEOUT_MS), MAP_TIMEOUT_MS);
    if (answer === 'late')
      return null;
    const unsure = mappingFlags(answer.confidence);
    if (unsure.length)
      grok.shell.warning(`Auto-map: ${unsure[0].text}`);
    return answer.mapped || null;
  } catch (e) {
    console.error('Chem | Auto-map could not map a step', e);
    return null;
  }
}

/** The hook for a reaction sketcher's options (`mapReaction`), or null when Chem's Crux file has no mapper (the Map
 * tool's menu then offers Clear all alone). */
export function sketchReactionMapper(): SketchReactionMapper | null {
  return testMapper || mapperExported() ? mapStep : null;
}
