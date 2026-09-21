/// The ordered extractors `grok kg build` runs (build-plan.md WO-1); what they share is `context.ts`.
import {Emitter} from './emitter';
import {homesExtractor} from './extract/homes';
import {packagesExtractor} from './extract/ts/packages';
import {functionsExtractor} from './extract/ts/functions';
import {testsExtractor} from './extract/ts/tests';
import {nodeTestsExtractor} from './extract/ts/node-tests';
import {samplesExtractor} from './extract/ts/samples';
import {changelogExtractor} from './extract/ts/changelog';
import {docsExtractor} from './extract/docs';
import {mediaExtractor} from './extract/media';
import {landingExtractor} from './extract/landing';
import {dartExtractor} from './extract/dart';
import {declarationsExtractor} from './extract/ts/declarations';
import {importsExtractor} from './extract/ts/imports';
import {usesExtractor} from './extract/ts/uses';
import {inlineMarkersExtractor} from './extract/ts/inline-markers';
import {processExtractor} from './extract/process';
import {membershipExtractor} from './extract/membership';
import {Mode, Extractor, BuildContext} from './context';

/** In order; membership resolution reads the claims of all the others, so it stays last. */
export const EXTRACTORS: Extractor[] = [homesExtractor, packagesExtractor, functionsExtractor, declarationsExtractor, importsExtractor, usesExtractor,
  testsExtractor, nodeTestsExtractor, samplesExtractor, changelogExtractor, docsExtractor, landingExtractor, mediaExtractor, inlineMarkersExtractor, dartExtractor, processExtractor, membershipExtractor];

/** What the sources of [extractors] contribute, by source name, for the manifest. */
export function provides(extractors: Extractor[]): Record<string, string> {
  return Object.fromEntries(extractors.flatMap((e) => Object.entries(e.describes)).sort(([a], [b]) => a < b ? -1 : a > b ? 1 : 0));
}

/** The extractors for [mode], narrowed by `--only`; names that match nothing are returned for the caller to refuse. */
export function selectExtractors(mode: Mode, only?: string[]): {selected: Extractor[], unknown: string[]} {
  const known = EXTRACTORS.filter((e) => e.modes.includes(mode));
  if (!only) return {selected: known, unknown: []};
  return {selected: known.filter((e) => only.includes(e.name)), unknown: only.filter((n) => !EXTRACTORS.some((e) => e.name === n))};
}

/** Runs the extractors in order; an extractor that reports no status is `ok`, one that throws is `missing`. */
export async function runExtractors(extractors: Extractor[], ctx: BuildContext, emitter: Emitter): Promise<void> {
  for (const extractor of extractors) {
    emitter.scope(extractor.name);
    try {
      await extractor.run(ctx, emitter);
      emitter.sources[extractor.name] ??= 'ok';
    } catch (e: any) {
      emitter.source(extractor.name, 'missing');
      emitter.problem('extractor_errors', `${extractor.name}: ${e.stack ?? e.message}`);
      console.error(`${extractor.name}: ${e.message}`);
    }
  }
  emitter.scope('');
}
