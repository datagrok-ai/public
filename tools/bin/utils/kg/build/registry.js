"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.EXTRACTORS = void 0;
exports.runExtractors = runExtractors;
exports.selectExtractors = selectExtractors;
var _homes = require("./extract/homes");
/// The ordered extractors `grok kg build` runs (build-plan.md WO-1) and the context they share.

const EXTRACTORS = exports.EXTRACTORS = [_homes.homesExtractor];

/** The extractors for [mode], narrowed by `--only`; names that match nothing are returned for the caller to refuse. */
function selectExtractors(mode, only) {
  const known = EXTRACTORS.filter(e => e.modes.includes(mode));
  if (!only) return {
    selected: known,
    unknown: []
  };
  return {
    selected: known.filter(e => only.includes(e.name)),
    unknown: only.filter(n => !EXTRACTORS.some(e => e.name === n))
  };
}

/** Runs the extractors in order; an extractor that reports no status is `ok`, one that throws is `missing`. */
async function runExtractors(extractors, ctx, emitter) {
  for (const extractor of extractors) {
    try {
      await extractor.run(ctx, emitter);
      emitter.sources[extractor.name] ??= 'ok';
    } catch (e) {
      emitter.source(extractor.name, 'missing');
      emitter.problem('extractor_errors', `${extractor.name}: ${e.stack ?? e.message}`);
      console.error(`${extractor.name}: ${e.message}`);
    }
  }
}