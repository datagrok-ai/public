/// Import statements resolved to files, libraries and packages (build-plan.md WO-3b): one `imports` edge per
/// (file, target) with the union of the symbols named; bare npm specifiers and asset files stay out of the graph.
import {Emitter} from '../../emitter';
import {BuildContext, Extractor} from '../../context';
import {fileId} from '../../../ids';
import {tsSources} from './declarations';

export const importsExtractor: Extractor = {
  name: 'ts-imports',
  describes: {'ts-imports': 'imports between source files'},
  modes: ['full'],
  run(ctx: BuildContext, emitter: Emitter): void {
    const sources = tsSources(ctx, emitter);
    for (const file of sources.files) {
      const from = fileId(file.path);
      const targets = new Map<string, {symbols: Set<string>, reexport: boolean}>();
      for (const imp of file.imports) {
        const {to, unresolved} = sources.resolveImport(file, imp.specifier);
        if (!to || to === from) {
          if (unresolved) emitter.problem('unresolved_ids', `${file.path}: import '${imp.specifier}' resolves to no file, library or package`);
          continue;
        }
        let target = targets.get(to);
        if (!target) targets.set(to, target = {symbols: new Set(), reexport: true});
        for (const s of imp.symbols) target.symbols.add(s);
        target.reexport = target.reexport && imp.reexport === true;
      }
      for (const [to, {symbols, reexport}] of targets)
        emitter.edge({type: 'imports', from, to, symbols: symbols.size ? [...symbols].sort() : undefined, reexport: reexport ? true : undefined, derived_by: 'ast', confidence: 1, evidence: [file.path]});
    }
    emitter.source('ts-imports', sources.failed ? 'partial' : 'ok');
  },
};
