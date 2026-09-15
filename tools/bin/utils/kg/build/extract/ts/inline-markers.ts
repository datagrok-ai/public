/// Inline feature markers in TypeScript sources (conventions.md §6, build-plan.md WO-4): `// ~id` on a line of its
/// own makes the file participate in that feature, never own it. The id resolves through the home documents like
/// every other marker, and one no home declares is counted under `unresolved_ids` and draws nothing (review 3 #8).
import * as fs from 'fs';
import * as path from 'path';
import {Emitter} from '../../emitter';
import {BuildContext, Extractor} from '../../registry';
import {homesOf, resolveMention, MARKER_LINE} from '../markers';
import {tsSources} from './declarations';

export const inlineMarkersExtractor: Extractor = {
  name: 'ts-markers',
  layer: 'public',
  modes: ['full'],
  run(ctx: BuildContext, emitter: Emitter): void {
    const homes = homesOf(ctx);
    let unresolved = 0;
    for (const file of tsSources(ctx, emitter).files) {
      const lines = fs.readFileSync(path.join(ctx.repoRoot, file.path), 'utf8').split(/\r?\n/);
      const seen = new Set<string>();
      for (let i = 0; i < lines.length; i++) {
        const token = MARKER_LINE.exec(lines[i])?.[1];
        if (token === undefined || seen.has(token)) continue;
        seen.add(token);
        const target = resolveMention(emitter, homes, token, `${file.path}:${i + 1}`);
        if (!target || target.root !== 'feature') {
          unresolved++;
          continue;
        }
        emitter.claim({file: file.path, feature: target.id, rung: 1, source: 'marker', mode: 'participates', props: {}, line: i + 1});
      }
    }
    emitter.source('ts-markers', unresolved ? 'partial' : 'ok');
  },
};
