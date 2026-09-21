/// What every extractor sees (build-plan.md WO-1): the context a build shares and the descriptor an extractor is.
/// Type-only imports, so that no extractor pulls the registry that lists it.
import type {TypeSystem} from '../types';
import type {HomeSet} from '../homes';
import type {Emitter} from './emitter';
import type {People} from './extract/process';

export type Mode = 'full' | 'public';

export interface BuildContext {
  system: TypeSystem;
  /** Absolute paths. */
  kgRoot: string;
  repoRoot: string;
  mode: Mode;
  backlogDir?: string;
  /** The marketing site's checkout (`--landing`), whose paths carry the `landing:` prefix; absent, the site is a missing source. */
  landingDir?: string;
  /** The home documents: the caller's, when check already loaded them, else loaded once per build by whichever
   * extractor asks first (`homesOf`). */
  homes?: HomeSet;
  /** The people the process layer resolves, shared with the packages extractor that runs before it (`peopleOf`). */
  people?: People;
}

export interface Extractor {
  name: string;
  /** What each source this extractor reports contributes, by source name (`emitter.source`); the manifest carries
   * the union as `provides`, so a caveat can say what is missing. */
  describes: Record<string, string>;
  modes: Mode[];
  run(ctx: BuildContext, emitter: Emitter): void | Promise<void>;
}
