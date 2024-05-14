import {ItemPath} from '../PipelineConfiguration';
import {normalizePaths, pathToKey} from '../config-processing-utils';


export class ControllerConfig {
  public from: ItemPath[] = [];
  public to: ItemPath[] = [];
  public fromKeys = new Set<string>();
  public toKeys = new Set<string>();

  constructor(
    public pipelinePath: ItemPath,
    from?: ItemPath | ItemPath[],
    to?: ItemPath | ItemPath[],
  ) {
    if (from) {
      this.from = normalizePaths(from);
      this.fromKeys = new Set(this.from.map((path) => pathToKey(path)));
    }

    if (to) {
      this.to = normalizePaths(to);
      this.toKeys = new Set(this.to.map((path) => pathToKey(path)));
    }
  }
}
