import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {PipelineLinkConfigurationBase} from '../config/PipelineConfiguration';
import {LinkParsed} from '../config/LinkSpec';
import {BaseTree, NodePath} from '../data/BaseTree';
import {PipelineState} from '../config/PipelineInstance';
import {PipelineConfigurationProcessed} from '../config/config-processing-utils';
import {StateTree} from './StateTree';


export class LinksState {
  constructor(private config: PipelineConfigurationProcessed) {}

  public updateTree(state: StateTree, mutationPath: NodePath) {
    this.updateLinks(state);
  }

  private updateLinks(state: StateTree) {
    state.traverse(state.getRoot(), (acc, node, path) => {
      return acc;
    }, undefined as any);
  }
}
