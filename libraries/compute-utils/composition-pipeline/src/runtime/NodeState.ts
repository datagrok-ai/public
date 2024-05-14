import {Subject} from 'rxjs';
import {ItemName, ItemPath, PipelineStepConfiguration, PipelinePopupConfiguration, PipelineActionConfiguraion} from '../PipelineConfiguration';
import {NodeConf, NodeConfTypes} from './NodeConf';
import {NodeItemState} from './NodeItemState';
import {PipelineGlobalState} from './PipelineGlobalState';
import {ControllerConfig} from './ControllerConfig';

export class NodeState {
  public controllerConfig?: ControllerConfig;
  public states = new Map<ItemName, NodeItemState>();
  public notifier = new Subject<true>();

  constructor(
    public conf: NodeConf,
    public type: NodeConfTypes,
    public pipelinePath: ItemPath,
    public pipelineState: PipelineGlobalState,
  ) {
    if (type === 'pipeline') {
      const states = (conf as PipelineStepConfiguration).states;
      for (const state of states ?? [])
        this.states.set(state.id, new NodeItemState(state, this.pipelineState, conf.id));
    }

    if (type === 'step') {
      const states = (conf as PipelineStepConfiguration).states;
      for (const state of states ?? [])
        this.states.set(state.id, new NodeItemState(state, this.pipelineState, conf.id));
    }

    if (type === 'popup') {
      const states = (conf as PipelinePopupConfiguration).states;
      for (const state of states ?? [])
        this.states.set(state.id, new NodeItemState(state, this.pipelineState, conf.id, this.notifier));
    }

    if (type === 'action') {
      const link = conf as PipelineActionConfiguraion;
      this.controllerConfig = new ControllerConfig(pipelinePath, link.from, link.to);
    }
  }
}
