import {Subject} from 'rxjs';
import {PipelineStepConfiguration, PipelineActionConfiguraion} from '../config/PipelineConfiguration';
import {ItemName, ItemPath} from '../config/CommonTypes';
import {NodeConf, NodeConfTypes} from './NodeConf';
import {NodeIOState} from './NodeIOState';
import {PipelineDriverState} from './PipelineDriverState';
import {ControllerConfig} from './ControllerConfig';

export class NodeState {
  public controllerConfig?: ControllerConfig;
  public states = new Map<ItemName, NodeIOState>();
  public notifier = new Subject<true>();

  constructor(
    public conf: NodeConf,
    public type: NodeConfTypes,
    public pipelinePath: ItemPath,
    public pipelineState: PipelineDriverState,
  ) {
    if (type === 'pipeline') {
      const states = (conf as PipelineStepConfiguration).states;
      for (const state of states ?? [])
        this.states.set(state.id, new NodeIOState(state, this.pipelineState, conf.id));
    }

    if (type === 'step') {
      const states = (conf as PipelineStepConfiguration).states;
      for (const state of states ?? [])
        this.states.set(state.id, new NodeIOState(state, this.pipelineState, conf.id));
    }

    if (type === 'action') {
      const link = conf as PipelineActionConfiguraion;
      this.controllerConfig = new ControllerConfig(pipelinePath, link.from, link.to);
    }
  }
}
