import {Subject} from 'rxjs';

export class PipelineDriverState {
  public pipelineConfig: Pipe;

  closed = new Subject<true>();
}
