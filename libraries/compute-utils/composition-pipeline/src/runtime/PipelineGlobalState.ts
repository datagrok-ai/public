import {Subject} from 'rxjs';

//
// Internal runtime implementation
//

export class PipelineGlobalState {
  closed = new Subject<true>();
}
