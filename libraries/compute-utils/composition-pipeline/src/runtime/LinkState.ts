import {BehaviorSubject, Observable, merge, of, Subject} from 'rxjs';
import {withLatestFrom, filter, switchMap, mapTo} from 'rxjs/operators';
import {PipelineLinkConfiguration} from '../config/PipelineConfiguration';
import {ItemPath} from '../config/CommonTypes';
import {PipelineDriverState} from './PipelineDriverState';
import {ControllerConfig} from './ControllerConfig';

export class LinkState {
  public controllerConfig: ControllerConfig;
  public enabled = new BehaviorSubject(true);
  private currentSource = new BehaviorSubject<Observable<any>>(of());
  private externalTrigger = new Subject<true>();
  private valueChanges = this.currentSource.pipe(
    switchMap((source) => source),
  );

  constructor(
    public conf: PipelineLinkConfiguration,
    public pipelinePath: ItemPath,
    public pipelineState: PipelineDriverState,
  ) {
    this.controllerConfig = new ControllerConfig(pipelinePath, conf.from, conf.to);
  }

  public setSources(sources: Observable<any>[]) {
    this.currentSource.next(merge(...sources, this.externalTrigger));
  }

  public trigger() {
    this.externalTrigger.next(true);
  }

  public getValuesChanges() {
    return this.valueChanges.pipe(
      withLatestFrom(this.enabled),
      filter(([, enabled]) => {
        if (!enabled)
          return false;

        return true;
      }),
      mapTo(true),
    );
  }
}
