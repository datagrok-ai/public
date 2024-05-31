import {BehaviorSubject, EMPTY, Observable, Subject} from 'rxjs';
import {PropertiesGatherer} from './PropertiesGatherer';
import {switchMap, takeUntil} from 'rxjs/operators';

export interface IEventsBase<EvType = any, EvData = any> {
  name: EvType;
  data: EvData;
}

// Base class for creating a bridge between underlying instances and a
// webcomponent properties and events.
export abstract class InstancesContainer<
Instance extends any,
CArgs extends object = any,
IProps extends object = any,
IEvents extends IEventsBase = any
> extends HTMLElement {
  private gatherer: PropertiesGatherer<CArgs, IProps>;
  private knowProps = new Set<string>();
  private instance$ = new BehaviorSubject<Instance | undefined>(undefined);

  protected destroyed$ = new Subject<true>();

  constructor() {
    super();

    this.gatherer = new PropertiesGatherer(
      this.isConstructrorProp.bind(this),
      this.argsCompleteChecker.bind(this),
    );

    this.gatherer.latestConstructorArguments$.pipe(
      takeUntil(this.destroyed$),
    ).subscribe((args) => {
      const instance = this.makeInstance(args!);
      this.instance$.next(instance);
    });

    this.gatherer.latestInstanceProperties$.pipe(
      takeUntil(this.destroyed$),
    ).subscribe((props) => {
      const instance = this.instance$.value;
      if (instance)
        this.applyInstanceProps(instance, props);
    });

    this.instance$.pipe(
      takeUntil(this.destroyed$),
    ).subscribe((instance) => {
      if (instance) {
        const allProps = this.getInstanceProperties(instance);
        const {cargs, props} = Object.entries(allProps).reduce((acc, [k, v]) => {
          if (this.isConstructrorProp(k))
            acc.cargs[k as keyof CArgs] = v as any;
          else
            acc.props[k as keyof IProps] = v as any;
          return acc;
        }, {cargs: {} as CArgs, props: {} as Partial<IProps>} as const);
        this.gatherer.provideInstanceData(cargs, props);
      } else
        this.gatherer.provideInstanceData(undefined, {});

      this.dispatchEvent(new CustomEvent('instance-changed', {detail: instance}));
    });

    this.instance$.pipe(
      switchMap((instance) => instance ? this.getInstanceEvents(instance) : EMPTY),
      takeUntil(this.destroyed$),
    ).subscribe((ev: IEvents) => {
      this.dispatchEvent(new CustomEvent(ev.name, {detail: ev.data}));
    });

    return this.makeProxy();
  }

  public setProp(name: keyof (CArgs | IProps), value: any) {
    this.gatherer.setProperty(name, value);
  }

  public getProp(name: keyof (CArgs | IProps)) {
    return this.instance$.value ? this.getInstanceProperty(this.instance$.value, name) : this.gatherer.getProperty(name);
  }

  public getCurrentMetadata() {
    const instance = this.instance$.value;
    if (instance)
      return this.getInstancePropertiesDesriptions(instance);
  }

  protected abstract isConstructrorProp(name: string | number | symbol): name is keyof CArgs;

  protected abstract argsCompleteChecker(args: Partial<CArgs>): boolean;

  protected abstract makeInstance(args: CArgs): Instance;

  protected abstract applyInstanceProps(instance: Instance, props: Partial<IProps>): void;

  protected abstract getInstancePropertiesDesriptions(instance: Instance): Record<string, any>;

  protected abstract getInstanceProperties(instance: Instance): Partial<IProps> & CArgs;

  protected abstract getInstanceProperty(instance: Instance, name: keyof (CArgs | IProps)): any;

  protected abstract getInstanceEvents(instance: Instance): Observable<IEvents>;

  private makeProxy() {
    return new Proxy(
      this,
      {
        set(target, name, value) {
          if (name === 'instance')
            target.instance$.next(value);
          else if (target.knowProps.has(name as any))
            target.setProp(name as any, value);
          else
            (target as any)[name] = value;
          return true;
        },
        get(target, name) {
          if (name === 'instance')
            return target.instance$.value;
          if (target.knowProps.has(name as any))
            return target.getProp(name as any);
          return (target as any)[name];
        },
      },
    );
  }
}
