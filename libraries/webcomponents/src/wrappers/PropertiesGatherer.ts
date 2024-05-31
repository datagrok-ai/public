import {BehaviorSubject, Subject, merge, combineLatest} from 'rxjs';
import {filter, scan, startWith, switchMap, takeUntil, map} from 'rxjs/operators';

export type GathererArgumentsSelector<A> = (name: string | number | symbol) => name is keyof A;
export type GathererArgumentsComplete<T> = (args: Partial<T>) => boolean;

// This class is responsible for gathering web component properties
// and splitting them into two categories: constructor arguments and
// underlying instance properties.
export class PropertiesGatherer<CArgs extends object, IProps extends object, Props = CArgs & IProps> {
  private destroyed$ = new Subject;
  private propertiesUpdates$ = new Subject<{name: keyof Props, value: Props[keyof Props]}>();
  private properties$ = new BehaviorSubject<Partial<Props> | Props>({});
  private providedConstructorArguments$ = new Subject<CArgs>();
  private providedInstanceProperties$ = new Subject<Partial<IProps>>();

  public latestConstructorArguments$ = new BehaviorSubject<CArgs | undefined>(undefined);
  public latestInstanceProperties$ = new BehaviorSubject<Partial<IProps>>({});

  constructor(
    private isConstructorProp: GathererArgumentsSelector<CArgs>,
    private argsCompleteChecker?: GathererArgumentsComplete<CArgs>,
  ) {
    this.providedConstructorArguments$.pipe(
      startWith(undefined),
      switchMap((args = {} as CArgs) =>
        this.propertiesUpdates$.pipe(
          filter(({name}) => this.isConstructorProp(name)),
          scan((acc, arg) => ({...acc, ...arg}), {...args}),
        )),
      filter((args) => this.argsCompleteChecker ? this.argsCompleteChecker(args) : true),
      takeUntil(this.destroyed$),
    ).subscribe(this.latestConstructorArguments$);

    this.providedInstanceProperties$.pipe(
      startWith({} as Partial<IProps>),
      switchMap((props = {}) =>
        this.propertiesUpdates$.pipe(
          filter(({name}) => !this.isConstructorProp(name)),
          scan((acc, prop) => ({...acc, ...prop}), {...props}),
        )),
      takeUntil(this.destroyed$),
    ).subscribe(this.latestInstanceProperties$);

    combineLatest([
      merge(this.providedConstructorArguments$, this.latestConstructorArguments$),
      merge(this.providedInstanceProperties$, this.latestInstanceProperties$),
    ]).pipe(
      map(([cargs, iprops]) => ({...(cargs ? cargs : {} as CArgs), ...iprops})),
      takeUntil(this.destroyed$),
    ).subscribe(this.properties$);
  }

  setProperty(name: keyof Props, value: Props[typeof name]) {
    const payload = {name, value};
    this.propertiesUpdates$.next(payload);
  }

  getProperty(name: keyof Props) {
    return this.properties$.value?.[name];
  }

  provideInstanceData(args: CArgs | undefined, props: Partial<IProps>) {
    this.providedConstructorArguments$.next(args);
    this.providedInstanceProperties$.next(props);
  }

  destroy() {
    this.destroyed$.next(true);
  }
}
