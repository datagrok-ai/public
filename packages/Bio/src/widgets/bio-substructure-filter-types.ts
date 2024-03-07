import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {Observable, Subject, Unsubscribable} from 'rxjs';
import {_package} from '../package';

export interface IFilterProps {
  get onChanged(): Observable<void>;

  save(): object;
  apply(propsObj: object): void;
}

/** Fasta and Helm */
export class BioFilterProps implements IFilterProps {
  private _onChanged: Subject<void> = new Subject<void>();

  get onChanged(): Observable<void> { return this._onChanged; }

  constructor(
    public substructure: string
  ) {
    return new Proxy(this, {
      set: (target: any, key: string | symbol, value: any) => {
        _package.logger.debug(`BioFilterProps.set ${key.toString()}( '${value}' )`);
        target[key] = value;
        this._onChanged.next();
        return true;
      }
    });
  }

  save(): object {
    const propsObj = {};
    for (const [key, value] of Object.entries(this)) {
      if (key !== '_onChanged') {
        // @ts-ignore
        propsObj[key] = this[key];
      }
    }
    return propsObj;
  }

  apply(propsObj: object) {
    for (const [key, value] of Object.entries(this)) {
      if (key !== '_onChanged') {
        // @ts-ignore
        this[key] = propsObj[key];
      }
    }
  }
}

export interface IBioFilter {
  get type(): string;

  get props(): IFilterProps;
  set props(value: IFilterProps);

  get onChanged(): Observable<void>;
  get filterPanel(): HTMLElement;
  get filterSummary(): string;
  get isFiltering(): boolean;

  attach(): Promise<void>;
  detach(): Promise<void>;
  resetFilter(): void;
  substructureSearch(col: DG.Column): Promise<DG.BitSet | null>;
}

/** Encapsulates input controls handling */
export abstract class BioFilterBase<TProps extends BioFilterProps> implements IBioFilter {
  abstract get filterPanel(): HTMLElement;

  abstract get emptyProps(): TProps;

  onChanged: Subject<void> = new Subject<void>();

  private _props: TProps;
  protected _propsChanging: boolean = false;
  private _propsOnChangedSub: Unsubscribable | null = null;

  abstract get type(): string;

  get props(): TProps {
    if (!this._props) this._props = this.emptyProps;
    return this._props;
  };

  set props(value: TProps) {
    this._propsChanging = true;
    try {
      if (this._propsOnChangedSub) {
        this._propsOnChangedSub.unsubscribe();
        this._propsOnChangedSub = null;
      }
      this._props = value;
      this.applyProps();
      this.onChanged.next();
      this._propsOnChangedSub = this._props.onChanged
        .subscribe(() => {
          this.onChanged.next();
        });
    } finally {
      this._propsChanging = false;
    }
  };

  abstract attach(): Promise<void>;

  async detach(): Promise<void> {
    if (this._propsOnChangedSub) {
      this._propsOnChangedSub.unsubscribe();
      this._propsOnChangedSub = null;
    }
  }

  abstract applyProps(): void;

  get filterSummary(): string { return this.props.substructure; };

  get isFiltering(): boolean { return this.props.substructure !== ''; }

  abstract substructureSearch(_column: DG.Column): Promise<DG.BitSet | null>;

  resetFilter(): void {
    this.props = this.emptyProps;
    this.onChanged.next();
  }
}
