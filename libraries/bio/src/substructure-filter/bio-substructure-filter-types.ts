import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {Observable, Subject, Unsubscribable} from 'rxjs';

export interface IFilterProps {
}

/** Fasta and Helm */
export class BioFilterProps implements IFilterProps {
  constructor(
    public readonly substructure: string,
    /** Pass false from an inheritors constructor, at the end set true. */ protected readOnly: boolean = true,
    protected logger?: DG.PackageLogger,
  ) {
    return new Proxy(this, {
      set: (target: any, key: string | symbol, value: any) => {
        this.logger?.debug(`BioFilterProps.set ${key.toString()}( '${value}' )`);
        if (this.readOnly)
          throw new Error('Properties are immutable.');
        target[key] = value;
        return true;
      }
    });
  }
}

export interface IBioFilter {
  get type(): string;

  get props(): IFilterProps;
  set props(value: IFilterProps);

  applyProps(props: IFilterProps): void;
  saveProps(): IFilterProps;

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

  private _props: TProps | null = null;
  protected _propsChanging: boolean = false;

  abstract get type(): string;

  get props(): TProps {
    if (!this._props) this._props = this.emptyProps;
    return this._props;
  };

  set props(value: TProps) {
    this._propsChanging = true;
    try {
      this._props = value;
      this.applyProps();
      this.onChanged.next();
    } finally {
      this._propsChanging = false;
    }
  };

  saveProps(): IFilterProps {
    const propsObj = {};
    for (const [key, value] of Object.entries(this.props)) {
      if (key !== '_onChanged' && key !== 'logger') {
        // @ts-ignore
        propsObj[key] = this.props[key];
      }
    }
    return propsObj;
  }

  abstract applyProps(): void;

  abstract attach(): Promise<void>;

  async detach(): Promise<void> { }

  get filterSummary(): string { return this.props.substructure; };

  get isFiltering(): boolean { return this.props.substructure !== ''; }

  abstract substructureSearch(_column: DG.Column): Promise<DG.BitSet | null>;

  resetFilter(): void {
    this.props = this.emptyProps;
    this.onChanged.next();
  }
}
