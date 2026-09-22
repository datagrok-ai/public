import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {Observable, Subject} from 'rxjs';

export class SearchBaseViewer extends DG.JsViewer {
  private _name: string = '';
  get name(): string { return this._name; }
  set name(x: string | undefined) { this._name = x ?? ''; }
  semType: string = '';
  limit: number;
  targetColumn?: DG.Column<string>;
  targetColumnName: string;
  initialized: boolean = false;
  gridSelect: boolean = false;
  protected maxLimit: number = 100;
  protected recomputeOnCurrentRowChange: boolean = true;
  protected skipRecomputingProperies: string[] = [];
  private _renderPending = 0;
  private _onRendered = new Subject<void>();

  /** Fires when the render queue drains — what automation settles on with [isRenderPending]. */
  get onRendered(): Observable<void> {return this._onRendered;}

  /** Whether a queued render still owes a frame. The viewer renders through a promise chain and
   * refingerprints off the main thread, so its DOM being non-empty says nothing about the frame. */
  get isRenderPending(): boolean {return this._renderPending > 0;}

  constructor(name: string, semType: string) {
    super();
    this.limit = this.int('limit', 10, {min: 1, max: this.maxLimit});
    this.targetColumnName = this.string('targetColumnName', null, {...(semType ? {semType: semType} : {}), nullable: false});
    this.name = name;
    this.semType = semType;
  }

  init(): void {
    this.initialized = true;
  }

  detach(): void {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  async onTableAttached(): Promise<void> {
    this.init();

    if (this.dataFrame) {
      this.subs.push(DG.debounce(this.dataFrame.onRowsRemoved, 50)
        .subscribe((_: any) => this.render(true)));
      const compute = this.name !== 'diversity';
      this.subs.push(DG.debounce(this.dataFrame.onCurrentRowChanged, 50)
        .subscribe((_: any) => {
          if (!this.gridSelect)
            this.render(compute);
        }));
      this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50)
        .subscribe((_: any) => this.render(false)));
      this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50)
        .subscribe((_: any) => this.render(false)));
      this.targetColumnName ??= this.dataFrame.columns.bySemType(this.semType)?.name ?? '';
      this.targetColumn = this.targetColumnName ? this.dataFrame.col(this.targetColumnName) ?? undefined : undefined;
      this.getProperty('limit')!.fromOptions({min: 1, max: this.maxLimit});
    }
    this.render();
  }

  onPropertyChanged(property: DG.Property): void {
    super.onPropertyChanged(property);
    if (!this.initialized)
      return;
    if (property.name === 'targetColumnName') {
      const col = this.dataFrame.col(property.get(this))!;
      if (col.semType === this.semType)
        this.targetColumn = col;
    }
    this.debouncedRender();
  }

  private debounceTimer: any = null;
  private debouncedRender(computeData = true) {
    if (this.debounceTimer)
      clearTimeout(this.debounceTimer);
    this.debounceTimer = setTimeout(() => {
      this.render(computeData);
      this.debounceTimer = null;
    }, 100);
  }

  /** For tests */ public computeRequested: boolean = false;
  public renderPromise: Promise<void> = Promise.resolve();

  protected render(computeData = true): void {
    this._renderPending++;
    const settle = () => {
      if (--this._renderPending === 0)
        this._onRendered.next();
    };
    this.renderPromise = this.renderPromise.then(async () => {
      if (this.dataFrame && !this.targetColumn) {
        ui.empty(this.root);
        this.root.appendChild(ui.divText(`No ${this.semType} column available in the table.`, 'd4-viewer-error'));
        return;
      }
      this.computeRequested = this.computeRequested || computeData;
      await this.renderInt(computeData);
    // settle on both outcomes, or one failed render leaves the viewer pending for good — and a
    // rejected chain would skip every later render, so the failure is logged rather than rethrown
    }).then(settle, (e) => {
      settle();
      console.error(e);
    });
  }

  async renderInt(_computeData: boolean): Promise<void> {

  }

  beforeRender() {
    if (!this.initialized || !this.targetColumnName)
      return false;
    if (this.dataFrame && this.targetColumnName &&
          this.dataFrame.col(this.targetColumnName)!.semType !== this.semType) {
      grok.shell.error(`${this.targetColumnName} is not ${this.semType} type`);
      return false;
    }
    return true;
  }
}

/** What a search viewer's result cards report about themselves — the `search-results` status a
 * viewer registers with `addStatusProvider`. The viewer marks each card with the row it shows
 * (`data-row`), and `d4-current` / `d4-selected` carry its state, so the areas and the counts hold
 * for any card layout. A card with no row of its own, such as the sketched reference molecule, is
 * not reported. */
export function searchResultsStatus(root: HTMLElement, targetColumn: string, limit: number):
  {hitAreas: {[name: string]: DG.IRectBounds}, values: {[name: string]: number | string}} {
  const origin = root.getBoundingClientRect();
  const hitAreas: {[name: string]: DG.IRectBounds} = {};
  let current = -1;
  let selected = 0;
  const cards = Array.from(root.querySelectorAll('[data-row]')) as HTMLElement[];
  for (const card of cards) {
    const row = Number(card.getAttribute('data-row'));
    const box = card.getBoundingClientRect();
    hitAreas[`card ${row}`] = {
      x: box.left - origin.left, y: box.top - origin.top, width: box.width, height: box.height,
    };
    // the similarity viewer marks both the current row and the reference molecule `d4-current`,
    // so this reports the first of them rather than whichever came last
    if (current === -1 && card.classList.contains('d4-current'))
      current = row;
    if (card.classList.contains('d4-selected'))
      selected++;
  }
  return {
    hitAreas,
    values: {
      'target column': targetColumn ?? '',
      'limit': limit,
      'cards': cards.length,
      'current card': current,
      'selected cards': selected,
    },
  };
}
