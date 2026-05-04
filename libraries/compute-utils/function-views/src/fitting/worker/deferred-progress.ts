import * as DG from 'datagrok-api/dg';

// A delayed `DG.TaskBarProgressIndicator`: the real indicator is created
// only after `delayMs` has elapsed. If `close()` is called sooner, no UI
// is ever shown — avoids flicker on sub-`delayMs` fits. API matches the
// subset the worker executor uses: `update`, `close`, `canceled`.

export interface DeferredProgressOptions {
  cancelable?: boolean;
  /** Default: 500. */
  delayMs?: number;
}

export class DeferredProgressIndicator {
  private real: DG.TaskBarProgressIndicator | null = null;
  private timer: ReturnType<typeof setTimeout> | null;
  private closed = false;
  private latestPercent = 0;
  private latestDescription: string;

  static create(name: string, opts: DeferredProgressOptions = {}): DeferredProgressIndicator {
    return new DeferredProgressIndicator(name, opts);
  }

  private constructor(
    private readonly name: string,
    private readonly opts: DeferredProgressOptions,
  ) {
    this.latestDescription = name;
    this.timer = setTimeout(() => this.materialize(), opts.delayMs ?? 500);
  }

  private materialize(): void {
    this.timer = null;
    if (this.closed) return;
    this.real = DG.TaskBarProgressIndicator.create(this.name,
      {cancelable: this.opts.cancelable ?? false});
    this.real.update(this.latestPercent, this.latestDescription);
  }

  get canceled(): boolean { return this.real?.canceled ?? false; }

  update(percent: number, description: string): void {
    this.latestPercent = percent;
    this.latestDescription = description;
    this.real?.update(percent, description);
  }

  close(): void {
    if (this.closed) return;
    this.closed = true;
    if (this.timer != null) {
      clearTimeout(this.timer);
      this.timer = null;
      return;
    }
    this.real?.close();
  }
}
