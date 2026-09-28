import * as DG from 'datagrok-api/dg';
import {Observable, Subject, Subscription} from 'rxjs';

export class ViewerRenderState {
  readonly rendered = new Subject<void>();
  private readonly children = new Map<DG.Viewer, Subscription>();
  private readonly timers = new Map<string, ReturnType<typeof setTimeout>>();
  private running = 0;
  private immediate = false;
  private disposed = false;

  get isPending(): boolean {
    return this.timers.size > 0 || this.running > 0 ||
      Array.from(this.children.keys()).some((viewer) => viewer.root.isConnected && (viewer as any).isRenderPending);
  }

  get immediateRendering(): boolean { return this.immediate; }
  set immediateRendering(value: boolean) {
    this.immediate = value;
    for (const viewer of this.children.keys())
      (viewer as any).immediateRendering = value;
  }

  add(viewer: DG.Viewer, rendered: Observable<unknown> = (viewer as any).onAfterDrawScene): void {
    (viewer as any).immediateRendering = this.immediate;
    this.children.set(viewer, rendered.subscribe(() => this.rendered.next()));
  }

  remove(viewer: DG.Viewer): void {
    this.children.get(viewer)?.unsubscribe();
    this.children.delete(viewer);
    viewer.detach();
  }

  defer(key: string, action: () => void | Promise<void>, delay: number): void {
    this.cancel(key);
    if (this.disposed)
      return;
    this.timers.set(key, setTimeout(() => {
      this.timers.delete(key);
      void this.run(action);
    }, this.immediate ? 0 : delay));
  }

  cancel(key: string): void {
    clearTimeout(this.timers.get(key));
    this.timers.delete(key);
  }

  async run(action: () => void | Promise<void>): Promise<void> {
    this.running++;
    try {
      await action();
    } finally {
      this.running--;
      this.rendered.next();
    }
  }

  dispose(): void {
    this.disposed = true;
    for (const key of this.timers.keys())
      this.cancel(key);
    for (const viewer of this.children.keys())
      this.remove(viewer);
    this.rendered.complete();
  }
}
