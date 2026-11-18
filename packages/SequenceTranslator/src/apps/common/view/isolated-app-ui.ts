/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import type {Subscription} from 'rxjs';

import {_package} from '../../../package';
import {AppUIBase} from './app-ui-base';

export abstract class IsolatedAppUIBase extends AppUIBase {
  constructor(appName: string) {
    super(appName);
    this.view = DG.View.create();
    this.configureView();
    this.releaseSubscriptionsWith(this.view);
  }

  protected view: DG.View;
  private removalWatch?: Subscription;

  /** The subscriptions of [view]'s `subs` (the app's own, as its table controls' on the workspace's tables) end when
   * [shown] leaves the shell: this view, or the toolkit view it is a tab of. A DG.View does not end them itself. */
  releaseSubscriptionsWith(shown: DG.ViewBase): void {
    this.removalWatch?.unsubscribe();
    this.removalWatch = grok.events.onViewRemoved.subscribe((removed) => {
      if (removed?.dart !== shown.dart)
        return;
      this.removalWatch?.unsubscribe();
      this.removalWatch = undefined;
      for (const sub of this.view.subs.splice(0))
        sub.unsubscribe();
    });
  }
  async constructView(): Promise<DG.ViewBase> {
    await this.initView();
    const name = this.parentAppName ? this.parentAppName + '/' + this.appName : this.appName;
    this.view.path = `/apps/${_package.name}/${name.replace(/\s/g, '')}/`;
    return this.view;
  }

  protected abstract getContent(): Promise<HTMLDivElement>;
  async initView(): Promise<void> {
    const content = await this.getContent();
    this.view.append(content);
  }

  protected configureView(): void {
    this.view.box = true;
    this.view.name = this.appName;

    const windows = grok.shell.windows;
    windows.showContextPanel = false;
    windows.showToolbox = false;
    windows.showHelp = false;
  }

  getView(): DG.View {
    return this.view;
  }
}

