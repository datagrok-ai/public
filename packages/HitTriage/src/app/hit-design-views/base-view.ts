import * as DG from 'datagrok-api/dg';
import {divText} from 'datagrok-api/ui';
import {HitDesignApp} from '../hit-design-app';
import {HitDesignTemplate} from '../types';
export class HitDesignBaseView extends DG.ViewBase {
  app: HitDesignApp;

  constructor(app: HitDesignApp) {
    super();

    this.app = app;

    this.root.classList.add('grok-hit-triage-view');
    this.root.style.display = 'flex';
    this.statusBarPanels = [divText('Hit Design')];
  }

  get template(): HitDesignTemplate {return this.app.template!;}

  /** Override to initialize the view based on the session. */
  onActivated(): void {
  }

  /** Gets called when any of the previous views change. */
  onReset() {
  }

  async process(): Promise<any> {
  }
}
