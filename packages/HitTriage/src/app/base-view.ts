import * as DG from 'datagrok-api/dg';
import {HitDesignTemplate, HitTriageTemplate} from './types';
import {HitDesignApp} from './hit-design-app';
import {HitTriageApp} from './hit-triage-app';

export class HitBaseView<Ttemplate extends HitDesignTemplate | HitTriageTemplate,
    Tapp extends HitDesignApp | HitTriageApp> extends DG.ViewBase {
  app: Tapp;

  constructor(app: Tapp) {
    super();

    this.app = app;

    this.root.classList.add('grok-hit-triage-view');
    this.root.style.display = 'flex';
  }

  get template(): Ttemplate {return this.app.template! as Ttemplate;}

  /** Override to initialize the view based on the session. */
  onActivated(): void {
  }

  /** Gets called when any of the previous views change. */
  onReset() {
  }

  async process(): Promise<any> {
  }
}
