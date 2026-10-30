import * as grok from 'datagrok-api/grok';
import {Control, signal, span, subscribe} from '@datagrok-libraries/u2';

const viewName = () => grok.shell.v?.name || 'No view';

class CurrentViewBadge extends Control {
  private readonly shown = signal(viewName());

  constructor() {
    super();
    this.root.append(this.build(() => span(this.shown)));
    this.subscribe(grok.events.onCurrentViewChanged, () => this.shown.value = viewName());
  }
}

export function currentViewBadge(): Control {
  return new CurrentViewBadge();
}
