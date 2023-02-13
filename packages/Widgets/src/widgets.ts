import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as rxjs from 'rxjs';
import $ from 'cash-dom';

import './styles/widgets.css';

/**
 * Widget demo
 * */
export class TimeWidget extends DG.Widget {
  caption: string;

  constructor() {
    super(ui.panel([], 'welcome-time-widget'));

    // properties
    this.caption = this.addProperty('caption', DG.TYPE.STRING, 'Current time');

    // external subscriptions that will be unsubscribed from in the detach() method
    this.sub(rxjs.interval(1000).subscribe((_) => this.render()));

    this.render();
  }

  onPropertyChanged(_: DG.Property) {this.render();}

  render() {
    $(this.root).empty();
    this.root.appendChild(ui.render([new Date().toTimeString()]));
  }
}
