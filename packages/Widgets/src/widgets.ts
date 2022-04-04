import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as rxjs from 'rxjs';
import $ from 'cash-dom';

import './styles/widget.css';

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


export class SmilesLengthWidget extends DG.Widget {
  smiles: string;

  constructor() {
    super(ui.div());

    // properties
    this.smiles = this.addProperty('smiles', DG.TYPE.STRING, null, {semType: DG.SEMTYPE.MOLECULE});
    this.render();
  }

  onPropertyChanged(_: DG.Property) {this.render();}

  render() {
    this.root.innerText = `Length: ${this.smiles === null ? 'none' : this.smiles.length}`;
  }
}
