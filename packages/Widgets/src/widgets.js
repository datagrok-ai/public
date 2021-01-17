/**
 * Widget demo
 * */
class TimeWidget extends DG.Widget {

  constructor() {
    super(ui.div());

    // properties
    this.caption = this.addProperty('caption', DG.TYPE.STRING, 'current time: ');

    // external subscriptions that will be unsubscribed from in the detach() method
    this.sub(rxjs.interval(1000).subscribe((_) => this.render()));

    this.render();
  }

  onPropertyChanged(property) { this.render(); }

  render() {
    $(this.root).empty();
    this.root.appendChild(ui.render([this.caption, new Date().toTimeString()]));
  }
}


class SmilesLengthWidget extends DG.Widget {
  constructor() {
    super(ui.div());

    // properties
    this.smiles = this.addProperty('smiles', DG.TYPE.STRING, null, { semType: DG.SEMTYPE.MOLECULE });
    this.render();
  }

  onPropertyChanged(property) { this.render(); }

  render() {
    this.root.innerText = `Length: ${this.smiles === null ? 'none' : this.smiles.length}`;
  }
}