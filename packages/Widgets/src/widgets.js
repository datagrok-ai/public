/**
 * Widget demo
 * */
class TimeWidget extends DG.Widget {

  constructor() {
    super();
    this.root = ui.divV(null, 'd4-radio-button-filter');

    // properties
    this.caption = this.addProperty('caption', DG.TYPE.STRING, 'current time: ');

    // external subscriptions that will be unsubscribed from in the detach() method
    this.subs = [];
    this.subs.push(rxjs.interval(1000).subscribe((_) => this.render()));

    this.render();
  }

  detach() {
    grok.shell.info('timer stopped.');
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  onPropertyChanged(property) { this.render(); }

  render() {
    $(this.root).empty();
    this.root.appendChild(ui.render([this.caption, new Date().toTimeString()]));
  }
}
