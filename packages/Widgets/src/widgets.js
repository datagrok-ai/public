/**
 * Widget demo
 * */
class TimeWidget extends DG.Widget {

  constructor() {
    super();
    this.root = ui.divV(null, 'd4-radio-button-filter');

    // properties
    this.caption = this.props.

    // external subscriptions that will be unsubscribed from in the detach() method
    this.subs = [];
  }

  attach(dataFrame) {
    this.dataFrame = dataFrame;
    this.column = DG.Utils.firstOrNull(this.dataFrame.columns.categorical);

    this.subs.push(this.dataFrame.onRowsFiltering.subscribe((_) => this.applyFilter()));

    this.render();
  }

  detach() {
    this.subs.forEach((s) => s.unsubscribe());
  }
}
