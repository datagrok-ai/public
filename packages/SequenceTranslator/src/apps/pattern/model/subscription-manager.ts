import * as DG from 'datagrok-api/dg';
import {Subscription} from 'rxjs';

export class SubscriptionManager {
  private rxjsSubscription = new Subscription();
  private dgSubscriptions = [] as DG.StreamSubscription[];

  constructor() { }

  add(subscription: Subscription | DG.StreamSubscription) {
    if (subscription instanceof Subscription)
      this.rxjsSubscription.add(subscription);
    else
      this.dgSubscriptions.push(subscription);
  }

  unsubscribeAll() {
    this.rxjsSubscription.unsubscribe();
    this.dgSubscriptions.forEach((s) => s.unsubscribe());
  }
}
