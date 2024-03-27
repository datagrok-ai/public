import * as DG from 'datagrok-api/dg';
import {Subscription} from 'rxjs';

export class SubscriptionManager {
  private rxjsSubscriptions = [] as Subscription[];
  private dgSubscriptions = [] as DG.StreamSubscription[];

  constructor() { }

  add(subscription: Subscription | DG.StreamSubscription) {
    if (subscription instanceof Subscription)
      this.rxjsSubscriptions.push(subscription);
    else
      this.dgSubscriptions.push(subscription);
  }

  unsubscribeAll() {
    for (const subs of [this.rxjsSubscriptions, this.dgSubscriptions]) {
      subs.forEach((sub) => sub.unsubscribe());
      subs.length = 0;
    }
  }
}
