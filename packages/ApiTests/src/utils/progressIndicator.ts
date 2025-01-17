import {category, expect, test, before} from '@datagrok-libraries/utils/src/test';
import * as DG from 'datagrok-api/dg';

category('ProgressIndicator', () => {
    let progressIndicator : DG.ProgressIndicator;

    before(async () => {
        progressIndicator = DG.ProgressIndicator.create();
    })

    test('updateEvent', async () => {
        let subscriptionPassed = false;
        let subscription = progressIndicator.onProgressUpdated.subscribe(() => {
            subscriptionPassed = true;
        });
        progressIndicator.update(5, 'test');
        expect(subscriptionPassed, true);
        subscription.unsubscribe();
    });
});