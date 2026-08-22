import {category, expect, test, before} from '@datagrok-libraries/test/src/test';
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

    test('onLogUpdated', async () => {
        const events: any[] = [];
        const subscription = progressIndicator.onLogUpdated.subscribe((m) => events.push(m));
        try {
            const call = DG.Script.create('//language: javascript\n//output: int x\nx = 1 + 1;').prepare({});
            call.options['debug'] = true;
            await call.call(false, progressIndicator, {processed: false, report: false});
        } finally {
            subscription.unsubscribe();
        }
        expect(events.length > 0, true, 'no debug events received');
        const m = events[0];
        expect(typeof m.message, 'string', 'message');
        expect(typeof m.time, 'number', 'time');
        expect(m.level, 'info', 'level');
        expect(m.flag, 'CALL DURATION', 'flag');
        expect(m.params.source, 'Client', 'params.source');
        expect(m.params.EVENT_STAGE_KEY, 'START', 'params.EVENT_STAGE_KEY');
    });
});