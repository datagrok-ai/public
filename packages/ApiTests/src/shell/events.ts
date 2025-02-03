import { before, category, delay, expect, expectArray, test } from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

category('Events', () => {
    test('onContextMenu', async () => {
        let subscriptionPassed = false;
        let tv = grok.shell.addTableView(grok.data.demo.demog())
        let barChart = tv.barChart();
        let subscription = barChart.onContextMenu.subscribe(menu => subscriptionPassed = true);
        const rightClickEvent = new MouseEvent("contextmenu", {
            bubbles: true,
            cancelable: true,
            clientX: 100,
            clientY: 100,
        });
        barChart.root.dispatchEvent(rightClickEvent);
        await delay(500);
        subscription.unsubscribe();
        expect(subscriptionPassed, true);
    })

    test('onEvent', async () => {
        let subscriptionPassed = false;
        let subscription = grok.events.onCustomEvent("test").subscribe((e) => { subscriptionPassed = true; });
        grok.events.fireCustomEvent("test", null);
        await delay(500);
        subscription.unsubscribe();
        expect(subscriptionPassed, true);
    });

    test('onCustomEvent/fireCustomEvent', async () => {
        let subscriptionPassed = false;
        let subscription = grok.events.onCustomEvent("test").subscribe((e) => { subscriptionPassed = true; });
        grok.events.fireCustomEvent("test", null);
        await delay(500);
        subscription.unsubscribe();
        expect(subscriptionPassed, true);
    });

    test('onTableAdded', async () => {
        let subscriptionPassed = false;
        let subscription = grok.events.onTableAdded.subscribe(() => { subscriptionPassed = true; });
        grok.shell.addTable(grok.data.demo.demog(100));
        await delay(500);
        subscription.unsubscribe();
        expect(subscriptionPassed, true);
    });

    test('onTableRemoved', async () => {
        let subscriptionPassed = false;
        let subscription = grok.events.onTableRemoved.subscribe(() => { subscriptionPassed = true; });
        let table = grok.shell.addTable(grok.data.demo.demog(100));
        grok.shell.closeTable(table);
        await delay(500);
        subscription.unsubscribe();
        expect(subscriptionPassed, true);
    });

    test('debounce', async () => {
        let subscriptionPassed = false;
        let tv = grok.shell.addTableView(grok.data.demo.demog())
        let subscription = DG.debounce(tv.dataFrame.onColumnNameChanged, 500).subscribe(() => {
            subscriptionPassed = true;
        });

        tv.dataFrame.columns.toList()[0].name = "debounceTest"
        await delay(2000);
        subscription.unsubscribe();
        expect(subscriptionPassed, true);
    });
});