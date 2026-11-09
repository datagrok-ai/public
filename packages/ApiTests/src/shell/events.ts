import { category, test, waitUntil } from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

category('Events', () => {
    test('onContextMenu', async () => {
        let subscriptionPassed = false;
        let tv = grok.shell.addTableView(grok.data.demo.demog())
        let barChart = tv.barChart();
        let subscription = barChart.onContextMenu.subscribe((menu : DG.Menu) => subscriptionPassed = true);
        const rightClickEvent = new MouseEvent("contextmenu", {
            bubbles: true,
            cancelable: true,
            clientX: 100,
            clientY: 100,
        });
        try {
            barChart.root.dispatchEvent(rightClickEvent);
            await waitUntil(() => subscriptionPassed, 'onContextMenu did not fire');
        } finally {
            subscription.unsubscribe();
            tv.close();
        }
    })

    test('onContextMenuShown / onContextMenuClosed', async () => {
        let shown = false;
        let closed = false;
        let tv = grok.shell.addTableView(grok.data.demo.demog());
        let barChart = tv.barChart();
        let subShown = grok.events.onContextMenuShown.subscribe(() => {
            shown = document.querySelectorAll('.d4-menu-popup').length > 0;
        });
        let subClosed = grok.events.onContextMenuClosed.subscribe(() => closed = true);
        try {
            barChart.root.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, clientX: 100, clientY: 100}));
            await waitUntil(() => shown, 'onContextMenuShown did not fire with a menu on screen');
            document.body.click();
            await waitUntil(() => closed, 'onContextMenuClosed did not fire');
        } finally {
            subShown.unsubscribe();
            subClosed.unsubscribe();
            tv.close();
        }
    });

    test('onEvent', async () => {
        let subscriptionPassed = false;
        let subscription = grok.events.onCustomEvent("test").subscribe((e: any) => { subscriptionPassed = true; });
        try {
            grok.events.fireCustomEvent("test", null);
            await waitUntil(() => subscriptionPassed, 'onCustomEvent did not fire');
        } finally {
            subscription.unsubscribe();
        }
    });

    test('onCustomEvent/fireCustomEvent', async () => {
        let subscriptionPassed = false;
        let subscription = grok.events.onCustomEvent("test").subscribe((e: any) => { subscriptionPassed = true; });
        try {
            grok.events.fireCustomEvent("test", null);
            await waitUntil(() => subscriptionPassed, 'onCustomEvent did not fire');
        } finally {
            subscription.unsubscribe();
        }
    });

    test('onTableAdded', async () => {
        let subscriptionPassed = false;
        let subscription = grok.events.onTableAdded.subscribe(() => { subscriptionPassed = true; });
        try {
            grok.shell.addTable(grok.data.demo.demog(100));
            await waitUntil(() => subscriptionPassed, 'onTableAdded did not fire');
        } finally {
            subscription.unsubscribe();
        }
    });

    test('onTableRemoved', async () => {
        let subscriptionPassed = false;
        let subscription = grok.events.onTableRemoved.subscribe(() => { subscriptionPassed = true; });
        try {
            let table = grok.shell.addTable(grok.data.demo.demog(100));
            grok.shell.closeTable(table);
            await waitUntil(() => subscriptionPassed, 'onTableRemoved did not fire');
        } finally {
            subscription.unsubscribe();
        }
    });

    test('onPanelVisibilityChanged', async () => {
        let fired = false;
        const initial = grok.shell.windows.showStatusBar;
        let subscription = grok.shell.windows.onPanelVisibilityChanged.subscribe(() => { fired = true; });
        try {
            grok.shell.windows.showStatusBar = !initial;
            await waitUntil(() => fired, 'onPanelVisibilityChanged did not fire');
        } finally {
            subscription.unsubscribe();
            grok.shell.windows.showStatusBar = initial;
        }
    });

    test('onFileDragEnter', async () => {
        let fired = false;
        const tv = grok.shell.addTableView(grok.data.demo.demog(10));
        const subscription = grok.events.onFileDragEnter.subscribe((ev) => {
            fired = true;
            ev.preventDefault(); // suppresses the platform's drop overlay
        });
        const dt = new DataTransfer();
        dt.items.add(new File(['a,b\n1,2'], 'drag-enter-test.csv'));
        try {
            tv.root.dispatchEvent(new DragEvent('dragenter', {bubbles: true, cancelable: true, dataTransfer: dt}));
            await waitUntil(() => fired, 'onFileDragEnter did not fire');
        } finally {
            subscription.unsubscribe();
            grok.shell.closeTable(tv.dataFrame);
        }
    });

    test('debounce', async () => {
        let subscriptionPassed = false;
        let tv = grok.shell.addTableView(grok.data.demo.demog())
        let subscription = DG.debounce(tv.dataFrame.onColumnNameChanged, 500).subscribe(() => {
            subscriptionPassed = true;
        });

        try {
            tv.dataFrame.columns.toList()[0].name = "debounceTest"
            await waitUntil(() => subscriptionPassed, 'The debounced onColumnNameChanged did not fire');
        } finally {
            subscription.unsubscribe();
            tv.close();
        }
    });
});