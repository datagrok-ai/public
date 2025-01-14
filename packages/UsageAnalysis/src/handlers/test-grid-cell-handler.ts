import { category } from "@datagrok-libraries/utils/src/test";
import { delay, DockManager, GridCell } from "datagrok-api/dg";
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

export class TestGridCellHandler extends DG.ObjectHandler {
    get type(): string {
        return 'test';
    }

    isApplicable(x: any): boolean {
        return x instanceof DG.SemanticValue && x.semType === 'test';
    }

    renderProperties(semValue: DG.SemanticValue, context: any = null): HTMLElement {
        let panel = ui.accordion('testData');
        const testData = semValue.cell.value.split(':');

        const buttonsData = ui.button('Run', async () => {
            let df = semValue.cell.dataFrame;
            let progressBar = DG.TaskBarProgressIndicator.create(`Running test ${testData.join(':')}`);
            let category = testData[1];
            for (let pathPart of testData)
                if (testData[0] != pathPart && testData[1] != pathPart && testData[testData.length - 1] != pathPart)
                    category = [category, pathPart.trim()].join(':');
            let testResult: DG.DataFrame = await grok.functions.call(`${testData[0].trim()}:test`, { test: testData[testData.length - 1].trim(), category: category.trim() })
            progressBar.close();
            let wasError = false;
            let wasSkipped = false;
            for (let a of testResult.rows) {
                if (a.get('skipped'))
                    wasSkipped = true;
                else if (!a.get('success'))
                    wasError = true;
            }
            grok.shell.closeAll();
            let tv = grok.shell.addTableView(df);
            if (wasSkipped)
                grok.shell.warning(`${testData.join(':')} was skipped`);
            else if (!wasError)
                grok.shell.info(`${testData.join(':')} successfully run`);
            else
                grok.shell.error(`${testData.join(':')} failed`);

            await delay(500); 
            grok.shell.o = (semValue); 
        });
        buttonsData.classList.add('ui-btn-raised');
        buttonsData.classList.add('ui-btn-test-run');

        const packageDiv = ui.divH([ui.p('package:'), ui.h3(testData[0])]);
        packageDiv.classList.add('ui-test-data');
        panel.addPane('Run Test', () => ui.divV([ui.h1(`${testData[1]}: ${testData[2]}`), packageDiv, buttonsData]));
        return panel.root;
    }
}