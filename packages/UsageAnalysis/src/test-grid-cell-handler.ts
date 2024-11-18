import { category } from "@datagrok-libraries/utils/src/test";
import { delay, DockManager, GridCell } from "datagrok-api/dg";
import * as DG from 'datagrok-api/dg';

export class TestGridCellHandler extends DG.ObjectHandler {
    get type(): string {
        return 'test';
    }

    isApplicable(x: any): boolean {
        return x instanceof DG.SemanticValue && x.semType === 'test';
    }

    renderProperties(gridCell: GridCell, context: any = null): HTMLElement {
        let panel = ui.accordion('testData');
        const testData = gridCell.cell.value.split(':');

        const buttonsData = ui.button('Run', async () => {
            let df = gridCell.cell.dataFrame;
            let progressBar = DG.TaskBarProgressIndicator.create(`Running test ${testData.join(':')}`);
            let category = testData[1];
            for (let pathPart of testData)
                if (testData[0] != pathPart && testData[1] != pathPart && testData[testData.length - 1] != pathPart)
                    category = [category, pathPart.trim()].join(':');
            let testResult: DG.DataFrame = await grok.functions.call(`${testData[0].trim()}:test`, { test: testData[testData.length - 1].trim(), category: category.trim() })
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

            await delay(100);
            tv.dataFrame.currentRowIdx = gridCell.cell.rowIndex ?? 0;
            await delay(100);
            progressBar.close();
            df.currentCell = df.cell(gridCell.cell.rowIndex ?? 0, gridCell.cell.column.name);
        });
        buttonsData.classList.add('ui-btn-raised');
        buttonsData.classList.add('ui-btn-test-run');

        const packageDiv = ui.divH([ui.p('package:'), ui.h3(testData[0])]);
        packageDiv.classList.add('ui-test-data');
        panel.addPane('Run Test', () => ui.divV([ui.h1(`${testData[1]}: ${testData[2]}`), packageDiv, buttonsData]));
        return panel.root;
    }
}