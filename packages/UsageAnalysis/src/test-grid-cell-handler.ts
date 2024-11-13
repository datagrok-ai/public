import { category } from "@datagrok-libraries/utils/src/test";
import { delay, DockManager, GridCell } from "datagrok-api/dg";

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
            await grok.functions.call(`${testData[0].trim()}:test`, { test: testData[2].trim(), category: testData[1].trim() })
            let tv = grok.shell.addTableView(df);

            await delay(100);
            tv.dataFrame.currentRowIdx =  gridCell.cell.rowIndex ?? 0;
            await delay(100);
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