import { category } from "@datagrok-libraries/utils/src/test";
import { GridCell } from "datagrok-api/dg";

export class TestGridCellHandler extends DG.ObjectHandler {
    get type(): string {
        return 'test';
    }

    isApplicable(x: any): boolean {
        return x instanceof DG.SemanticValue && x.semType === 'test';
    }

    renderProperties(gridCell: any, context: any = null): HTMLElement {
        let panel = ui.accordion('testData'); 
        const testData = gridCell.cell.value.split(':');
        
        const buttonsData =   ui.button('Run', async ()=>{await grok.functions.call(`${testData[0].trim()}:test`, {test:testData[2].trim(), category:testData[1].trim()})});
        buttonsData.classList.add('ui-btn-raised');
        panel.addPane('Run Test', ()=>ui.divV([ui.h3(gridCell.cell.value), buttonsData]));
        return panel.root;
    }
}