import { after, before, category, delay, expect, test } from "@datagrok-libraries/utils/src/test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { checkHTMLElement } from './utils';

category('UI: Grid', () => {
    let v: DG.TableView;

    before(async () => {
        let demog = grok.data.demo.demog(1000);
        v = grok.shell.addTableView(demog);
    });

    test('grid.setOrder', async () => {
       v.grid.columns.setOrder(['race', 'age']);
       let firstCol = v.grid.columns.byIndex(1);
       let secondCol = v.grid.columns.byIndex(2);
       if(firstCol?.dart.columnName != 'race' || secondCol?.dart.columnName != 'age'){
           throw 'grid.setOrder does not work'
       }
    });

    test('grid.resizeColumn', async () => {
        let check = false;

        v.grid.onColumnResized.subscribe((_) => { check = true });
        v.grid.columns.byName('age')!.width = 200;

        if (check == false)
            throw 'Column Resize error';
    });

    /*after(async () => {
        v.close();
        grok.shell.closeAll();
    }); */

});