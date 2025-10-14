import { category, delay, expect, expectArray, test } from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok'; 
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui'; 


category('Balloon', () => {
    test('closeAll', async () => {
        let balloon  = new DG.Balloon();
        balloon.info(ui.div());
        balloon.warning(ui.div());
        balloon.error(ui.div());
        expect(document.querySelectorAll('.d4-balloon').length >= 3, true);
        DG.Balloon.closeAll();
        expect(document.querySelectorAll('.d4-balloon').length, 0);
    });
});

