import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';


category('Demo: Projests', () => {
     test('project.city_gps', async () => {
        await grok.dapi.projects.open('city_gps')
     
      });
    
});