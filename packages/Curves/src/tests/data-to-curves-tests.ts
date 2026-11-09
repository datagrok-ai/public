import * as DG from 'datagrok-api/dg';

import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {convertDataToCurves} from '../fit/data-to-curves';

const WELLS = `compound_id,assay_name,target_entity,batch_id,run_id,concentration_uM,inhibition
C1,A,T,B1,R1,0.1,5
C1,A,T,B1,R1,1,20
C1,A,T,B1,R1,10,60
C1,A,T,B1,R1,100,95
C1,A,T,B1,R2,0.1,4
C1,A,T,B1,R2,1,25
C1,A,T,B1,R2,10,55
C1,A,T,B1,R2,100,90
C2,A,T,B2,R1,0.1,2
C2,A,T,B2,R1,1,10
C2,A,T,B2,R1,10,30
C2,A,T,B2,R1,100,50`;

category('Data to Curves', () => {
  test('wells group into one curve per compound, assay, target and run', async () => {
    const wells = DG.DataFrame.fromCsv(WELLS);
    const col = (name: string) => wells.col(name)!;
    const curves = await convertDataToCurves(wells, col('concentration_uM'), col('inhibition'), col('batch_id'),
      col('assay_name'), col('run_id'), col('compound_id'), col('target_entity'));
    const values = (name: string) => curves.col(name)!.toList().join(', ');
    expect(curves.rowCount, 3);
    expect(curves.col('Fitted Curve')!.semType, 'fit');
    expect(values('Compound|Assay|Target'), 'C1 | A | T, C1 | A | T, C2 | A | T');
    expect(values('Run ID'), 'R1, R2, R1');
    expect(values('Max Percent Inhibition'), '95, 90, 50');
    expect(values('Max Percent Inhibition (Compound | Assay | Target)'), '95, 95, 50');
  });
});
