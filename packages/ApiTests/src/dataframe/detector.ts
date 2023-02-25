/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/utils/src/test';


category('Detector', () => {
  test('SingleEmptyCategory', async () => {
    const df: DG.DataFrame = DG.DataFrame.fromObjects([
      {col1: null},
      {col1: null},
      {col1: null},
      {col1: null},
      {col1: null},
    ])!;
    const col: DG.Column = df.col('col1')!;
    const res: boolean = DG.Detector.sampleCategories(col, (v) => false, 1);
    if (res)
      throw new Error('DG.Detector.sampleCategories() always returns true on single empty category.');
  });
});

const testData = DG.DataFrame.fromCsv(`countries,fasta,smiles,molregno,LON,Zip Code,Street Address Line 1,ImageUrl,user_id,error_message,xray,flag,magnitude,CS-id,pdb_id,accel_a,time_offset,empty_number,empty_string
Belgium,MSNFHNEHVMQFYRNNLKTKGVFGRQ,CC(C(=O)OCCCc1cccnc1)c2cccc(c2)C(=O)c3ccccc3,1480014,36.276729583740234,995042300,14016 ROUTE 31W,https://datagrok.ai/img/slides/access-db-connect.png,id,ErrorMessage,COMPND ATOM END,flag,1,1,1QBS,1,1.23,100,abc
Burundi,MDYKETLLMPKTDFPMRGGLPNKEPQIQEKW,COc1ccc2cc(ccc2c1)C(C)C(=O)Oc3ccc(C)cc3OC,1480015,36.276729583740234,995073444,80 STATE HIGHWAY 310,https://datagrok.ai/img/slides/access-db-connect.png,id,ErrorMessage,COMPND ATOM END,flag,2,2,1ZP8,2,1.23,,
Cameroon,MIEVFLFGIVLGLIPITLAGLFVTAYLQYRRGDQLDL,COc1ccc2cc(ccc2c1)C(C)C(=O)OCCCc3cccnc3,1480016,36.26095962524414,995153596,30-56 WHITESTONE EXPY,https://datagrok.ai/img/slides/access-db-connect.png,id,ErrorMessage,COMPND ATOM END,flag,3,3,2BDJ,3,1.23,,
Canada,MMELVLKTIIGPIVVGVVLRIVDKWLNKDK,CC(C(=O)NCCS)c1cccc(c1)C(=O)c2ccccc2,1480017,36.26095962524414,99515,30-56 WHITESTONE EXPY,https://datagrok.ai/img/slides/access-db-connect.png,id,ErrorMessage,COMPND ATOM END,flag,4,4,1IAN,4,1.23,,
Colombia,MDRTDEVSNHTHDKPTLTWFEEIFEEYHSPFHN,FC(F)(F)c1ccc(OC2CCNCC2)cc1,1480029,36.3309440612793,995152050,1 COURT HOUSE SQUARE,https://datagrok.ai/img/slides/access-db-connect.png,id,ErrorMessage,COMPND ATOM END,flag,5,5,4UJ1,5,1.23,,
Costa Rica,MKSTKEEIQTIKTLLKDSRTAKYHKRLQIVL,CC(C)Cc1ccc(cc1)C(C)C(=O)N2CCCC2C(=O)OCCCc3ccccc3,1480018,36.3309440612793,995084218,4041 SOUTHWESTERN BLVD,https://datagrok.ai/img/slides/access-db-connect.png,id,ErrorMessage,COMPND ATOM END,flag,6,6,2BPW,6,1.23,,
Cuba,MHAILRYFIRRLFYHIFYKIYSLISKKHQSLPSDVRQF,COc1ccc2c(c1)c(CC(=O)N3CCCC3C(=O)Oc4ccc(C)cc4OC)c(C)n2C(=O)c5ccc(Cl)cc5,1480019,36.33115768432617,995081928,1227 US HIGHWAY 11,https://datagrok.ai/img/slides/access-db-connect.png,id,ErrorMessage,COMPND ATOM END,flag,7,7,1QBS,7,1.23,,
Italy,MSNFHNEHVMQFYRNNLKTKGVFGRQ,CC(C)Cc1ccc(cc1)C(C)C(=O)N2CCCC2C(=O)OCCO[N+](=O)[O-],1480020,36.33115768432617,99502,"168-46 91ST AVE., 2ND FLR",https://datagrok.ai/img/slides/access-db-connect.png,id,ErrorMessage,COMPND ATOM END,flag,8,8,1ZP8,8,1.23,,
Rwanda,MPNSEPASLLELFNSIATQGELVRSLKAGNASK,CC(C)Cc1ccc(cc1)C(C)C(=O)N2CCCC2C(=O)OCCO,1480021,36.33137130737305,995037247,"168-46 91ST AVE., 2ND FLR",https://datagrok.ai/img/slides/access-db-connect.png,id,ErrorMessage,COMPND ATOM END,flag,9,9,2BDJ,9,1.23,,
Switzerland,IRVVGRYLIEVWKAAGMDMDKVLFLWSSDEI,CN1CCC(CC1)Oc2ccc(cc2)C(F)(F)F,1480028,36.33137130737305,99504,92-11 179TH PLACE,https://datagrok.ai/img/slides/access-db-connect.png,id,ErrorMessage,COMPND ATOM END,flag,9,10,1IAN,10,1.23,,
,,,,,,,,,,,,,,,,,,`);
testData.columns.add(DG.Column.fromList(DG.TYPE.BYTE_ARRAY, 'BinaryImage', Array.from(new Uint8Array(11))));

category('Detector: All Detectors', () => {
  const detectors = DG.Func.find({tags: ['semTypeDetector']});
  for (const detector of detectors) {
    test(detector.friendlyName, async () => {
      const arr = [];
      const cols = testData.clone().columns;
      for (const col of cols) {
        const res = await detector.apply([col]);
        arr.push(res || col.semType);
      }
      expect(arr.filter((i) => i).length, 1);
    });
  }
});
