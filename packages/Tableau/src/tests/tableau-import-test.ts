import {category, test, expect} from '@datagrok-libraries/test/src/test';
import * as DG from 'datagrok-api/dg';

import {parseTwbFile} from '../tableau/tableau-parser';
import {twbDatasourceToDataFrame} from '../tableau/tableau-to-dataframe';


function makeTestTwb(): string {
  return `<?xml version='1.0' encoding='utf-8' ?>
<workbook version='18.1'>
  <datasources>
    <datasource caption='demog' inline='true' name='ds1' version='18.1'>
      <connection class='federated'>
        <named-connections>
          <named-connection caption='demog' name='nc1'>
            <connection class='textscan' directory='/data' filename='demog.csv' />
          </named-connection>
        </named-connections>
        <relation connection='nc1' name='demog.csv' table='[demog#csv]' type='table'>
          <columns character-set='UTF-8' header='yes'>
            <column datatype='string' name='USUBJID' ordinal='0' />
            <column datatype='integer' name='AGE' ordinal='1' />
            <column datatype='real' name='WEIGHT' ordinal='2' />
            <column datatype='boolean' name='CONTROL' ordinal='3' />
            <column datatype='date' name='STARTED' ordinal='4' />
          </columns>
        </relation>
        <metadata-records>
          <metadata-record class='column'>
            <remote-name>AGE</remote-name>
            <aggregation>Sum</aggregation>
            <contains-null>true</contains-null>
          </metadata-record>
        </metadata-records>
      </connection>
      <column caption='Age' datatype='integer' name='[AGE]' role='measure' type='quantitative' />
      <column caption='Usubjid' datatype='string' name='[USUBJID]' role='dimension' type='nominal' />
    </datasource>
  </datasources>
  <worksheets></worksheets>
</workbook>`;
}


category('Tableau: Import', () => {
  test('DataFrame has datasource columns', async () => {
    const twb = parseTwbFile(makeTestTwb());
    const df = twbDatasourceToDataFrame(twb.datasources[0]);
    expect(df.name, 'demog');
    expect(df.columns.length, 5);
    expect(df.rowCount, 0);
    expect(df.columns.byName('USUBJID') !== null, true);
    expect(df.columns.byName('AGE') !== null, true);
    expect(df.columns.byName('WEIGHT') !== null, true);
    expect(df.columns.byName('CONTROL') !== null, true);
    expect(df.columns.byName('STARTED') !== null, true);
  });

  test('column types mapped correctly', async () => {
    const twb = parseTwbFile(makeTestTwb());
    const df = twbDatasourceToDataFrame(twb.datasources[0]);
    expect(df.columns.byName('USUBJID').type, DG.COLUMN_TYPE.STRING);
    expect(df.columns.byName('AGE').type, DG.COLUMN_TYPE.INT);
    expect(df.columns.byName('WEIGHT').type, DG.COLUMN_TYPE.FLOAT);
    expect(df.columns.byName('CONTROL').type, DG.COLUMN_TYPE.BOOL);
    expect(df.columns.byName('STARTED').type, DG.COLUMN_TYPE.DATE_TIME);
  });

  test('column tags from metadata', async () => {
    const twb = parseTwbFile(makeTestTwb());
    const df = twbDatasourceToDataFrame(twb.datasources[0]);

    const ageCol = df.columns.byName('AGE');
    expect(ageCol.getTag('tableau.role'), 'measure');
    expect(ageCol.getTag('tableau.type'), 'quantitative');
    expect(ageCol.getTag('tableau.aggregation'), 'Sum');
    expect(ageCol.getTag('tableau.containsNull'), 'true');
    expect(ageCol.getTag('tableau.caption'), 'Age');
    expect(ageCol.getTag('description'), 'Age');

    const idCol = df.columns.byName('USUBJID');
    expect(idCol.getTag('tableau.role'), 'dimension');
    expect(idCol.getTag('tableau.type'), 'nominal');
    expect(idCol.getTag('tableau.caption'), 'Usubjid');
  });

  test('table-level tags', async () => {
    const twb = parseTwbFile(makeTestTwb());
    const df = twbDatasourceToDataFrame(twb.datasources[0]);
    expect(df.getTag('source.format'), 'tableau-twb');
    expect(df.getTag('tableau.sourceFile'), 'demog.csv');
    expect(df.getTag('tableau.sourceDirectory'), '/data');
  });
});
