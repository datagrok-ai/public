import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {parseTwbFile} from '../tableau/tableau-parser';


function makeMinimalTwb(datasourceXml: string = '', worksheetXml: string = ''): string {
  return `<?xml version='1.0' encoding='utf-8' ?>
<workbook version='18.1'>
  <datasources>${datasourceXml}</datasources>
  <worksheets>${worksheetXml}</worksheets>
</workbook>`;
}


category('Tableau: Parser', () => {
  test('single datasource with columns', async () => {
    const xml = makeMinimalTwb(`
      <datasource caption='demog' inline='true' name='ds1' version='18.1'>
        <connection class='federated'>
          <named-connections>
            <named-connection caption='demog' name='nc1'>
              <connection class='textscan' directory='/data' filename='demog.csv' />
            </named-connection>
          </named-connections>
          <relation connection='nc1' name='demog.csv' table='[demog#csv]' type='table'>
            <columns character-set='UTF-8' header='yes'>
              <column datatype='string' name='NAME' ordinal='0' />
              <column datatype='integer' name='AGE' ordinal='1' />
              <column datatype='real' name='WEIGHT' ordinal='2' />
            </columns>
          </relation>
        </connection>
        <column datatype='integer' name='[AGE]' role='measure' type='quantitative' />
        <column datatype='string' name='[NAME]' role='dimension' type='nominal' />
      </datasource>
    `);
    const result = parseTwbFile(xml);
    expect(result.datasources.length, 1);
    expect(result.datasources[0].caption, 'demog');
    expect(result.datasources[0].sourceFile, 'demog.csv');
    expect(result.datasources[0].columns.length, 3);

    const ageCol = result.datasources[0].columns.find((c) => c.name === 'AGE')!;
    expect(ageCol.datatype, 'integer');
    expect(ageCol.role, 'measure');
    expect(ageCol.type, 'quantitative');

    const nameCol = result.datasources[0].columns.find((c) => c.name === 'NAME')!;
    expect(nameCol.role, 'dimension');
    expect(nameCol.type, 'nominal');
  });

  test('multiple datasources', async () => {
    const xml = makeMinimalTwb(`
      <datasource caption='ds_a' name='a'>
        <connection class='federated'>
          <named-connections>
            <named-connection caption='a' name='nc1'>
              <connection class='textscan' directory='/data' filename='a.csv' />
            </named-connection>
          </named-connections>
          <relation connection='nc1' name='a.csv' table='[a#csv]' type='table'>
            <columns><column datatype='string' name='X' ordinal='0' /></columns>
          </relation>
        </connection>
      </datasource>
      <datasource caption='ds_b' name='b'>
        <connection class='federated'>
          <named-connections>
            <named-connection caption='b' name='nc2'>
              <connection class='textscan' directory='/data' filename='b.csv' />
            </named-connection>
          </named-connections>
          <relation connection='nc2' name='b.csv' table='[b#csv]' type='table'>
            <columns><column datatype='integer' name='Y' ordinal='0' /></columns>
          </relation>
        </connection>
      </datasource>
    `);
    const result = parseTwbFile(xml);
    expect(result.datasources.length, 2);
    expect(result.datasources[0].caption, 'ds_a');
    expect(result.datasources[1].caption, 'ds_b');
  });

  test('worksheets', async () => {
    const xml = makeMinimalTwb(
      `<datasource caption='data' name='ds1'>
        <connection class='federated'>
          <named-connections>
            <named-connection caption='data' name='nc1'>
              <connection class='textscan' directory='/data' filename='data.csv' />
            </named-connection>
          </named-connections>
          <relation connection='nc1' name='data.csv' table='[data#csv]' type='table'>
            <columns>
              <column datatype='string' name='A' ordinal='0' />
              <column datatype='integer' name='B' ordinal='1' />
            </columns>
          </relation>
        </connection>
      </datasource>`,
      `<worksheet name='Sheet 1'>
        <table>
          <view>
            <datasource-dependencies datasource='ds1'>
              <column name='[A]' datatype='string' role='dimension' type='nominal' />
              <column name='[B]' datatype='integer' role='measure' type='quantitative' />
            </datasource-dependencies>
          </view>
          <panes>
            <pane><mark class='Bar' /></pane>
          </panes>
          <rows>[ds1].[none:A:nk]</rows>
          <cols>[ds1].[sum:B:qk]</cols>
        </table>
      </worksheet>`
    );
    const result = parseTwbFile(xml);
    expect(result.worksheets.length, 1);
    expect(result.worksheets[0].name, 'Sheet 1');
    expect(result.worksheets[0].datasourceName, 'ds1');
    expect(result.worksheets[0].markClass, 'Bar');
    expect(result.worksheets[0].usedColumns.length, 2);
  });

  test('empty workbook', async () => {
    const xml = makeMinimalTwb();
    const result = parseTwbFile(xml);
    expect(result.datasources.length, 0);
    expect(result.worksheets.length, 0);
    expect(result.version, '18.1');
  });

  test('internal columns filtered out', async () => {
    const xml = makeMinimalTwb(`
      <datasource caption='test' name='ds1'>
        <connection class='federated'>
          <named-connections>
            <named-connection caption='test' name='nc1'>
              <connection class='textscan' directory='/data' filename='test.csv' />
            </named-connection>
          </named-connections>
          <relation connection='nc1' name='test.csv' table='[test#csv]' type='table'>
            <columns>
              <column datatype='string' name='COL1' ordinal='0' />
            </columns>
          </relation>
        </connection>
        <column datatype='table' name='[__tableau_internal_object_id__].[test_ABC]' role='measure' type='quantitative' />
      </datasource>
    `);
    const result = parseTwbFile(xml);
    expect(result.datasources[0].columns.length, 1);
    expect(result.datasources[0].columns[0].name, 'COL1');
  });

  test('metadata enrichment', async () => {
    const xml = makeMinimalTwb(`
      <datasource caption='test' name='ds1'>
        <connection class='federated'>
          <named-connections>
            <named-connection caption='test' name='nc1'>
              <connection class='textscan' directory='/data' filename='test.csv' />
            </named-connection>
          </named-connections>
          <relation connection='nc1' name='test.csv' table='[test#csv]' type='table'>
            <columns>
              <column datatype='integer' name='AGE' ordinal='0' />
              <column datatype='string' name='NAME' ordinal='1' />
            </columns>
          </relation>
          <metadata-records>
            <metadata-record class='column'>
              <remote-name>AGE</remote-name>
              <aggregation>Sum</aggregation>
              <contains-null>true</contains-null>
            </metadata-record>
            <metadata-record class='column'>
              <remote-name>NAME</remote-name>
              <aggregation>Count</aggregation>
              <contains-null>false</contains-null>
            </metadata-record>
          </metadata-records>
        </connection>
      </datasource>
    `);
    const result = parseTwbFile(xml);
    const ageCol = result.datasources[0].columns.find((c) => c.name === 'AGE')!;
    expect(ageCol.aggregation, 'Sum');
    expect(ageCol.containsNull, true);

    const nameCol = result.datasources[0].columns.find((c) => c.name === 'NAME')!;
    expect(nameCol.aggregation, 'Count');
    expect(nameCol.containsNull, false);
  });
});
