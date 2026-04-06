import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {parsePzfxXml, pzfxTableToFitChartData, pzfxToFitDataFrame, pzfxTableToDataFrame} from '../formats/pzfx/pzfx-parser';


const DOSE_RESPONSE_PZFX = `<?xml version="1.0" encoding="UTF-8"?>
<GraphPadPrismFile xmlns="http://graphpad.com/prism/Prism.htm" PrismXMLVersion="5.00">
  <Table ID="Table0" XFormat="numbers" TableType="XY">
    <Title>Test Dose Response</Title>
    <XColumn Width="81" Decimals="2" Subcolumns="1">
      <Title>Concentration</Title>
      <Subcolumn>
        <d>1e-9</d>
        <d>1e-8</d>
        <d>1e-7</d>
        <d>1e-6</d>
        <d>1e-5</d>
      </Subcolumn>
    </XColumn>
    <YColumn Width="81" Decimals="2" Subcolumns="2">
      <Title>Compound X</Title>
      <Subcolumn>
        <d>5.1</d>
        <d>25.3</d>
        <d>52.7</d>
        <d>85.4</d>
        <d>97.2</d>
      </Subcolumn>
      <Subcolumn>
        <d>4.8</d>
        <d>23.1</d>
        <d>55.2</d>
        <d>83.9</d>
        <d>96.8</d>
      </Subcolumn>
    </YColumn>
    <YColumn Width="81" Decimals="2" Subcolumns="1">
      <Title>Compound Y</Title>
      <Subcolumn>
        <d>2.1</d>
        <d>8.5</d>
        <d>30.2</d>
        <d>65.7</d>
        <d>90.1</d>
      </Subcolumn>
    </YColumn>
  </Table>
</GraphPadPrismFile>`;

const MIXED_PZFX = `<?xml version="1.0" encoding="UTF-8"?>
<GraphPadPrismFile xmlns="http://graphpad.com/prism/Prism.htm" PrismXMLVersion="5.00">
  <Table ID="Table0" XFormat="numbers" TableType="XY">
    <Title>XY Data</Title>
    <XColumn Width="81" Decimals="0" Subcolumns="1">
      <Title>X</Title>
      <Subcolumn><d>1</d><d>2</d><d>3</d></Subcolumn>
    </XColumn>
    <YColumn Width="81" Decimals="0" Subcolumns="1">
      <Title>Y</Title>
      <Subcolumn><d>10</d><d>20</d><d>30</d></Subcolumn>
    </YColumn>
  </Table>
  <Table ID="Table1" XFormat="none" TableType="OneWay">
    <Title>Column Data</Title>
    <YColumn Width="81" Decimals="0" Subcolumns="1">
      <Title>Group A</Title>
      <Subcolumn><d>5</d><d>6</d><d>7</d></Subcolumn>
    </YColumn>
  </Table>
</GraphPadPrismFile>`;

const EXCLUDED_PZFX = `<?xml version="1.0" encoding="UTF-8"?>
<GraphPadPrismFile xmlns="http://graphpad.com/prism/Prism.htm" PrismXMLVersion="5.00">
  <Table ID="Table0" XFormat="numbers" TableType="XY">
    <Title>With Exclusions</Title>
    <XColumn Width="81" Decimals="0" Subcolumns="1">
      <Title>X</Title>
      <Subcolumn><d>1</d><d>2</d><d>3</d><d>4</d></Subcolumn>
    </XColumn>
    <YColumn Width="81" Decimals="0" Subcolumns="1">
      <Title>Y</Title>
      <Subcolumn>
        <d>10</d>
        <d Excluded="1">99</d>
        <d>30</d>
        <d></d>
      </Subcolumn>
    </YColumn>
  </Table>
</GraphPadPrismFile>`;


category('PZFX Parser', () => {
  test('parsePzfxXml: parses XY table', async () => {
    const tables = parsePzfxXml(DOSE_RESPONSE_PZFX);
    expect(tables.length, 1);
    expect(tables[0].tableType, 'XY');
    expect(tables[0].title, 'Test Dose Response');
    expect(tables[0].xValues.length, 5);
    expect(tables[0].yColumns.length, 2);
    expect(tables[0].yColumns[0].title, 'Compound X');
    expect(tables[0].yColumns[0].subcolumns.length, 2);
    expect(tables[0].yColumns[1].title, 'Compound Y');
    expect(tables[0].yColumns[1].subcolumns.length, 1);
  });

  test('parsePzfxXml: parses multiple tables', async () => {
    const tables = parsePzfxXml(MIXED_PZFX);
    expect(tables.length, 2);
    expect(tables[0].tableType, 'XY');
    expect(tables[1].tableType, 'OneWay');
    expect(tables[1].title, 'Column Data');
  });

  test('parsePzfxXml: handles excluded and missing data', async () => {
    const tables = parsePzfxXml(EXCLUDED_PZFX);
    const subcol = tables[0].yColumns[0].subcolumns[0];
    expect(subcol.values[0], 10);
    expect(subcol.excluded[1], true);
    expect(subcol.values[3], null);
    expect(subcol.excluded[0], false);
  });

  test('parsePzfxXml: handles trailing binary data', async () => {
    const withTrailing = DOSE_RESPONSE_PZFX + '\x00\x01\x02binary garbage';
    const tables = parsePzfxXml(withTrailing);
    expect(tables.length, 1);
  });

  test('parsePzfxXml: returns empty for malformed XML', async () => {
    const tables = parsePzfxXml('<not valid xml');
    expect(tables.length, 0);
  });

  test('pzfxTableToFitChartData: converts XY table', async () => {
    const tables = parsePzfxXml(DOSE_RESPONSE_PZFX);
    const chartData = pzfxTableToFitChartData(tables[0]);
    expect(chartData !== null, true);
    expect(chartData!.series!.length, 2);
    expect(chartData!.series![0].name, 'Compound X');
    // 5 x-values * 2 replicates = 10 points
    expect(chartData!.series![0].points.length, 10);
    // 5 x-values * 1 replicate = 5 points
    expect(chartData!.series![1].points.length, 5);
    expect(chartData!.chartOptions!.logX, true);
  });

  test('pzfxTableToFitChartData: returns null for non-XY', async () => {
    const tables = parsePzfxXml(MIXED_PZFX);
    const chartData = pzfxTableToFitChartData(tables[1]);
    expect(chartData, null);
  });

  test('pzfxTableToFitChartData: maps outliers', async () => {
    const tables = parsePzfxXml(EXCLUDED_PZFX);
    const chartData = pzfxTableToFitChartData(tables[0]);
    expect(chartData !== null, true);
    // 4 x-values, 1 replicate, but one is null → 3 points
    expect(chartData!.series![0].points.length, 3);
    // second point (x=2) should be an outlier
    const outlierPoint = chartData!.series![0].points.find((p) => p.x === 2);
    expect(outlierPoint!.outlier, true);
  });

  test('pzfxToFitDataFrame: creates fit DataFrame', async () => {
    const tables = parsePzfxXml(DOSE_RESPONSE_PZFX);
    const df = pzfxToFitDataFrame(tables);
    expect(df.rowCount, 1);
    expect(df.col('Fitted Curve') !== null, true);
    expect(df.col('Fitted Curve')!.semType, 'fit');
  });

  test('pzfxTableToDataFrame: creates plain DataFrame', async () => {
    const tables = parsePzfxXml(MIXED_PZFX);
    const df = pzfxTableToDataFrame(tables[1]);
    expect(df.rowCount, 3);
    expect(df.col('Group A') !== null, true);
  });
});
