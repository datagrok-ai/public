/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {FitConstants} from '@datagrok-libraries/statistics/src/fit/const';
import {IFitChartData} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {convertXmlCurveToJson} from '../fit/converters/xml-converter';
import {convertCompactDrToJson} from '../fit/converters/compact-dr-converter';
import {convertPzfxToJson} from '../fit/converters/pzfx-converter';
import {parseCellValue, registerCurveConverter, isNativeFormat, hasConverter} from '../fit/curve-converter';


// Sample data for each format
const SAMPLE_XML = '<chart><grid><xAxis min="0.0" max="5" /><yAxis min="-100" max="100" /></grid><settings xLabel="Concentration" yLabel="Activity" logX="true" /><seriesCollection><series name="Series1"><x>7.2e-12,7.1e-12,7e-12,6e-8,6e-5,6e-4</x><y>112,115,116,100,123,124</y><mask>111101</mask><fitLine><function type="sigif"><params>2e-08,-1.5,125,-20</params></function></fitLine><settings color="Green" marker="SquareBorder" markerColor="Yellow" markerBorderColor="Black" errorBarColor="Red" showErrorBars="True" drawLine="True" bound="True" /></series></seriesCollection></chart>';

const SAMPLE_COMPACT_DR = '{"rdt":"2025-05-31 00:00:00.0","mnr":-4.28,"mxr":109.134,"yu":"%","hs":-1.5498,"poi":0.0064758,"xu":"uM","ctp":"dose-response","p":[[0.001,3],[0.003,4],[0.01,7],[0.03,15],[0.1,35],[0.3,65],[1,85],[3,98],[10,105]]}';

const SAMPLE_COMPACT_DR_WITH_OUTLIER = '{"rdt":"2025-05-31 00:00:00.0","mnr":-4.28,"mxr":109.134,"yu":"%","hs":-1.5498,"poi":0.0064758,"xu":"uM","ctp":"dose-response","p":[[0.001,3],[0.003,4],[0.01,7,1],[0.03,15],[0.1,35]]}';

const SAMPLE_NATIVE_JSON = '{"series":[{"fitFunction":"sigmoid","parameters":[100,-1.5,0.01,-5],"fitLineColor":"#1f77b4","pointColor":"#1f77b4","showFitLine":true,"showPoints":"points","points":[{"x":0.001,"y":3},{"x":0.01,"y":7},{"x":0.1,"y":35},{"x":1,"y":85},{"x":10,"y":105}]}],"chartOptions":{"logX":true,"xAxisName":"Concentration","yAxisName":"Response"}}';

const SAMPLE_PZFX = `<?xml version="1.0" encoding="UTF-8"?>
<GraphPadPrismFile xmlns="http://graphpad.com/prism/Prism.htm" PrismXMLVersion="5.00">
  <Table ID="Table0" XFormat="numbers" TableType="XY">
    <Title>IC50 Curve</Title>
    <XColumn Width="81" Decimals="2" Subcolumns="1">
      <Title>Dose</Title>
      <Subcolumn>
        <d>1e-9</d>
        <d>1e-8</d>
        <d>1e-7</d>
        <d>1e-6</d>
        <d>1e-5</d>
      </Subcolumn>
    </XColumn>
    <YColumn Width="81" Decimals="2" Subcolumns="1">
      <Title>Response A</Title>
      <Subcolumn>
        <d>5.1</d>
        <d>25.3</d>
        <d>52.7</d>
        <d>85.4</d>
        <d>97.2</d>
      </Subcolumn>
    </YColumn>
  </Table>
</GraphPadPrismFile>`;

const SAMPLE_PZFX_WITH_EXCLUSION = `<?xml version="1.0" encoding="UTF-8"?>
<GraphPadPrismFile xmlns="http://graphpad.com/prism/Prism.htm" PrismXMLVersion="5.00">
  <Table ID="Table0" XFormat="numbers" TableType="XY">
    <Title>With Outlier</Title>
    <XColumn Width="81" Decimals="0" Subcolumns="1">
      <Title>X</Title>
      <Subcolumn><d>1</d><d>2</d><d>3</d></Subcolumn>
    </XColumn>
    <YColumn Width="81" Decimals="0" Subcolumns="1">
      <Title>Y</Title>
      <Subcolumn>
        <d>10</d>
        <d Excluded="1">99</d>
        <d>30</d>
      </Subcolumn>
    </YColumn>
  </Table>
</GraphPadPrismFile>`;


category('converters', () => {
  test('xmlConverter.basicConversion', async () => {
    const result = convertXmlCurveToJson(SAMPLE_XML);
    const chartData: IFitChartData = JSON.parse(result);
    expect(chartData.series !== undefined && chartData.series.length > 0, true);
    expect(chartData.series![0].points.length, 6);
    expect(chartData.chartOptions?.xAxisName, 'Concentration');
    expect(chartData.chartOptions?.yAxisName, 'Activity');
    expect(chartData.chartOptions?.logX, true);
  });

  test('xmlConverter.parameterReordering', async () => {
    const result = convertXmlCurveToJson(SAMPLE_XML);
    const chartData: IFitChartData = JSON.parse(result);
    // XML params: [IC50, tan, max, min] = [2e-08, -1.5, 125, -20]
    // Native params: [max, tan, IC50, min] = [125, -1.5, 2e-08, -20]
    const params = chartData.series![0].parameters!;
    expect(params[0], 125);
    expect(params[1], -1.5);
    expect(params[2], 2e-08);
    expect(params[3], -20);
  });

  test('xmlConverter.outlierMask', async () => {
    const result = convertXmlCurveToJson(SAMPLE_XML);
    const chartData: IFitChartData = JSON.parse(result);
    // mask "111101": '1' = included (not outlier), '0' = excluded (outlier)
    expect(chartData.series![0].points.length, 6);
    expect(chartData.series![0].points[0].outlier, false);
    expect(chartData.series![0].points[3].outlier, false);
    expect(chartData.series![0].points[4].outlier, true);
    expect(chartData.series![0].points[5].outlier, false);
  });

  test('compactDrConverter.basicConversion', async () => {
    const result = convertCompactDrToJson(SAMPLE_COMPACT_DR);
    const chartData: IFitChartData = JSON.parse(result);
    expect(chartData.series !== undefined && chartData.series.length === 1, true);
    expect(chartData.series![0].points.length, 9);
    expect(chartData.chartOptions?.xAxisName, 'uM');
    expect(chartData.chartOptions?.yAxisName, '%');
    expect(chartData.chartOptions?.logX, true);
  });

  test('compactDrConverter.parameters', async () => {
    const result = convertCompactDrToJson(SAMPLE_COMPACT_DR);
    const chartData: IFitChartData = JSON.parse(result);
    // Parameters should be [mxr, hs, poi, mnr] = [max, slope, IC50, min]
    const params = chartData.series![0].parameters!;
    expect(params[0], 109.134);
    expect(params[1], -1.5498);
    expect(params[2], 0.0064758);
    expect(params[3], -4.28);
  });

  test('compactDrConverter.outlierFlag', async () => {
    const result = convertCompactDrToJson(SAMPLE_COMPACT_DR_WITH_OUTLIER);
    const chartData: IFitChartData = JSON.parse(result);
    expect(chartData.series![0].points[0].outlier, false);
    expect(chartData.series![0].points[1].outlier, false);
    expect(chartData.series![0].points[2].outlier, true);
    expect(chartData.series![0].points[3].outlier, false);
  });

  test('compactDrConverter.fitFunction', async () => {
    const result = convertCompactDrToJson(SAMPLE_COMPACT_DR);
    const chartData: IFitChartData = JSON.parse(result);
    expect(chartData.series![0].fitFunction, '4pl-dose-response');
    expect(chartData.series![0].showFitLine, true);
    expect(chartData.series![0].showPoints, 'points');
  });

  test('parseCellValue.nativeJson', async () => {
    const chartData = parseCellValue(SAMPLE_NATIVE_JSON, null);
    expect(chartData.series !== undefined && chartData.series.length === 1, true);
    expect(chartData.series![0].points.length, 5);
    expect(chartData.chartOptions?.logX, true);
  });

  test('parseCellValue.withConverter', async () => {
    // Register the converter
    registerCurveConverter('compact-dr', convertCompactDrToJson);

    // Create a column with the converter tag
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('curves', [SAMPLE_COMPACT_DR])]);
    const col = df.col('curves')!;
    col.setTag(FitConstants.TAG_CURVE_FORMAT, 'compact-dr');

    const chartData = parseCellValue(SAMPLE_COMPACT_DR, col);
    expect(chartData.series !== undefined && chartData.series.length === 1, true);
    expect(chartData.series![0].points.length, 9);
    expect(chartData.chartOptions?.logX, true);
  });

  test('isNativeFormat.nativeColumn', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('curves', [SAMPLE_NATIVE_JSON])]);
    const col = df.col('curves')!;
    expect(isNativeFormat(col), true);
    expect(hasConverter(col), false);
  });

  test('isNativeFormat.converterColumn', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('curves', [SAMPLE_COMPACT_DR])]);
    const col = df.col('curves')!;
    col.setTag(FitConstants.TAG_CURVE_FORMAT, 'compact-dr');
    expect(isNativeFormat(col), false);
    expect(hasConverter(col), true);
  });

  test('parseCellValue.emptyValue', async () => {
    const chartData = parseCellValue('', null);
    expect(chartData.series !== undefined && chartData.series.length === 0, true);
  });

  test('converterCaching', async () => {
    registerCurveConverter('compact-dr', convertCompactDrToJson);

    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('curves', [SAMPLE_COMPACT_DR])]);
    const col = df.col('curves')!;
    col.setTag(FitConstants.TAG_CURVE_FORMAT, 'compact-dr');

    // Call twice with the same value - second call should use cache
    const result1 = parseCellValue(SAMPLE_COMPACT_DR, col);
    const result2 = parseCellValue(SAMPLE_COMPACT_DR, col);
    expect(result1.series![0].points.length, result2.series![0].points.length);
  });
});


category('converters: rendering integration', () => {
  test('xmlColumnRendering', async () => {
    registerCurveConverter('3dx', convertXmlCurveToJson);

    const testCsv = `xmlCurve
"<chart><grid><xAxis min=""0.0"" max=""5"" /><yAxis min=""-100"" max=""100"" /></grid><settings xLabel=""Concentration"" yLabel=""Activity"" logX=""true"" /><seriesCollection><series name=""S1""><x>7.2e-12,7.1e-12,7e-12,6e-8,6e-5,6e-4</x><y>112,115,116,100,123,124</y><mask>111111</mask><fitLine><function type=""sigif""><params>2e-08,-1.5,125,-20</params></function></fitLine><settings color=""Green"" marker=""SquareBorder"" markerColor=""Yellow"" markerBorderColor=""Black"" errorBarColor=""Red"" showErrorBars=""True"" drawLine=""True"" bound=""True"" /></series></seriesCollection></chart>"`;
    const df = DG.DataFrame.fromCsv(testCsv);
    const col = df.columns.byName('xmlCurve');
    await grok.data.detectSemanticTypes(df);

    expect(col.semType, FitConstants.FIT_SEM_TYPE);

    // Verify that parseCellValue works through the converter
    const chartData = parseCellValue(col.get(0), col);
    expect(chartData.series !== undefined && chartData.series.length > 0, true);
    expect(chartData.series![0].points.length, 6);
  });

  test('compactDrColumnRendering', async () => {
    registerCurveConverter('compact-dr', convertCompactDrToJson);

    const testCsv = `compactDR
"{""rdt"":""2025-05-31 00:00:00.0"",""mnr"":-4.28,""mxr"":109.134,""yu"":""%"",""hs"":-1.5498,""poi"":0.0064758,""xu"":""uM"",""ctp"":""dose-response"",""p"":[[0.001,3],[0.003,4],[0.01,7],[0.03,15],[0.1,35],[0.3,65],[1,85],[3,98],[10,105]]}"`;
    const df = DG.DataFrame.fromCsv(testCsv);
    const col = df.columns.byName('compactDR');
    await grok.data.detectSemanticTypes(df);

    expect(col.semType, FitConstants.FIT_SEM_TYPE);

    const chartData = parseCellValue(col.get(0), col);
    expect(chartData.series !== undefined && chartData.series.length === 1, true);
    expect(chartData.series![0].parameters![0], 109.134);
  });

  test('nativeJsonColumnRendering', async () => {
    const nativeJson = '{"series":[{"fitFunction":"sigmoid","parameters":[100,-1.5,0.01,-5],"fitLineColor":"#1f77b4","pointColor":"#1f77b4","showFitLine":true,"showPoints":"points","points":[{"x":0.001,"y":3},{"x":0.01,"y":7},{"x":0.1,"y":35},{"x":1,"y":85},{"x":10,"y":105}]}],"chartOptions":{"logX":true}}';
    const testCsv = `nativeJSON\n"${nativeJson.replace(/"/g, '""')}"`;
    const df = DG.DataFrame.fromCsv(testCsv);
    const col = df.columns.byName('nativeJSON');
    await grok.data.detectSemanticTypes(df);

    expect(col.semType, FitConstants.FIT_SEM_TYPE);
    // No converter tag for native format
    expect(isNativeFormat(col), true);

    const chartData = parseCellValue(col.get(0), col);
    expect(chartData.series![0].points.length, 5);
  });

  test('pzfxColumnRendering', async () => {
    registerCurveConverter('pzfx', convertPzfxToJson);

    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('pzfx', [SAMPLE_PZFX])]);
    const col = df.col('pzfx')!;
    col.setTag(FitConstants.TAG_CURVE_FORMAT, 'pzfx');
    col.semType = FitConstants.FIT_SEM_TYPE;

    const chartData = parseCellValue(col.get(0), col);
    expect(chartData.series !== undefined && chartData.series.length === 1, true);
    expect(chartData.series![0].points.length, 5);
    expect(chartData.chartOptions?.logX, true);
  });
});


category('converters: pzfx', () => {
  test('pzfxConverter.basicConversion', async () => {
    const result = convertPzfxToJson(SAMPLE_PZFX);
    const chartData: IFitChartData = JSON.parse(result);
    expect(chartData.series !== undefined && chartData.series.length === 1, true);
    expect(chartData.series![0].points.length, 5);
    expect(chartData.series![0].name, 'Response A');
    expect(chartData.chartOptions?.xAxisName, 'Dose');
    expect(chartData.chartOptions?.logX, true);
  });

  test('pzfxConverter.pointValues', async () => {
    const result = convertPzfxToJson(SAMPLE_PZFX);
    const chartData: IFitChartData = JSON.parse(result);
    const points = chartData.series![0].points;
    expect(points[0].x, 1e-9);
    expect(points[0].y, 5.1);
    expect(points[4].x, 1e-5);
    expect(points[4].y, 97.2);
  });

  test('pzfxConverter.outlierExclusion', async () => {
    const result = convertPzfxToJson(SAMPLE_PZFX_WITH_EXCLUSION);
    const chartData: IFitChartData = JSON.parse(result);
    expect(chartData.series![0].points.length, 3);
    expect(chartData.series![0].points[0].outlier, false);
    expect(chartData.series![0].points[1].outlier, true);
    expect(chartData.series![0].points[2].outlier, false);
  });

  test('pzfxConverter.withConverterTag', async () => {
    registerCurveConverter('pzfx', convertPzfxToJson);

    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('curves', [SAMPLE_PZFX])]);
    const col = df.col('curves')!;
    col.setTag(FitConstants.TAG_CURVE_FORMAT, 'pzfx');

    const chartData = parseCellValue(SAMPLE_PZFX, col);
    expect(chartData.series !== undefined && chartData.series.length === 1, true);
    expect(chartData.series![0].points.length, 5);
  });
});
