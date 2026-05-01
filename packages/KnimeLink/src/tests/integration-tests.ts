import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {KnimeDeployment} from '../types';
import {PackageFunctions} from '../package';
import {_package} from '../package-test';
import { funcs } from '../package-api';

const TEST_DEPLOYMENT: KnimeDeployment = {
  id: 'rest:aadd942d-1528-4e1c-9cfa-20434ea7c677',
  name: 'test_inputs',
  type: 'rest',
};

const EXPECTED_INPUT_NAMES = [
  'input_file',
  'json_input_my_int_value',
  'json_input_my_double_value',
  'json_input_my_string_value',
  'table_input',
  'variable_input_race',
  'variable_input_weight',
  'arbitary_variable_input',
];

const EXPECTED_OUTPUT_NAMES = ['json_output', 'table_output', 'file_output'];


category('KnimeLink: integration', () => {

  test('lists deployments and finds test_inputs', async () => {
    const deployments: KnimeDeployment[] = await funcs.knimeListDeployments('rest');
    expect(deployments.length > 0, true);
    const dep = deployments.find((d) => d.id === TEST_DEPLOYMENT.id);
    expect(dep !== undefined, true);
    expect(dep!.name.toLowerCase(), 'test_inputs');
  });

  test('registers function with correct inputs and outputs', async () => {
    const func: DG.Func = await funcs.knimeGetOrRegisterFunc(
      TEST_DEPLOYMENT.id, TEST_DEPLOYMENT.name, TEST_DEPLOYMENT.type);
    expect(func !== null && func !== undefined, true);

    const inputNames = func.inputs.map((p) => p.name);
    expect(inputNames.length, EXPECTED_INPUT_NAMES.length);
    for (const expected of EXPECTED_INPUT_NAMES)
      expect(inputNames.includes(expected), true);

    const outputNames = func.outputs.map((p) => p.name);
    expect(outputNames.length, EXPECTED_OUTPUT_NAMES.length);
    for (const expected of EXPECTED_OUTPUT_NAMES)
      expect(outputNames.includes(expected), true);
  });

  test('returns workflow image URL for the test pipeline', async () => {
    const deployments: KnimeDeployment[] = await funcs.knimeListDeployments('rest');
    const dep = deployments.find((d) => d.id === TEST_DEPLOYMENT.id);
    expect(dep !== undefined, true);
    expect(dep!.workflowId !== undefined && dep!.workflowId !== null, true);

    const url: string | null = await funcs.knimeGetWorkflowImageUrl(dep!.workflowId!);
    expect(url !== null && url !== undefined, true);
    expect(typeof url, 'string');
    expect(url!.includes('connectors/proxy'), true);

    // The image must actually load via <img> (fetch() fails for the KNIME SVG endpoint
    // per the JSDoc on getWorkflowImageUrl — that's why production uses an <img> element).
    await new Promise<void>((resolve, reject) => {
      const img = document.createElement('img');
      const timeoutId = setTimeout(() => reject(new Error('Image load timed out')), 15000);
      img.onload = () => { clearTimeout(timeoutId); resolve(); };
      img.onerror = () => { clearTimeout(timeoutId); reject(new Error('Image failed to load')); };
      img.src = url!;
    });
  }, {timeout: 30000});

  test('preview renders for the function', async () => {
    const func = await funcs.knimeGetOrRegisterFunc(
      TEST_DEPLOYMENT.id, TEST_DEPLOYMENT.name, TEST_DEPLOYMENT.type);
    const handler = DG.ObjectHandler.forEntity(func);
    expect(handler !== null && handler !== undefined, true);
    const preview = await handler!.renderPreview(func);
    expect(preview !== null && preview !== undefined, true);
  });

  test('runs pipeline and echoes inputs back', async () => {
    await funcs.knimeGetOrRegisterFunc(
      TEST_DEPLOYMENT.id, TEST_DEPLOYMENT.name, TEST_DEPLOYMENT.type);

    const csvText = await _package.files.readAsText('test_pipeline.csv');
    const inputDf = DG.DataFrame.fromCsv(csvText);
    const fileName = '111111.txt';
    const fileBytes = new TextEncoder().encode('test content');
    const fileInfo = DG.FileInfo.fromBytes(fileName, fileBytes);

    const result = await grok.functions.call('KnimeLink:Knime_Test_inputs', {
      input_file: fileInfo,
      json_input_my_int_value: 5,
      json_input_my_double_value: 5.5,
      json_input_my_string_value: 'my test string 111',
      table_input: inputDf,
      variable_input_race: 'Caucasian',
      variable_input_weight: 70,
      arbitary_variable_input: '{}',
    });
    expect(result !== null && result !== undefined, true);

    // table_output: same dataframe back, matching row count and key columns
    const outputDf = result.table_output as DG.DataFrame;
    expect(outputDf instanceof DG.DataFrame, true);
    expect(outputDf.rowCount, inputDf.rowCount);
    for (const col of ['USUBJID', 'AGE', 'SEX', 'RACE', 'HEIGHT', 'WEIGHT'])
      expect(outputDf.columns.contains(col), true);
    expect(outputDf.col('USUBJID')!.get(0), inputDf.col('USUBJID')!.get(0));
    expect(outputDf.col('AGE')!.get(0), inputDf.col('AGE')!.get(0));

    // json_output: the JSON-input values echoed back (format may be JSON or printable)
    const jsonOut = typeof result.json_output === 'string'
      ? result.json_output
      : JSON.stringify(result.json_output);
    expect(jsonOut.includes('my test string 111'), true);
    expect(jsonOut.includes('5.5'), true);
    expect(jsonOut.includes('5'), true);

    // file_output: the input file name (KNIME returns filename string per workflow design)
    const fileOut = String(result.file_output);
    expect(fileOut.length > 0, true);
    expect(fileOut.includes('111111') || fileOut.includes('input-file') || fileOut.includes('.txt'), true);
  }, {timeout: 60000});
}, {timeout: 60000});
