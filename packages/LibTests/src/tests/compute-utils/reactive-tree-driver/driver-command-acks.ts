import {category, test} from '@datagrok-libraries/test/src/test';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {Driver} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/Driver';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {InitPipeline} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/view/ViewCommunication';
import {FuncCallNode} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';

const config: PipelineConfiguration = {
  id: 'root',
  type: 'static',
  steps: [
    {id: 'step1', nqName: 'LibTests:TestAdd2'},
    {id: 'step2', nqName: 'LibTests:TestMul2'},
  ],
};

async function initMsg(): Promise<InitPipeline> {
  return {event: 'initPipeline', provider: '', config: await getProcessedConfig(config)};
}

category('ComputeUtils: Driver command acks', async () => {
  test('Failed command resolves null and logs the error', async () => {
    const driver = new Driver(true);
    const res = await driver.sendCommand({event: 'runStep', uuid: 'missing'});
    expectDeepEqual(res, null, {prefix: 'Resolved value'});
    expectDeepEqual(driver.logger.errors.length, 1, {prefix: 'Error count'});
    expectDeepEqual(driver.logger.errors[0].context, 'command:runStep', {prefix: 'Error context'});
    driver.close();
  });

  test('Queue survives a failed command', async () => {
    const driver = new Driver(true);
    await driver.sendCommand({event: 'runStep', uuid: 'missing'});
    const state = await driver.sendCommand(await initMsg());
    expectDeepEqual(state instanceof StateTree, true, {prefix: 'Init after failure'});
    expectDeepEqual(driver.currentState$.value != null, true, {prefix: 'State set'});
    driver.close();
  });

  test('Successful command resolves the last emitted value', async () => {
    const driver = new Driver(true);
    const state = await driver.sendCommand(await initMsg());
    expectDeepEqual(state instanceof StateTree, true, {prefix: 'Resolved value'});
    driver.close();
  });

  test('Acks follow the command order', async () => {
    const driver = new Driver(true);
    const acks: number[] = [];
    driver.commandAcks$.subscribe((ack) => acks.push(ack.cid));
    const p1 = driver.sendCommand({event: 'runStep', uuid: 'missing1'});
    const p2 = driver.sendCommand({event: 'runStep', uuid: 'missing2'});
    await Promise.all([p1, p2]);
    expectDeepEqual(acks, [1, 2], {prefix: 'Ack order'});
    driver.close();
  });

  test('Locked driver resolves null without executing', async () => {
    const driver = new Driver(true);
    const msg = await initMsg();
    driver.globalROLocked$.next(true);
    const res = await driver.sendCommand(msg);
    expectDeepEqual(res, null, {prefix: 'Resolved value'});
    expectDeepEqual(driver.currentState$.value == null, true, {prefix: 'No state set'});
    expectDeepEqual(driver.logger.errors.length, 0, {prefix: 'No errors'});
    driver.close();
  });

  test('Save on a mock tree acks the error and resolves null', async () => {
    const driver = new Driver(true);
    await driver.sendCommand(await initMsg());
    const res = await driver.sendCommand({event: 'savePipeline'});
    expectDeepEqual(res, null, {prefix: 'Resolved value'});
    expectDeepEqual(driver.logger.errors.length >= 1, true, {prefix: 'Error logged'});
    driver.close();
  });

  test('Run step command resolves after the run completes', async () => {
    const driver = new Driver(true);
    const tree = await driver.sendCommand(await initMsg()) as StateTree;
    const step = tree.nodeTree.getNode([{idx: 0}]).getItem() as FuncCallNode;
    step.getStateStore().setState('a', 1);
    step.getStateStore().setState('b', 2);
    const stepUuid = (driver.currentState$.value as any).steps[0].uuid;
    await driver.sendCommand({event: 'runStep', uuid: stepUuid, mockResults: {res: 5}});
    const firstError = driver.logger.errors[0];
    expectDeepEqual(firstError ? `${firstError.context}: ${firstError.message}` : 'none', 'none', {prefix: 'No errors'});
    driver.close();
  });
});
