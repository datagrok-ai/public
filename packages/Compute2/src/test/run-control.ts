import {category, test} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {pressControl, runControl, RunState, READONLY_REASON} from '../utils';

const state = (patch: Partial<RunState> = {}): RunState => ({
  isOutputOutdated: true, allowRerun: false, autorun: false,
  isRunnable: true, isRunning: false, isReadonly: false, isBlocked: false, ...patch,
});

category('RFV: run control', () => {
  test('Run shows for outdated outputs or allowed reruns, and never on autorun', async () => {
    expectDeepEqual(runControl(state()), {shown: true, disabledReason: undefined}, {prefix: 'Outdated'});
    expectDeepEqual(runControl(state({isOutputOutdated: false})).shown, false, {prefix: 'Up to date'});
    expectDeepEqual(runControl(state({isOutputOutdated: false, allowRerun: true})),
      {shown: true, disabledReason: undefined}, {prefix: 'Rerun'});
    expectDeepEqual(runControl(state({autorun: true, allowRerun: true})).shown, false, {prefix: 'Autorun'});
  });

  test('A disabled Run gives the first reason: read-only, running, busy, invalid', async () => {
    const reason = (patch: Partial<RunState>) => runControl(state(patch)).disabledReason;
    const all = {isReadonly: true, isRunning: true, isBlocked: true, isRunnable: false};
    expectDeepEqual(reason(all), READONLY_REASON, {prefix: 'Read-only first'});
    expectDeepEqual(reason({...all, isReadonly: false}), 'It is already running.', {prefix: 'Running'});
    expectDeepEqual(reason({isBlocked: true, isRunnable: false})?.startsWith('The workflow is busy'), true,
      {prefix: 'Busy'});
    expectDeepEqual(reason({isRunnable: false}), 'Fix the invalid inputs first.', {prefix: 'Invalid'});
  });

  test('Pressing a control acts, warns instead when disabled, and reports a hidden one', async () => {
    let acted = 0;
    const act = () => acted++;
    expectDeepEqual(pressControl({shown: false}, act), false, {prefix: 'Hidden'});
    expectDeepEqual(pressControl({shown: true, disabledReason: 'No'}, act), true, {prefix: 'Disabled'});
    expectDeepEqual(pressControl({shown: true}, act), true, {prefix: 'Enabled'});
    expectDeepEqual(acted, 1, {prefix: 'Acted once'});
  });
});
