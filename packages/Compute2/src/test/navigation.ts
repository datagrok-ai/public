import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {
  PipelineState,
  PipelineInstanceRuntimeData,
  PipelineStateStatic,
  PipelineStateDynamic,
  StepFunCallState,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {findNextStep, findPrevStep, findNextSubStep} from '../utils';

function mockFuncCall(uuid: string, opts?: {isReadonly?: boolean}): StepFunCallState {
  return {
    type: 'funccall',
    uuid,
    configId: uuid,
    isReadonly: opts?.isReadonly ?? false,
  };
}

const defaultPipelineRuntimeData: PipelineInstanceRuntimeData = {
  actions: undefined,
  approversGroup: undefined,
  disableHistory: false,
  customExports: undefined,
  forceNavigate: false,
};

function mockStaticPipeline(
  uuid: string,
  steps: PipelineState[],
  opts?: {isReadonly?: boolean; isActionStep?: boolean; forceNavigate?: boolean},
): PipelineStateStatic<StepFunCallState, PipelineInstanceRuntimeData> {
  return {
    type: 'static',
    uuid,
    configId: uuid,
    friendlyName: undefined,
    version: undefined,
    nqName: undefined,
    isReadonly: opts?.isReadonly ?? false,
    steps,
    isActionStep: opts?.isActionStep,
    ...defaultPipelineRuntimeData,
    forceNavigate: opts?.forceNavigate ?? false,
  };
}

function mockDynamicPipeline(
  uuid: string,
  steps: PipelineState[],
  opts?: {isReadonly?: boolean; forceNavigate?: boolean; type?: 'dynamic' | 'parallel' | 'sequential'},
): PipelineStateDynamic<StepFunCallState, PipelineInstanceRuntimeData> {
  return {
    type: opts?.type ?? 'dynamic',
    uuid,
    configId: uuid,
    friendlyName: undefined,
    version: undefined,
    nqName: undefined,
    isReadonly: opts?.isReadonly ?? false,
    steps,
    stepTypes: [],
    ...defaultPipelineRuntimeData,
    forceNavigate: opts?.forceNavigate ?? false,
  };
}

function collectSequence(startUuid: string, tree: PipelineState, direction: 'forward' | 'backward'): string[] {
  const fn = direction === 'forward' ? findNextStep : findPrevStep;
  const result: string[] = [];
  let uuid = startUuid;
  for (let i = 0; i < 20; i++) {
    const next = fn(uuid, tree);
    if (!next) break;
    result.push(next.state.uuid);
    uuid = next.state.uuid;
  }
  return result;
}

function expectSequencesMatch(forward: string[], backward: string[]) {
  expect(JSON.stringify(forward), JSON.stringify(backward));
}

// ============================================================
// Default: pipelines with children are SKIPPED
// ============================================================
//
//   Root (static, has children -> skipped)
//     |- Step1 (funccall)
//     |- PipelineA (static, has children -> skipped)
//     |   |- Step2 (funccall)
//     |   +- Step3 (funccall)
//     |- ActionStep (static, no children -> navigable)
//     +- Step4 (funccall)
//
// Navigable sequence: step1, step2, step3, action1, step4

function buildDefaultTree(): PipelineState {
  return mockStaticPipeline('root', [
    mockFuncCall('step1'),
    mockStaticPipeline('pipeA', [
      mockFuncCall('step2'),
      mockFuncCall('step3'),
    ]),
    mockStaticPipeline('action1', [], {isActionStep: true}),
    mockFuncCall('step4'),
  ]);
}

category('Navigation: default (pipelines skipped)', () => {
  test('forward: full sequence skips pipelines with children', async () => {
    const tree = buildDefaultTree();
    expectSequencesMatch(
      collectSequence('step1', tree, 'forward'),
      ['step2', 'step3', 'action1', 'step4'],
    );
  });

  test('forward: action step is navigable', async () => {
    const tree = buildDefaultTree();
    expect(findNextStep('step3', tree)?.state.uuid, 'action1');
  });

  test('forward: last step returns undefined', async () => {
    const tree = buildDefaultTree();
    expect(findNextStep('step4', tree) == null, true);
  });

  test('forward: from skipped pipeline finds first navigable child', async () => {
    const tree = buildDefaultTree();
    expect(findNextStep('pipeA', tree)?.state.uuid, 'step2');
  });

  test('forward: from skipped root finds first navigable child', async () => {
    const tree = buildDefaultTree();
    expect(findNextStep('root', tree)?.state.uuid, 'step1');
  });

  test('backward: full sequence skips pipelines with children', async () => {
    const tree = buildDefaultTree();
    expectSequencesMatch(
      collectSequence('step4', tree, 'backward'),
      ['action1', 'step3', 'step2', 'step1'],
    );
  });

  test('backward: from first navigable step returns undefined', async () => {
    const tree = buildDefaultTree();
    expect(findPrevStep('step1', tree) == null, true);
  });

  test('backward: from skipped pipeline finds previous navigable step', async () => {
    const tree = buildDefaultTree();
    expect(findPrevStep('pipeA', tree)?.state.uuid, 'step1');
  });

  test('backward: from root returns undefined', async () => {
    const tree = buildDefaultTree();
    expect(findPrevStep('root', tree) == null, true);
  });

  test('backward is exact reverse of forward', async () => {
    const tree = buildDefaultTree();
    const forward = ['step1', ...collectSequence('step1', tree, 'forward')];
    const backward = collectSequence('step4', tree, 'backward');
    backward.reverse();
    backward.push('step4');
    expectSequencesMatch(forward, backward);
  });

  test('readonly steps are skipped in forward', async () => {
    const tree = mockStaticPipeline('root', [
      mockFuncCall('step1'),
      mockFuncCall('step2', {isReadonly: true}),
      mockFuncCall('step3'),
    ]);
    expect(findNextStep('step1', tree)?.state.uuid, 'step3');
  });

  test('readonly steps are skipped in backward', async () => {
    const tree = mockStaticPipeline('root', [
      mockFuncCall('step1'),
      mockFuncCall('step2', {isReadonly: true}),
      mockFuncCall('step3'),
    ]);
    expect(findPrevStep('step3', tree)?.state.uuid, 'step1');
  });
});

// ============================================================
// forceNavigate: pipelines become non-skippable
// ============================================================
//
//   Root (static, forceNavigate -> navigable)
//     |- Step1 (funccall)
//     |- PipelineA (static, forceNavigate -> navigable)
//     |   |- Step2 (funccall)
//     |   +- Step3 (funccall)
//     |- ActionStep (static, no children -> always navigable)
//     +- Step4 (funccall)
//
// Navigable sequence: root, step1, pipeA, step2, step3, action1, step4

function buildForceNavTree(): PipelineState {
  return mockStaticPipeline('root', [
    mockFuncCall('step1'),
    mockStaticPipeline('pipeA', [
      mockFuncCall('step2'),
      mockFuncCall('step3'),
    ], {forceNavigate: true}),
    mockStaticPipeline('action1', [], {isActionStep: true}),
    mockFuncCall('step4'),
  ], {forceNavigate: true});
}

category('Navigation: forceNavigate (pipelines included)', () => {
  test('forward: full sequence includes pipelines with forceNavigate', async () => {
    const tree = buildForceNavTree();
    expectSequencesMatch(
      collectSequence('root', tree, 'forward'),
      ['step1', 'pipeA', 'step2', 'step3', 'action1', 'step4'],
    );
  });

  test('forward: step before pipeline goes to pipeline', async () => {
    const tree = buildForceNavTree();
    expect(findNextStep('step1', tree)?.state.uuid, 'pipeA');
  });

  test('forward: pipeline to first child', async () => {
    const tree = buildForceNavTree();
    expect(findNextStep('pipeA', tree)?.state.uuid, 'step2');
  });

  test('forward: root to first child', async () => {
    const tree = buildForceNavTree();
    expect(findNextStep('root', tree)?.state.uuid, 'step1');
  });

  test('backward: full sequence includes pipelines with forceNavigate', async () => {
    const tree = buildForceNavTree();
    expectSequencesMatch(
      collectSequence('step4', tree, 'backward'),
      ['action1', 'step3', 'step2', 'pipeA', 'step1', 'root'],
    );
  });

  test('backward: first child to parent pipeline', async () => {
    const tree = buildForceNavTree();
    expect(findPrevStep('step2', tree)?.state.uuid, 'pipeA');
  });

  test('backward: pipeline to previous sibling', async () => {
    const tree = buildForceNavTree();
    expect(findPrevStep('pipeA', tree)?.state.uuid, 'step1');
  });

  test('backward: first step to root', async () => {
    const tree = buildForceNavTree();
    expect(findPrevStep('step1', tree)?.state.uuid, 'root');
  });

  test('backward is exact reverse of forward', async () => {
    const tree = buildForceNavTree();
    const forward = ['root', ...collectSequence('root', tree, 'forward')];
    const backward = collectSequence('step4', tree, 'backward');
    backward.reverse();
    backward.push('step4');
    expectSequencesMatch(forward, backward);
  });
});

// ============================================================
// Mixed: only some pipelines have forceNavigate
// ============================================================
//
//   Root (static, no flag -> skipped)
//     |- Step1 (funccall)
//     |- PipelineA (static, forceNavigate -> navigable)
//     |   |- Step2 (funccall)
//     |   +- Step3 (funccall)
//     |- PipelineB (static, no flag -> skipped)
//     |   +- Step4 (funccall)
//     +- Step5 (funccall)
//
// Navigable sequence: step1, pipeA, step2, step3, step4, step5

function buildMixedTree(): PipelineState {
  return mockStaticPipeline('root', [
    mockFuncCall('step1'),
    mockStaticPipeline('pipeA', [
      mockFuncCall('step2'),
      mockFuncCall('step3'),
    ], {forceNavigate: true}),
    mockStaticPipeline('pipeB', [
      mockFuncCall('step4'),
    ]),
    mockFuncCall('step5'),
  ]);
}

category('Navigation: mixed forceNavigate', () => {
  test('forward: includes only pipelines with forceNavigate', async () => {
    const tree = buildMixedTree();
    expectSequencesMatch(
      collectSequence('step1', tree, 'forward'),
      ['pipeA', 'step2', 'step3', 'step4', 'step5'],
    );
  });

  test('backward: includes only pipelines with forceNavigate', async () => {
    const tree = buildMixedTree();
    expectSequencesMatch(
      collectSequence('step5', tree, 'backward'),
      ['step4', 'step3', 'step2', 'pipeA', 'step1'],
    );
  });

  test('backward: from skipped pipelineB finds previous navigable', async () => {
    const tree = buildMixedTree();
    expect(findPrevStep('pipeB', tree)?.state.uuid, 'step3');
  });

  test('forward: from skipped pipelineB finds first navigable child', async () => {
    const tree = buildMixedTree();
    expect(findNextStep('pipeB', tree)?.state.uuid, 'step4');
  });

  test('backward is exact reverse of forward', async () => {
    const tree = buildMixedTree();
    const forward = ['step1', ...collectSequence('step1', tree, 'forward')];
    const backward = collectSequence('step5', tree, 'backward');
    backward.reverse();
    backward.push('step5');
    expectSequencesMatch(forward, backward);
  });
});

// ============================================================
// Empty pipelines and action steps: always navigable
// ============================================================

category('Navigation: empty pipelines and action steps', () => {
  test('action step is always navigable without forceNavigate', async () => {
    const tree = mockStaticPipeline('root', [
      mockFuncCall('step1'),
      mockStaticPipeline('action1', [], {isActionStep: true}),
      mockFuncCall('step2'),
    ]);
    expectSequencesMatch(
      collectSequence('step1', tree, 'forward'),
      ['action1', 'step2'],
    );
  });

  test('empty pipeline (no steps) is always navigable without forceNavigate', async () => {
    const tree = mockStaticPipeline('root', [
      mockFuncCall('step1'),
      mockStaticPipeline('empty', []),
      mockFuncCall('step2'),
    ]);
    expectSequencesMatch(
      collectSequence('step1', tree, 'forward'),
      ['empty', 'step2'],
    );
  });

  test('empty dynamic pipeline is navigable', async () => {
    const tree = mockStaticPipeline('root', [
      mockFuncCall('step1'),
      mockDynamicPipeline('dynEmpty', []),
      mockFuncCall('step2'),
    ]);
    expectSequencesMatch(
      collectSequence('step1', tree, 'forward'),
      ['dynEmpty', 'step2'],
    );
  });

  test('backward through action step', async () => {
    const tree = mockStaticPipeline('root', [
      mockFuncCall('step1'),
      mockStaticPipeline('action1', [], {isActionStep: true}),
      mockFuncCall('step2'),
    ]);
    expectSequencesMatch(
      collectSequence('step2', tree, 'backward'),
      ['action1', 'step1'],
    );
  });
});

// ============================================================
// Deep nesting
// ============================================================

category('Navigation: deep nesting', () => {
  test('forward: deeply nested with forceNavigate', async () => {
    const tree = mockStaticPipeline('root', [
      mockStaticPipeline('pipeA', [
        mockFuncCall('step1'),
        mockStaticPipeline('pipeB', [
          mockFuncCall('step2'),
          mockFuncCall('step3'),
        ], {forceNavigate: true}),
      ], {forceNavigate: true}),
      mockFuncCall('step4'),
    ]);
    expectSequencesMatch(
      collectSequence('pipeA', tree, 'forward'),
      ['step1', 'pipeB', 'step2', 'step3', 'step4'],
    );
  });

  test('backward: deeply nested with forceNavigate', async () => {
    const tree = mockStaticPipeline('root', [
      mockStaticPipeline('pipeA', [
        mockFuncCall('step1'),
        mockStaticPipeline('pipeB', [
          mockFuncCall('step2'),
          mockFuncCall('step3'),
        ], {forceNavigate: true}),
      ], {forceNavigate: true}),
      mockFuncCall('step4'),
    ]);
    expectSequencesMatch(
      collectSequence('step4', tree, 'backward'),
      ['step3', 'step2', 'pipeB', 'step1', 'pipeA'],
    );
  });

  test('forward: deeply nested without forceNavigate skips inner pipelines', async () => {
    const tree = mockStaticPipeline('root', [
      mockStaticPipeline('pipeA', [
        mockFuncCall('step1'),
        mockStaticPipeline('pipeB', [
          mockFuncCall('step2'),
        ]),
      ]),
      mockFuncCall('step3'),
    ]);
    expectSequencesMatch(
      collectSequence('step1', tree, 'forward'),
      ['step2', 'step3'],
    );
  });
});

// ============================================================
// findNextSubStep (entering a pipeline)
// ============================================================

category('Navigation: findNextSubStep', () => {
  test('finds first funccall child', async () => {
    const tree = buildDefaultTree();
    expect(findNextSubStep(tree)?.state.uuid, 'step1');
  });

  test('finds first child in nested pipeline', async () => {
    const pipeA = mockStaticPipeline('pipeA', [
      mockFuncCall('step2'),
      mockFuncCall('step3'),
    ]);
    expect(findNextSubStep(pipeA)?.state.uuid, 'step2');
  });

  test('action step returns itself (no children)', async () => {
    const action = mockStaticPipeline('action1', [], {isActionStep: true});
    expect(findNextSubStep(action)?.state.uuid, 'action1');
  });

  test('empty pipeline returns itself', async () => {
    const empty = mockStaticPipeline('empty', []);
    expect(findNextSubStep(empty)?.state.uuid, 'empty');
  });

  test('pipeline with forceNavigate returns itself', async () => {
    const pipe = mockStaticPipeline('pipe', [
      mockFuncCall('step1'),
    ], {forceNavigate: true});
    expect(findNextSubStep(pipe)?.state.uuid, 'pipe');
  });
});

// ============================================================
// Dynamic pipeline variants
// ============================================================

category('Navigation: dynamic pipelines', () => {
  test('dynamic pipeline with children is skipped by default', async () => {
    const tree = mockStaticPipeline('root', [
      mockFuncCall('step1'),
      mockDynamicPipeline('dynPipe', [
        mockFuncCall('step2'),
      ]),
      mockFuncCall('step3'),
    ]);
    expectSequencesMatch(
      collectSequence('step1', tree, 'forward'),
      ['step2', 'step3'],
    );
  });

  test('dynamic pipeline with forceNavigate is included', async () => {
    const tree = mockStaticPipeline('root', [
      mockFuncCall('step1'),
      mockDynamicPipeline('dynPipe', [
        mockFuncCall('step2'),
      ], {forceNavigate: true}),
      mockFuncCall('step3'),
    ]);
    expectSequencesMatch(
      collectSequence('step1', tree, 'forward'),
      ['dynPipe', 'step2', 'step3'],
    );
  });

  test('backward through dynamic pipeline with forceNavigate', async () => {
    const tree = mockStaticPipeline('root', [
      mockFuncCall('step1'),
      mockDynamicPipeline('dynPipe', [
        mockFuncCall('step2'),
      ], {forceNavigate: true}),
      mockFuncCall('step3'),
    ]);
    expect(findPrevStep('step2', tree)?.state.uuid, 'dynPipe');
  });
});
