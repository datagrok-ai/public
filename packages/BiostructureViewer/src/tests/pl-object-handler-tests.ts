import * as DG from 'datagrok-api/dg';

import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {PlDiagramObjectHandler, PROLIF_SOURCE_TAG} from '../utils/pl-object-handler';


// `PlDiagramObjectHandler` is the entry point for the Protein-Ligand
// Interactions context panel after the move off `@grok.decorators.panel`.
// `isApplicable` is the whole contract — it decides whether the handler
// claims a clicked cell. If it ever returns true for non-PL cells, the
// panel fires everywhere; if it returns false for our own cells, the
// panel never fires. The three tests below lock in: (1) the positive case
// keyed off the `.%prolif-source` tag set by `runPlBatch`; (2) the negative
// case for cells without the tag — critical because we share the
// `rawPng` semType with PowerGrid's generic renderer; (3) safe rejection
// of any non-SemanticValue input the platform may pass.
category('PLObjectHandler', () => {
  test('isApplicable: matches cells whose column has .%prolif-source tag', async () => {
    const handler = new PlDiagramObjectHandler();
    const df = DG.DataFrame.fromColumns([DG.Column.string('PL Diagram', 1)]);
    df.col('PL Diagram')!.tags[PROLIF_SOURCE_TAG] = 'protein';
    const sv = DG.SemanticValue.fromTableCell(df.cell(0, 'PL Diagram'));
    expect(handler.isApplicable(sv), true);
  });

  test('isApplicable: rejects cells without the tag', async () => {
    // Cell is a real Datagrok cell with a non-null dart and column — the
    // SemType could even be `rawPng` (shared with PowerGrid) — but the
    // absence of `.%prolif-source` is what keeps the handler from claiming
    // unrelated rawPng cells from other features.
    const handler = new PlDiagramObjectHandler();
    const df = DG.DataFrame.fromColumns([DG.Column.string('something else', 1)]);
    const sv = DG.SemanticValue.fromTableCell(df.cell(0, 'something else'));
    expect(handler.isApplicable(sv), false);
  });

  test('isApplicable: rejects non-SemanticValue inputs', async () => {
    // Datagrok calls isApplicable with many object shapes — strings, plain
    // objects, null, numbers — when resolving the right handler for an
    // arbitrary `grok.shell.o`. We must not throw and must return false
    // for anything that isn't a SemanticValue.
    const handler = new PlDiagramObjectHandler();
    expect(handler.isApplicable(null), false);
    expect(handler.isApplicable(undefined), false);
    expect(handler.isApplicable('a string'), false);
    expect(handler.isApplicable({some: 'plain object'}), false);
    expect(handler.isApplicable(42), false);
  });
});
