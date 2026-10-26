/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/matrix-plot/matrix-plot-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.matrix-plot, core#321]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/file-formats.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, boundTable, hasArea, hasNoArea, noErrors, readingAsRemembered, readingIs, readingReads, rememberReading, reportsNoError} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Matrix plot regressions of closed issues", () => {
  const session = feature(test, "features/viewers/matrix-plot/matrix-plot-regressions.feature", import.meta.url);
  test("The Table property re-tiles the matrix to a two-column table without slider dots, and back (core#321)", {tag: ["@viewers", "@realizes:viewers.matrix-plot", "@realizes:core#321"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens a table \"two-numbers\" with:", () => openTableOf(page, "two-numbers", [["a","b"],["1","10"],["2","30"],["3","20"],["4","40"]]), [["a","b"],["1","10"],["2","30"],["3","20"],["4","40"]]);
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user adds a matrix plot viewer", () => addViewer(page, "matrix plot"));
    await session.step(20, "Then matrix plot viewer should be bound to table \"demog-1000\"", () => boundTable(page, el("matrix plot viewer"), "demog-1000"));
    await session.step(21, "And the \"cells drawn\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells drawn", el("matrix plot viewer"), 16));
    await session.step(22, "And the \"x scroll shown\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "x scroll shown", el("matrix plot viewer"), "true"));
    await session.step(23, "And matrix plot viewer should have a \"x scroll slider\" area", () => hasArea(page, el("matrix plot viewer"), "x scroll slider"));
    await session.step(24, "When user remembers the \"x columns\" reading of matrix plot viewer", () => rememberReading(page, "x columns", el("matrix plot viewer")));
    await session.step(25, "And user clicks on settings icon of matrix plot viewer", () => clickOn(page, el("settings icon of matrix plot viewer")));
    await session.step(26, "And user selects \"two-numbers\" in \"Table\" property in context panel", () => selectIn(page, "two-numbers", el("\"Table\" property in context panel")));
    await session.step(27, "Then matrix plot viewer should be bound to table \"two-numbers\"", () => boundTable(page, el("matrix plot viewer"), "two-numbers"));
    await session.step(28, "And the \"x columns\" reading of matrix plot viewer should be 2", () => readingIs(page, "x columns", el("matrix plot viewer"), 2));
    await session.step(29, "And the \"cells drawn\" reading of matrix plot viewer should be 4", () => readingIs(page, "cells drawn", el("matrix plot viewer"), 4));
    await session.step(30, "And the \"blank cells\" reading of matrix plot viewer should be 0", () => readingIs(page, "blank cells", el("matrix plot viewer"), 0));
    await session.step(31, "And the \"x scroll shown\" reading of matrix plot viewer should be \"false\"", () => readingReads(page, "x scroll shown", el("matrix plot viewer"), "false"));
    await session.step(32, "And the \"y scroll shown\" reading of matrix plot viewer should be \"false\"", () => readingReads(page, "y scroll shown", el("matrix plot viewer"), "false"));
    await session.step(33, "And matrix plot viewer should not have a \"x scroll slider\" area", () => hasNoArea(page, el("matrix plot viewer"), "x scroll slider"));
    await session.step(34, "And matrix plot viewer should not have a \"y scroll slider\" area", () => hasNoArea(page, el("matrix plot viewer"), "y scroll slider"));
    await session.step(35, "And matrix plot viewer should report no error", () => reportsNoError(page, el("matrix plot viewer")));
    await session.step(36, "When user selects \"demog-1000\" in \"Table\" property in context panel", () => selectIn(page, "demog-1000", el("\"Table\" property in context panel")));
    await session.step(37, "Then matrix plot viewer should be bound to table \"demog-1000\"", () => boundTable(page, el("matrix plot viewer"), "demog-1000"));
    await session.step(38, "And the \"x columns\" reading of matrix plot viewer should be as remembered", () => readingAsRemembered(page, "x columns", el("matrix plot viewer")));
    await session.step(39, "And the \"cells drawn\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells drawn", el("matrix plot viewer"), 16));
    await session.step(40, "And matrix plot viewer should have a \"x scroll slider\" area", () => hasArea(page, el("matrix plot viewer"), "x scroll slider"));
    await session.step(41, "And no errors should have been logged", () => noErrors(page));
  });
});
