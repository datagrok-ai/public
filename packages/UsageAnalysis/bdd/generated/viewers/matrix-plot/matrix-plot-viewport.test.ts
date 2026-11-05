/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/matrix-plot/matrix-plot-viewport.feature
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
import {columnCount} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {addCalculated, removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, openTableOf, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, boundTable, hasArea, hasNoArea, noErrors, readingAsRemembered, readingBetween, readingHigher, readingIs, readingReads, readingsEqual, rememberReading, reportsNoError, setProperties} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {cellsWideTall, dragHandleToEnd} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Matrix plot — the viewport over a long column list, and the 250-cell cap", () => {
  const session = feature(test, "features/viewers/matrix-plot/matrix-plot-viewport.feature", import.meta.url);
  test("Matrix plot — the viewport over a long column list, and the 250-cell cap", {tag: ["@journey", "@viewers", "@realizes:viewers.matrix-plot", "@realizes:core#321"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(25, "And user adds a matrix plot viewer", () => addViewer(page, "matrix plot"));
    await session.step(26, "Then the \"cells\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells", el("matrix plot viewer"), 16));
    await session.step(27, "And the \"viewport limit\" reading of matrix plot viewer should be 250", () => readingIs(page, "viewport limit", el("matrix plot viewer"), 250));
    await session.step(28, "And the \"viewport rejected\" reading of matrix plot viewer should be \"false\"", () => readingReads(page, "viewport rejected", el("matrix plot viewer"), "false"));
    await session.step(29, "And matrix plot viewer should have an \"x scroll slider\" area", () => hasArea(page, el("matrix plot viewer"), "x scroll slider"));
    await session.step(30, "And matrix plot viewer should have a \"y scroll slider\" area", () => hasArea(page, el("matrix plot viewer"), "y scroll slider"));
    await run.scenario("With four columns the viewport already holds them all and a slider drag is inert", async () => {
      await session.step(39, "Then the cells of matrix plot viewer should be 4 wide and 4 tall", () => cellsWideTall(page, el("matrix plot viewer"), 4, 4));
      await session.step(40, "And matrix plot viewer should have an \"x label STARTED\" area", () => hasArea(page, el("matrix plot viewer"), "x label STARTED"));
      await session.step(41, "And matrix plot viewer should have a \"cell AGE x AGE\" area", () => hasArea(page, el("matrix plot viewer"), "cell AGE x AGE"));
      await session.step(42, "When user drags the max handle of the \"x\" scroll slider of matrix plot viewer to its end", () => dragHandleToEnd(page, "max", "x", el("matrix plot viewer"), "end"));
      await session.step(43, "Then the cells of matrix plot viewer should be 4 wide and 4 tall", () => cellsWideTall(page, el("matrix plot viewer"), 4, 4));
      await session.step(44, "And the \"cells\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells", el("matrix plot viewer"), 16));
      await session.step(45, "And the \"cells drawn\" and \"cells\" readings of matrix plot viewer should be the same", () => readingsEqual(page, "cells drawn", "cells", el("matrix plot viewer")));
      await session.step(46, "And the \"viewport rejected\" reading of matrix plot viewer should be \"false\"", () => readingReads(page, "viewport rejected", el("matrix plot viewer"), "false"));
      await session.step(47, "When user drags the max handle of the \"y\" scroll slider of matrix plot viewer to its end", () => dragHandleToEnd(page, "max", "y", el("matrix plot viewer"), "end"));
      await session.step(48, "Then the cells of matrix plot viewer should be 4 wide and 4 tall", () => cellsWideTall(page, el("matrix plot viewer"), 4, 4));
      await session.step(49, "And matrix plot viewer should have an \"x label STARTED\" area", () => hasArea(page, el("matrix plot viewer"), "x label STARTED"));
      await session.step(50, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Sixteen columns start at a five by five viewport and the sliders open it up", async () => {
      await session.step(53, "Given user adds a calculated column \"MP_FX_1\" with formula \"${AGE} + 1\"", () => addCalculated(page, "MP_FX_1", "${AGE} + 1"));
      await session.step(54, "And user adds a calculated column \"MP_FX_2\" with formula \"${AGE} + 2\"", () => addCalculated(page, "MP_FX_2", "${AGE} + 2"));
      await session.step(55, "And user adds a calculated column \"MP_FX_3\" with formula \"${AGE} + 3\"", () => addCalculated(page, "MP_FX_3", "${AGE} + 3"));
      await session.step(56, "And user adds a calculated column \"MP_FX_4\" with formula \"${AGE} + 4\"", () => addCalculated(page, "MP_FX_4", "${AGE} + 4"));
      await session.step(57, "And user adds a calculated column \"MP_FX_5\" with formula \"${AGE} + 5\"", () => addCalculated(page, "MP_FX_5", "${AGE} + 5"));
      await session.step(58, "And user adds a calculated column \"MP_FX_6\" with formula \"${AGE} + 6\"", () => addCalculated(page, "MP_FX_6", "${AGE} + 6"));
      await session.step(59, "And user adds a calculated column \"MP_FX_7\" with formula \"${AGE} + 7\"", () => addCalculated(page, "MP_FX_7", "${AGE} + 7"));
      await session.step(60, "And user adds a calculated column \"MP_FX_8\" with formula \"${AGE} + 8\"", () => addCalculated(page, "MP_FX_8", "${AGE} + 8"));
      await session.step(61, "And user adds a calculated column \"MP_FX_9\" with formula \"${AGE} + 9\"", () => addCalculated(page, "MP_FX_9", "${AGE} + 9"));
      await session.step(62, "And user adds a calculated column \"MP_FX_10\" with formula \"${AGE} + 10\"", () => addCalculated(page, "MP_FX_10", "${AGE} + 10"));
      await session.step(63, "And user adds a calculated column \"MP_FX_11\" with formula \"${AGE} + 11\"", () => addCalculated(page, "MP_FX_11", "${AGE} + 11"));
      await session.step(64, "And user adds a calculated column \"MP_FX_12\" with formula \"${AGE} + 12\"", () => addCalculated(page, "MP_FX_12", "${AGE} + 12"));
      await session.step(65, "When user sets properties of matrix plot viewer:", () => setProperties(page, el("matrix plot viewer"), [["xColumnNames","AGE, HEIGHT, WEIGHT, STARTED, MP_FX_1, MP_FX_2, MP_FX_3, MP_FX_4, MP_FX_5, MP_FX_6, MP_FX_7, MP_FX_8, MP_FX_9, MP_FX_10, MP_FX_11, MP_FX_12"],["yColumnNames","AGE, HEIGHT, WEIGHT, STARTED, MP_FX_1, MP_FX_2, MP_FX_3, MP_FX_4, MP_FX_5, MP_FX_6, MP_FX_7, MP_FX_8, MP_FX_9, MP_FX_10, MP_FX_11, MP_FX_12"]]), [["xColumnNames","AGE, HEIGHT, WEIGHT, STARTED, MP_FX_1, MP_FX_2, MP_FX_3, MP_FX_4, MP_FX_5, MP_FX_6, MP_FX_7, MP_FX_8, MP_FX_9, MP_FX_10, MP_FX_11, MP_FX_12"],["yColumnNames","AGE, HEIGHT, WEIGHT, STARTED, MP_FX_1, MP_FX_2, MP_FX_3, MP_FX_4, MP_FX_5, MP_FX_6, MP_FX_7, MP_FX_8, MP_FX_9, MP_FX_10, MP_FX_11, MP_FX_12"]]);
      await session.step(68, "Then the \"x columns\" reading of matrix plot viewer should be 16", () => readingIs(page, "x columns", el("matrix plot viewer"), 16));
      await session.step(69, "And the \"y columns\" reading of matrix plot viewer should be 16", () => readingIs(page, "y columns", el("matrix plot viewer"), 16));
      await session.step(70, "And the cells of matrix plot viewer should be 5 wide and 5 tall", () => cellsWideTall(page, el("matrix plot viewer"), 5, 5));
      await session.step(71, "And the \"cells\" reading of matrix plot viewer should be 25", () => readingIs(page, "cells", el("matrix plot viewer"), 25));
      await session.step(72, "And the \"viewport rejected\" reading of matrix plot viewer should be \"false\"", () => readingReads(page, "viewport rejected", el("matrix plot viewer"), "false"));
      await session.step(73, "When user drags the max handle of the \"x\" scroll slider of matrix plot viewer to its end", () => dragHandleToEnd(page, "max", "x", el("matrix plot viewer"), "end"));
      await session.step(74, "Then the \"columns\" reading of matrix plot viewer should be 16", () => readingIs(page, "columns", el("matrix plot viewer"), 16));
      await session.step(75, "And the \"cells\" reading of matrix plot viewer should be 80", () => readingIs(page, "cells", el("matrix plot viewer"), 80));
      await session.step(76, "And the \"cells drawn\" and \"cells\" readings of matrix plot viewer should be the same", () => readingsEqual(page, "cells drawn", "cells", el("matrix plot viewer")));
      await session.step(77, "And the \"viewport rejected\" reading of matrix plot viewer should be \"false\"", () => readingReads(page, "viewport rejected", el("matrix plot viewer"), "false"));
      await session.step(78, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The viewport stops short of the whole grid and records that it refused the rest", async () => {
      await session.step(86, "Then the \"cells\" reading of matrix plot viewer should be 80", () => readingIs(page, "cells", el("matrix plot viewer"), 80));
      await session.step(87, "When user drags the max handle of the \"y\" scroll slider of matrix plot viewer to its end", () => dragHandleToEnd(page, "max", "y", el("matrix plot viewer"), "end"));
      await session.step(88, "Then the \"cells\" reading of matrix plot viewer should be higher than before", () => readingHigher(page, "cells", el("matrix plot viewer")));
      await session.step(89, "And the \"cells\" reading of matrix plot viewer should be between 96 and 250", () => readingBetween(page, "cells", el("matrix plot viewer"), 96, 250));
      await session.step(90, "And the \"columns\" reading of matrix plot viewer should be 16", () => readingIs(page, "columns", el("matrix plot viewer"), 16));
      await session.step(91, "And the \"y columns\" reading of matrix plot viewer should be 16", () => readingIs(page, "y columns", el("matrix plot viewer"), 16));
      await session.step(92, "And the \"rows\" reading of matrix plot viewer should be between 6 and 15", () => readingBetween(page, "rows", el("matrix plot viewer"), 6, 15));
      await session.step(93, "And the \"viewport rejected\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "viewport rejected", el("matrix plot viewer"), "true"));
      await session.step(94, "And the \"cells drawn\" and \"cells\" readings of matrix plot viewer should be the same", () => readingsEqual(page, "cells drawn", "cells", el("matrix plot viewer")));
      await session.step(95, "And the \"error\" reading of matrix plot viewer should be \"\"", () => readingReads(page, "error", el("matrix plot viewer"), ""));
      await session.step(96, "When user sets properties of matrix plot viewer:", () => setProperties(page, el("matrix plot viewer"), [["xColumnNames","AGE, HEIGHT, WEIGHT, STARTED"],["yColumnNames","AGE, HEIGHT, WEIGHT, STARTED"]]), [["xColumnNames","AGE, HEIGHT, WEIGHT, STARTED"],["yColumnNames","AGE, HEIGHT, WEIGHT, STARTED"]]);
      await session.step(99, "Then the \"cells\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells", el("matrix plot viewer"), 16));
      await session.step(100, "When user removes \"MP_FX_1\" column", () => removeColumn(page, "MP_FX_1"));
      await session.step(101, "And user removes \"MP_FX_2\" column", () => removeColumn(page, "MP_FX_2"));
      await session.step(102, "And user removes \"MP_FX_3\" column", () => removeColumn(page, "MP_FX_3"));
      await session.step(103, "And user removes \"MP_FX_4\" column", () => removeColumn(page, "MP_FX_4"));
      await session.step(104, "And user removes \"MP_FX_5\" column", () => removeColumn(page, "MP_FX_5"));
      await session.step(105, "And user removes \"MP_FX_6\" column", () => removeColumn(page, "MP_FX_6"));
      await session.step(106, "And user removes \"MP_FX_7\" column", () => removeColumn(page, "MP_FX_7"));
      await session.step(107, "And user removes \"MP_FX_8\" column", () => removeColumn(page, "MP_FX_8"));
      await session.step(108, "And user removes \"MP_FX_9\" column", () => removeColumn(page, "MP_FX_9"));
      await session.step(109, "And user removes \"MP_FX_10\" column", () => removeColumn(page, "MP_FX_10"));
      await session.step(110, "And user removes \"MP_FX_11\" column", () => removeColumn(page, "MP_FX_11"));
      await session.step(111, "And user removes \"MP_FX_12\" column", () => removeColumn(page, "MP_FX_12"));
      await session.step(112, "Then the table should have 11 columns", () => columnCount(page, 11));
      await session.step(113, "And the \"cells\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells", el("matrix plot viewer"), 16));
      await session.step(114, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Table property re-tiles the matrix to a two-column table without slider dots, and back (core#321)", async () => {
      await session.step(118, "Given user opens a table \"two-numbers\" with:", () => openTableOf(page, "two-numbers", [["a","b"],["1","10"],["2","30"],["3","20"],["4","40"]]), [["a","b"],["1","10"],["2","30"],["3","20"],["4","40"]]);
      await session.step(124, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
      await session.step(125, "Then the \"cells drawn\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells drawn", el("matrix plot viewer"), 16));
      await session.step(126, "And the \"x scroll shown\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "x scroll shown", el("matrix plot viewer"), "true"));
      await session.step(127, "And matrix plot viewer should have a \"x scroll slider\" area", () => hasArea(page, el("matrix plot viewer"), "x scroll slider"));
      await session.step(128, "When user remembers the \"x columns\" reading of matrix plot viewer", () => rememberReading(page, "x columns", el("matrix plot viewer")));
      await session.step(129, "And user clicks on settings icon of matrix plot viewer", () => clickOn(page, el("settings icon of matrix plot viewer")));
      await session.step(130, "And user selects \"two-numbers\" in \"Table\" property in context panel", () => selectIn(page, "two-numbers", el("\"Table\" property in context panel")));
      await session.step(131, "Then matrix plot viewer should be bound to table \"two-numbers\"", () => boundTable(page, el("matrix plot viewer"), "two-numbers"));
      await session.step(132, "And the \"x columns\" reading of matrix plot viewer should be 2", () => readingIs(page, "x columns", el("matrix plot viewer"), 2));
      await session.step(133, "And the \"cells drawn\" reading of matrix plot viewer should be 4", () => readingIs(page, "cells drawn", el("matrix plot viewer"), 4));
      await session.step(134, "And the \"blank cells\" reading of matrix plot viewer should be 0", () => readingIs(page, "blank cells", el("matrix plot viewer"), 0));
      await session.step(135, "And the \"x scroll shown\" reading of matrix plot viewer should be \"false\"", () => readingReads(page, "x scroll shown", el("matrix plot viewer"), "false"));
      await session.step(136, "And the \"y scroll shown\" reading of matrix plot viewer should be \"false\"", () => readingReads(page, "y scroll shown", el("matrix plot viewer"), "false"));
      await session.step(137, "And matrix plot viewer should not have a \"x scroll slider\" area", () => hasNoArea(page, el("matrix plot viewer"), "x scroll slider"));
      await session.step(138, "And matrix plot viewer should not have a \"y scroll slider\" area", () => hasNoArea(page, el("matrix plot viewer"), "y scroll slider"));
      await session.step(139, "And matrix plot viewer should report no error", () => reportsNoError(page, el("matrix plot viewer")));
      await session.step(140, "When user selects \"demog-1000\" in \"Table\" property in context panel", () => selectIn(page, "demog-1000", el("\"Table\" property in context panel")));
      await session.step(141, "Then matrix plot viewer should be bound to table \"demog-1000\"", () => boundTable(page, el("matrix plot viewer"), "demog-1000"));
      await session.step(142, "And the \"x columns\" reading of matrix plot viewer should be as remembered", () => readingAsRemembered(page, "x columns", el("matrix plot viewer")));
      await session.step(143, "And the \"cells drawn\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells drawn", el("matrix plot viewer"), 16));
      await session.step(144, "And matrix plot viewer should have a \"x scroll slider\" area", () => hasArea(page, el("matrix plot viewer"), "x scroll slider"));
      await session.step(145, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
