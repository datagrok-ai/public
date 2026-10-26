/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/correlation-plot/correlation-plot-colors-and-chrome.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.correlation-plot]
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
import {shouldBe, shouldContainText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCalculated, removeColumn, selectFirstRows, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, areaColor, areaNotColor, areaRepainted, clickArea, hoverArea, noErrors, oneTooltip, painted, pickFromAreaContextMenu, pointerAway, propertyShouldBe, readingBetween, readingIs, readingReads, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Correlation plot — cell colours, the cell tooltip switch, the column labels and a selection", () => {
  const session = feature(test, "features/viewers/correlation-plot/correlation-plot-colors-and-chrome.feature", import.meta.url);
  test("Correlation plot — cell colours, the cell tooltip switch, the column labels and a selection", {tag: ["@journey", "@viewers", "@realizes:viewers.correlation-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a correlation plot viewer", () => addViewer(page, "correlation plot"));
    await session.step(18, "Then the \"cells\" reading of correlation plot viewer should be 16", () => readingIs(page, "cells", el("correlation plot viewer"), 16));
    await session.step(19, "And the \"rows shown\" reading of correlation plot viewer should be 1000", () => readingIs(page, "rows shown", el("correlation plot viewer"), 1000));
    await session.step(20, "And correlation plot viewer should be painted", () => painted(page, el("correlation plot viewer")));
    await run.scenario("A negative coefficient is blue, a positive one red, and one near zero almost white", async () => {
      await session.step(23, "Then the \"correlation of HEIGHT and AGE\" reading of correlation plot viewer should be between -0.2349 and -0.2348", () => readingBetween(page, "correlation of HEIGHT and AGE", el("correlation plot viewer"), -0.2349, -0.2348));
      await session.step(24, "And the \"correlation of WEIGHT and HEIGHT\" reading of correlation plot viewer should be between 0.4124 and 0.4125", () => readingBetween(page, "correlation of WEIGHT and HEIGHT", el("correlation plot viewer"), 0.4124, 0.4125));
      await session.step(25, "And the \"color of cell HEIGHT x AGE\" reading of correlation plot viewer should be \"#c3c3ff\"", () => readingReads(page, "color of cell HEIGHT x AGE", el("correlation plot viewer"), "#c3c3ff"));
      await session.step(26, "And the \"color of cell WEIGHT x HEIGHT\" reading of correlation plot viewer should be \"#ff9595\"", () => readingReads(page, "color of cell WEIGHT x HEIGHT", el("correlation plot viewer"), "#ff9595"));
      await session.step(27, "And the \"color of cell STARTED x AGE\" reading of correlation plot viewer should be \"#fcfcff\"", () => readingReads(page, "color of cell STARTED x AGE", el("correlation plot viewer"), "#fcfcff"));
      await session.step(28, "And the \"cell HEIGHT x AGE\" area of correlation plot viewer should contain the color \"#C3C3FF\"", () => areaColor(page, "cell HEIGHT x AGE", el("correlation plot viewer"), "#C3C3FF"));
      await session.step(29, "And the \"cell HEIGHT x AGE\" area of correlation plot viewer should not contain the color \"#FF9595\"", () => areaNotColor(page, "cell HEIGHT x AGE", el("correlation plot viewer"), "#FF9595"));
      await session.step(30, "And the \"cell WEIGHT x HEIGHT\" area of correlation plot viewer should contain the color \"#FF9595\"", () => areaColor(page, "cell WEIGHT x HEIGHT", el("correlation plot viewer"), "#FF9595"));
      await session.step(31, "And the \"cell WEIGHT x HEIGHT\" area of correlation plot viewer should not contain the color \"#C3C3FF\"", () => areaNotColor(page, "cell WEIGHT x HEIGHT", el("correlation plot viewer"), "#C3C3FF"));
      await session.step(32, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Tooltip > Visible in the cell menu turns Show Tooltip off and back on, and the tooltip returns", async () => {
      await session.step(38, "When user hovers over the \"cell HEIGHT x AGE\" area of correlation plot viewer", () => hoverArea(page, "cell HEIGHT x AGE", el("correlation plot viewer")));
      await session.step(39, "Then tooltip should contain the text \"Pearson R: -0.235\"", () => shouldContainText(page, el("tooltip"), "Pearson R: -0.235"));
      await session.step(40, "When user moves the pointer away from correlation plot viewer", () => pointerAway(page, el("correlation plot viewer")));
      await session.step(41, "And user picks \"Tooltip > Visible\" from the context menu of the \"cell HEIGHT x AGE\" area of correlation plot viewer", () => pickFromAreaContextMenu(page, "Tooltip > Visible", "cell HEIGHT x AGE", el("correlation plot viewer")));
      await session.step(42, "Then \"showTooltip\" property of correlation plot viewer should be \"false\"", () => propertyShouldBe(page, "showTooltip", el("correlation plot viewer"), "false"));
      await session.step(43, "When user clicks on the \"cell STARTED x AGE\" area of correlation plot viewer", () => clickArea(page, "cell STARTED x AGE", el("correlation plot viewer")));
      await session.step(44, "And user hovers over the \"cell WEIGHT x HEIGHT\" area of correlation plot viewer", () => hoverArea(page, "cell WEIGHT x HEIGHT", el("correlation plot viewer")));
      await session.step(45, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(46, "When user moves the pointer away from correlation plot viewer", () => pointerAway(page, el("correlation plot viewer")));
      await session.step(47, "And user picks \"Tooltip > Visible\" from the context menu of the \"cell HEIGHT x AGE\" area of correlation plot viewer", () => pickFromAreaContextMenu(page, "Tooltip > Visible", "cell HEIGHT x AGE", el("correlation plot viewer")));
      await session.step(48, "Then \"showTooltip\" property of correlation plot viewer should be \"true\"", () => propertyShouldBe(page, "showTooltip", el("correlation plot viewer"), "true"));
      await session.step(49, "When user clicks on the \"cell STARTED x AGE\" area of correlation plot viewer", () => clickArea(page, "cell STARTED x AGE", el("correlation plot viewer")));
      await session.step(50, "And user hovers over the \"cell WEIGHT x HEIGHT\" area of correlation plot viewer", () => hoverArea(page, "cell WEIGHT x HEIGHT", el("correlation plot viewer")));
      await session.step(51, "Then tooltip should contain the text \"Pearson R: 0.412\"", () => shouldContainText(page, el("tooltip"), "Pearson R: 0.412"));
      await session.step(52, "When user moves the pointer away from correlation plot viewer", () => pointerAway(page, el("correlation plot viewer")));
      await session.step(53, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column with no variance hovers as N/A instead of a number (GROK-12586)", async () => {
      await session.step(56, "Given user adds a calculated column \"FLAT\" with formula \"0\"", () => addCalculated(page, "FLAT", "0"));
      await session.step(57, "When user sets properties of correlation plot viewer:", () => setProperties(page, el("correlation plot viewer"), [["xColumnNames","AGE, HEIGHT, WEIGHT, STARTED, FLAT"],["yColumnNames","AGE, HEIGHT, WEIGHT, STARTED, FLAT"]]), [["xColumnNames","AGE, HEIGHT, WEIGHT, STARTED, FLAT"],["yColumnNames","AGE, HEIGHT, WEIGHT, STARTED, FLAT"]]);
      await session.step(60, "Then the \"cells\" reading of correlation plot viewer should be 25", () => readingIs(page, "cells", el("correlation plot viewer"), 25));
      await session.step(61, "When user hovers over the \"cell FLAT x AGE\" area of correlation plot viewer", () => hoverArea(page, "cell FLAT x AGE", el("correlation plot viewer")));
      await session.step(62, "Then exactly one tooltip should be shown", () => oneTooltip(page));
      await session.step(63, "And tooltip should contain the text \"N/A\"", () => shouldContainText(page, el("tooltip"), "N/A"));
      await session.step(64, "And tooltip should not contain the text \"NaN\"", () => shouldNotContainText(page, el("tooltip"), "NaN"));
      await session.step(65, "When user moves the pointer away from correlation plot viewer", () => pointerAway(page, el("correlation plot viewer")));
      await session.step(66, "And user sets properties of correlation plot viewer:", () => setProperties(page, el("correlation plot viewer"), [["xColumnNames","AGE, HEIGHT, WEIGHT, STARTED"],["yColumnNames","AGE, HEIGHT, WEIGHT, STARTED"]]), [["xColumnNames","AGE, HEIGHT, WEIGHT, STARTED"],["yColumnNames","AGE, HEIGHT, WEIGHT, STARTED"]]);
      await session.step(69, "And user removes \"FLAT\" column", () => removeColumn(page, "FLAT"));
      await session.step(70, "Then the \"cells\" reading of correlation plot viewer should be 16", () => readingIs(page, "cells", el("correlation plot viewer"), 16));
      await session.step(71, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Col Header Font redraws the column labels and moves no coefficient (GROK-16818)", async () => {
      await session.step(74, "When user sets \"colHeaderFont\" property of correlation plot viewer to 'normal normal 18px \"Roboto\"'", () => setProperty(page, "colHeaderFont", el("correlation plot viewer"), "normal normal 18px \"Roboto\""));
      await session.step(75, "Then the \"header HEIGHT\" area of correlation plot viewer should have repainted", () => areaRepainted(page, "header HEIGHT", el("correlation plot viewer")));
      await session.step(76, "And the \"text of cell HEIGHT x AGE\" reading of correlation plot viewer should be \"-0.23\"", () => readingReads(page, "text of cell HEIGHT x AGE", el("correlation plot viewer"), "-0.23"));
      await session.step(77, "And the \"cells\" reading of correlation plot viewer should be 16", () => readingIs(page, "cells", el("correlation plot viewer"), 16));
      await session.step(78, "When user sets \"colHeaderFont\" property of correlation plot viewer to 'normal normal 13px \"Roboto\"'", () => setProperty(page, "colHeaderFont", el("correlation plot viewer"), "normal normal 13px \"Roboto\""));
      await session.step(79, "Then the \"header HEIGHT\" area of correlation plot viewer should have repainted", () => areaRepainted(page, "header HEIGHT", el("correlation plot viewer")));
      await session.step(80, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A selection is drawn into the diagonal histograms", async () => {
      await session.step(83, "Then the \"cell AGE x AGE\" area of correlation plot viewer should not contain the color \"#FF8C00\"", () => areaNotColor(page, "cell AGE x AGE", el("correlation plot viewer"), "#FF8C00"));
      await session.step(84, "When user selects the first 20 rows", () => selectFirstRows(page, 20));
      await session.step(85, "Then 20 rows should be selected", () => selectedRowCount(page, 20));
      await session.step(86, "And the \"cell AGE x AGE\" area of correlation plot viewer should have repainted", () => areaRepainted(page, "cell AGE x AGE", el("correlation plot viewer")));
      await session.step(87, "And the \"cell AGE x AGE\" area of correlation plot viewer should contain the color \"#FF8C00\"", () => areaColor(page, "cell AGE x AGE", el("correlation plot viewer"), "#FF8C00"));
      await session.step(88, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
