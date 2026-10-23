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
import {shouldBe, shouldContainText, shouldHaveText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCalculated, removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, areaColor, areaNotColor, areaRepainted, clickArea, hoverArea, noErrors, oneTooltip, painted, pickFromAreaContextMenu, pointerAway, propertyShouldBe, readingBetween, readingIs, readingReads, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {descriptionAbove, descriptionBelow} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Correlation plot — cell colours, the cell tooltip switch, and the chrome around the matrix", () => {
  const session = feature(test, "features/viewers/correlation-plot/correlation-plot-colors-and-chrome.feature", import.meta.url);
  test("Correlation plot — cell colours, the cell tooltip switch, and the chrome around the matrix", {tag: ["@journey", "@viewers", "@realizes:viewers.correlation-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user adds a correlation plot viewer", () => addViewer(page, "correlation plot"));
    await session.step(20, "Then the \"cells\" reading of correlation plot viewer should be 16", () => readingIs(page, "cells", el("correlation plot viewer"), 16));
    await session.step(21, "And the \"rows shown\" reading of correlation plot viewer should be 1000", () => readingIs(page, "rows shown", el("correlation plot viewer"), 1000));
    await session.step(22, "And correlation plot viewer should be painted", () => painted(page, el("correlation plot viewer")));
    await run.scenario("A negative coefficient is blue, a positive one red, and one near zero almost white", async () => {
      await session.step(25, "Then the \"correlation of HEIGHT and AGE\" reading of correlation plot viewer should be between -0.2349 and -0.2348", () => readingBetween(page, "correlation of HEIGHT and AGE", el("correlation plot viewer"), -0.2349, -0.2348));
      await session.step(26, "And the \"correlation of WEIGHT and HEIGHT\" reading of correlation plot viewer should be between 0.4124 and 0.4125", () => readingBetween(page, "correlation of WEIGHT and HEIGHT", el("correlation plot viewer"), 0.4124, 0.4125));
      await session.step(27, "And the \"color of cell HEIGHT x AGE\" reading of correlation plot viewer should be \"#c3c3ff\"", () => readingReads(page, "color of cell HEIGHT x AGE", el("correlation plot viewer"), "#c3c3ff"));
      await session.step(28, "And the \"color of cell WEIGHT x HEIGHT\" reading of correlation plot viewer should be \"#ff9595\"", () => readingReads(page, "color of cell WEIGHT x HEIGHT", el("correlation plot viewer"), "#ff9595"));
      await session.step(29, "And the \"color of cell STARTED x AGE\" reading of correlation plot viewer should be \"#fcfcff\"", () => readingReads(page, "color of cell STARTED x AGE", el("correlation plot viewer"), "#fcfcff"));
      await session.step(30, "And the \"cell HEIGHT x AGE\" area of correlation plot viewer should contain the color \"#C3C3FF\"", () => areaColor(page, "cell HEIGHT x AGE", el("correlation plot viewer"), "#C3C3FF"));
      await session.step(31, "And the \"cell HEIGHT x AGE\" area of correlation plot viewer should not contain the color \"#FF9595\"", () => areaNotColor(page, "cell HEIGHT x AGE", el("correlation plot viewer"), "#FF9595"));
      await session.step(32, "And the \"cell WEIGHT x HEIGHT\" area of correlation plot viewer should contain the color \"#FF9595\"", () => areaColor(page, "cell WEIGHT x HEIGHT", el("correlation plot viewer"), "#FF9595"));
      await session.step(33, "And the \"cell WEIGHT x HEIGHT\" area of correlation plot viewer should not contain the color \"#C3C3FF\"", () => areaNotColor(page, "cell WEIGHT x HEIGHT", el("correlation plot viewer"), "#C3C3FF"));
      await session.step(34, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The title shows and the description sits above the matrix, then below it, then nowhere", async () => {
      await session.step(37, "When user sets properties of correlation plot viewer:", () => setProperties(page, el("correlation plot viewer"), [["showTitle","true"],["title","Pairwise r"],["description","Every numerical pair"]]), [["showTitle","true"],["title","Pairwise r"],["description","Every numerical pair"]]);
      await session.step(41, "Then title of correlation plot viewer should have text \"Pairwise r\"", () => shouldHaveText(page, el("title of correlation plot viewer"), "Pairwise r"));
      await session.step(42, "And description of correlation plot viewer should be visible", () => shouldBe(page, el("description of correlation plot viewer"), "visible"));
      await session.step(43, "And description of correlation plot viewer should have text \"Every numerical pair\"", () => shouldHaveText(page, el("description of correlation plot viewer"), "Every numerical pair"));
      await session.step(44, "And the description of correlation plot viewer should be above its content", () => descriptionAbove(page, el("correlation plot viewer")));
      await session.step(45, "When user sets \"descriptionPosition\" property of correlation plot viewer to \"Bottom\"", () => setProperty(page, "descriptionPosition", el("correlation plot viewer"), "Bottom"));
      await session.step(46, "Then the description of correlation plot viewer should be below its content", () => descriptionBelow(page, el("correlation plot viewer")));
      await session.step(47, "When user sets \"descriptionVisibilityMode\" property of correlation plot viewer to \"Never\"", () => setProperty(page, "descriptionVisibilityMode", el("correlation plot viewer"), "Never"));
      await session.step(48, "Then description of correlation plot viewer should be hidden", () => shouldBe(page, el("description of correlation plot viewer"), "hidden"));
      await session.step(49, "When user sets properties of correlation plot viewer:", () => setProperties(page, el("correlation plot viewer"), [["descriptionVisibilityMode","Auto"],["descriptionPosition","Top"],["description",""],["title",""],["showTitle","false"]]), [["descriptionVisibilityMode","Auto"],["descriptionPosition","Top"],["description",""],["title",""],["showTitle","false"]]);
      await session.step(55, "Then the \"cells\" reading of correlation plot viewer should be 16", () => readingIs(page, "cells", el("correlation plot viewer"), 16));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Tooltip > Visible in the cell menu switches the cell tooltip off and back on", async () => {
      await session.step(62, "When user hovers over the \"cell HEIGHT x AGE\" area of correlation plot viewer", () => hoverArea(page, "cell HEIGHT x AGE", el("correlation plot viewer")));
      await session.step(63, "Then tooltip should contain the text \"Pearson R: -0.235\"", () => shouldContainText(page, el("tooltip"), "Pearson R: -0.235"));
      await session.step(64, "When user moves the pointer away from correlation plot viewer", () => pointerAway(page, el("correlation plot viewer")));
      await session.step(65, "And user picks \"Tooltip > Visible\" from the context menu of the \"cell HEIGHT x AGE\" area of correlation plot viewer", () => pickFromAreaContextMenu(page, "Tooltip > Visible", "cell HEIGHT x AGE", el("correlation plot viewer")));
      await session.step(66, "Then \"showTooltip\" property of correlation plot viewer should be \"false\"", () => propertyShouldBe(page, "showTooltip", el("correlation plot viewer"), "false"));
      await session.step(67, "When user clicks on the \"cell STARTED x AGE\" area of correlation plot viewer", () => clickArea(page, "cell STARTED x AGE", el("correlation plot viewer")));
      await session.step(68, "And user hovers over the \"cell WEIGHT x HEIGHT\" area of correlation plot viewer", () => hoverArea(page, "cell WEIGHT x HEIGHT", el("correlation plot viewer")));
      await session.step(69, "Then tooltip should not contain the text \"Pearson R\"", () => shouldNotContainText(page, el("tooltip"), "Pearson R"));
      await session.step(70, "When user moves the pointer away from correlation plot viewer", () => pointerAway(page, el("correlation plot viewer")));
      await session.step(71, "And user picks \"Tooltip > Visible\" from the context menu of the \"cell HEIGHT x AGE\" area of correlation plot viewer", () => pickFromAreaContextMenu(page, "Tooltip > Visible", "cell HEIGHT x AGE", el("correlation plot viewer")));
      await session.step(72, "Then \"showTooltip\" property of correlation plot viewer should be \"true\"", () => propertyShouldBe(page, "showTooltip", el("correlation plot viewer"), "true"));
      await session.step(73, "When user clicks on the \"cell STARTED x AGE\" area of correlation plot viewer", () => clickArea(page, "cell STARTED x AGE", el("correlation plot viewer")));
      await session.step(74, "And user hovers over the \"cell WEIGHT x HEIGHT\" area of correlation plot viewer", () => hoverArea(page, "cell WEIGHT x HEIGHT", el("correlation plot viewer")));
      await session.step(75, "Then tooltip should contain the text \"Pearson R: 0.412\"", () => shouldContainText(page, el("tooltip"), "Pearson R: 0.412"));
      await session.step(76, "When user moves the pointer away from correlation plot viewer", () => pointerAway(page, el("correlation plot viewer")));
      await session.step(77, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column with no variance hovers as N/A instead of a number (GROK-12586)", async () => {
      await session.step(80, "Given user adds a calculated column \"FLAT\" with formula \"0\"", () => addCalculated(page, "FLAT", "0"));
      await session.step(81, "When user sets properties of correlation plot viewer:", () => setProperties(page, el("correlation plot viewer"), [["xColumnNames","AGE, HEIGHT, WEIGHT, STARTED, FLAT"],["yColumnNames","AGE, HEIGHT, WEIGHT, STARTED, FLAT"]]), [["xColumnNames","AGE, HEIGHT, WEIGHT, STARTED, FLAT"],["yColumnNames","AGE, HEIGHT, WEIGHT, STARTED, FLAT"]]);
      await session.step(84, "Then the \"cells\" reading of correlation plot viewer should be 25", () => readingIs(page, "cells", el("correlation plot viewer"), 25));
      await session.step(85, "When user hovers over the \"cell FLAT x AGE\" area of correlation plot viewer", () => hoverArea(page, "cell FLAT x AGE", el("correlation plot viewer")));
      await session.step(86, "Then exactly one tooltip should be shown", () => oneTooltip(page));
      await session.step(87, "And tooltip should contain the text \"N/A\"", () => shouldContainText(page, el("tooltip"), "N/A"));
      await session.step(88, "And tooltip should not contain the text \"NaN\"", () => shouldNotContainText(page, el("tooltip"), "NaN"));
      await session.step(89, "When user moves the pointer away from correlation plot viewer", () => pointerAway(page, el("correlation plot viewer")));
      await session.step(90, "And user sets properties of correlation plot viewer:", () => setProperties(page, el("correlation plot viewer"), [["xColumnNames","AGE, HEIGHT, WEIGHT, STARTED"],["yColumnNames","AGE, HEIGHT, WEIGHT, STARTED"]]), [["xColumnNames","AGE, HEIGHT, WEIGHT, STARTED"],["yColumnNames","AGE, HEIGHT, WEIGHT, STARTED"]]);
      await session.step(93, "And user removes \"FLAT\" column", () => removeColumn(page, "FLAT"));
      await session.step(94, "Then the \"cells\" reading of correlation plot viewer should be 16", () => readingIs(page, "cells", el("correlation plot viewer"), 16));
      await session.step(95, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Col Header Font redraws the column labels and moves no coefficient (GROK-16818)", async () => {
      await session.step(98, "When user sets \"colHeaderFont\" property of correlation plot viewer to 'normal normal 18px \"Roboto\"'", () => setProperty(page, "colHeaderFont", el("correlation plot viewer"), "normal normal 18px \"Roboto\""));
      await session.step(99, "Then the \"header HEIGHT\" area of correlation plot viewer should have repainted", () => areaRepainted(page, "header HEIGHT", el("correlation plot viewer")));
      await session.step(100, "And the \"text of cell HEIGHT x AGE\" reading of correlation plot viewer should be \"-0.23\"", () => readingReads(page, "text of cell HEIGHT x AGE", el("correlation plot viewer"), "-0.23"));
      await session.step(101, "And the \"cells\" reading of correlation plot viewer should be 16", () => readingIs(page, "cells", el("correlation plot viewer"), 16));
      await session.step(102, "When user sets \"colHeaderFont\" property of correlation plot viewer to 'normal normal 13px \"Roboto\"'", () => setProperty(page, "colHeaderFont", el("correlation plot viewer"), "normal normal 13px \"Roboto\""));
      await session.step(103, "Then the \"header HEIGHT\" area of correlation plot viewer should have repainted", () => areaRepainted(page, "header HEIGHT", el("correlation plot viewer")));
      await session.step(104, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
