/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/correlation-plot/correlation-plot-scope.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.correlation-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {correlationMatches} from '../../../bindings/correlation-plot.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCategoricalFilter, clearSelection, filterPasses, selectFirstRows, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, openDataset, openProject, saveAsProject, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, boundTable, loadLayout, noErrors, painted, readingBetween, readingDiffers, readingIs, readingReads, repainted, saveLayout, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {readingContains, readingNotContains} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Correlation plot — which rows it correlates, which table, and what survives a round-trip", () => {
  const session = feature(test, "features/viewers/correlation-plot/correlation-plot-scope.feature", import.meta.url);
  test("Correlation plot — which rows it correlates, which table, and what survives a round-trip", {tag: ["@journey", "@viewers", "@realizes:viewers.correlation-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a correlation plot viewer", () => addViewer(page, "correlation plot"));
    await session.step(21, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
    await session.step(22, "And the \"rows shown\" reading of correlation plot viewer should be 1000", () => readingIs(page, "rows shown", el("correlation plot viewer"), 1000));
    await session.step(23, "And the \"correlation of HEIGHT and AGE\" reading of correlation plot viewer should be between -0.2349 and -0.2348", () => readingBetween(page, "correlation of HEIGHT and AGE", el("correlation plot viewer"), -0.2349, -0.2348));
    await session.step(24, "And correlation plot viewer should be painted", () => painted(page, el("correlation plot viewer")));
    await run.scenario("The viewer's own filter narrows what it correlates and leaves the table alone", async () => {
      await session.step(27, "When user sets \"filter\" property of correlation plot viewer to \"${AGE} > 40\"", () => setProperty(page, "filter", el("correlation plot viewer"), "${AGE} > 40"));
      await session.step(28, "Then the \"rows shown\" reading of correlation plot viewer should be 635", () => readingIs(page, "rows shown", el("correlation plot viewer"), 635));
      await session.step(29, "And 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(30, "And the \"correlation of HEIGHT and AGE\" reading of correlation plot viewer should be between -0.2040 and -0.2039", () => readingBetween(page, "correlation of HEIGHT and AGE", el("correlation plot viewer"), -0.204, -0.2039));
      await session.step(31, "And the \"text of cell HEIGHT x AGE\" reading of correlation plot viewer should be \"-0.20\"", () => readingReads(page, "text of cell HEIGHT x AGE", el("correlation plot viewer"), "-0.20"));
      await session.step(32, "And correlation plot viewer should have repainted", () => repainted(page, el("correlation plot viewer")));
      await session.step(33, "When user sets \"filter\" property of correlation plot viewer to \"\"", () => setProperty(page, "filter", el("correlation plot viewer"), ""));
      await session.step(34, "Then the \"rows shown\" reading of correlation plot viewer should be 1000", () => readingIs(page, "rows shown", el("correlation plot viewer"), 1000));
      await session.step(35, "And the \"correlation of HEIGHT and AGE\" reading of correlation plot viewer should be between -0.2349 and -0.2348", () => readingBetween(page, "correlation of HEIGHT and AGE", el("correlation plot viewer"), -0.2349, -0.2348));
      await session.step(36, "And the \"text of cell HEIGHT x AGE\" reading of correlation plot viewer should be \"-0.23\"", () => readingReads(page, "text of cell HEIGHT x AGE", el("correlation plot viewer"), "-0.23"));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A filter on the table moves the coefficients with it", async () => {
      await session.step(40, "When user adds a categorical filter on \"SEX\" keeping \"M\"", () => addCategoricalFilter(page, "SEX", "M"));
      await session.step(41, "Then 447 rows should pass the filter", () => filterPasses(page, 447));
      await session.step(42, "And the \"rows shown\" reading of correlation plot viewer should be 447", () => readingIs(page, "rows shown", el("correlation plot viewer"), 447));
      await session.step(43, "And the \"correlation of HEIGHT and AGE\" reading of correlation plot viewer should differ from before", () => readingDiffers(page, "correlation of HEIGHT and AGE", el("correlation plot viewer")));
      await session.step(44, "And the correlation of \"HEIGHT\" and \"AGE\" of correlation plot viewer should match the Pearson coefficient of the table", () => correlationMatches(page, "HEIGHT", "AGE", el("correlation plot viewer"), "Pearson"));
      await session.step(45, "And correlation plot viewer should have repainted", () => repainted(page, el("correlation plot viewer")));
      await session.step(46, "When user hovers over \"SEX\" filter card", () => hoverOver(page, el("\"SEX\" filter card")));
      await session.step(47, "And user clicks on close of \"SEX\" filter card", () => clickOn(page, el("close of \"SEX\" filter card")));
      await session.step(48, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(49, "And the \"rows shown\" reading of correlation plot viewer should be 1000", () => readingIs(page, "rows shown", el("correlation plot viewer"), 1000));
      await session.step(50, "And the \"correlation of HEIGHT and AGE\" reading of correlation plot viewer should be between -0.2349 and -0.2348", () => readingBetween(page, "correlation of HEIGHT and AGE", el("correlation plot viewer"), -0.2349, -0.2348));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Row Source Selected correlates the selected rows and nothing else", async () => {
      await session.step(54, "When user selects the first 20 rows", () => selectFirstRows(page, 20));
      await session.step(55, "Then 20 rows should be selected", () => selectedRowCount(page, 20));
      await session.step(56, "When user sets \"rowSource\" property of correlation plot viewer to \"Selected\"", () => setProperty(page, "rowSource", el("correlation plot viewer"), "Selected"));
      await session.step(57, "Then the \"rows shown\" reading of correlation plot viewer should be 20", () => readingIs(page, "rows shown", el("correlation plot viewer"), 20));
      await session.step(58, "And the correlation of \"HEIGHT\" and \"AGE\" of correlation plot viewer should match the Pearson coefficient of the table", () => correlationMatches(page, "HEIGHT", "AGE", el("correlation plot viewer"), "Pearson"));
      await session.step(59, "And the \"correlation of HEIGHT and AGE\" reading of correlation plot viewer should be between -0.2161 and -0.2160", () => readingBetween(page, "correlation of HEIGHT and AGE", el("correlation plot viewer"), -0.2161, -0.216));
      await session.step(60, "And correlation plot viewer should have repainted", () => repainted(page, el("correlation plot viewer")));
      await session.step(61, "When user sets \"rowSource\" property of correlation plot viewer to \"Filtered\"", () => setProperty(page, "rowSource", el("correlation plot viewer"), "Filtered"));
      await session.step(62, "And user clears the row selection", () => clearSelection(page));
      await session.step(63, "Then the \"rows shown\" reading of correlation plot viewer should be 1000", () => readingIs(page, "rows shown", el("correlation plot viewer"), 1000));
      await session.step(64, "And the \"correlation of HEIGHT and AGE\" reading of correlation plot viewer should be between -0.2349 and -0.2348", () => readingBetween(page, "correlation of HEIGHT and AGE", el("correlation plot viewer"), -0.2349, -0.2348));
      await session.step(65, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Bound to another table the matrix becomes that table's numerical columns", async () => {
      await session.step(68, "Given user opens spgi dataset", () => openDataset(page, ds("spgi")));
      await session.step(69, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
      await session.step(70, "When user sets \"table\" property of correlation plot viewer to \"spgi-100\"", () => setProperty(page, "table", el("correlation plot viewer"), "spgi-100"));
      await session.step(71, "Then correlation plot viewer should be bound to table \"spgi-100\"", () => boundTable(page, el("correlation plot viewer"), "spgi-100"));
      await session.step(72, "And the \"rows shown\" reading of correlation plot viewer should be 100", () => readingIs(page, "rows shown", el("correlation plot viewer"), 100));
      await session.step(73, "And the \"cells\" reading of correlation plot viewer should be 529", () => readingIs(page, "cells", el("correlation plot viewer"), 529));
      await session.step(74, "And the \"numerical columns\" reading of correlation plot viewer should contain \"TPSA\"", () => readingContains(page, "numerical columns", el("correlation plot viewer"), "TPSA"));
      await session.step(75, "And the \"numerical columns\" reading of correlation plot viewer should not contain \"HEIGHT\"", () => readingNotContains(page, "numerical columns", el("correlation plot viewer"), "HEIGHT"));
      await session.step(76, "And the \"error\" reading of correlation plot viewer should be \"\"", () => readingReads(page, "error", el("correlation plot viewer"), ""));
      await session.step(77, "And correlation plot viewer should be painted", () => painted(page, el("correlation plot viewer")));
      await session.step(78, "When user sets \"table\" property of correlation plot viewer to \"demog-1000\"", () => setProperty(page, "table", el("correlation plot viewer"), "demog-1000"));
      await session.step(79, "Then correlation plot viewer should be bound to table \"demog-1000\"", () => boundTable(page, el("correlation plot viewer"), "demog-1000"));
      await session.step(80, "And the \"cells\" reading of correlation plot viewer should be 16", () => readingIs(page, "cells", el("correlation plot viewer"), 16));
      await session.step(81, "And the \"x columns\" reading of correlation plot viewer should be \"AGE, HEIGHT, WEIGHT, STARTED\"", () => readingReads(page, "x columns", el("correlation plot viewer"), "AGE, HEIGHT, WEIGHT, STARTED"));
      await session.step(82, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A saved layout brings the configured matrix back", async () => {
      await session.step(85, "When user sets properties of correlation plot viewer:", () => setProperties(page, el("correlation plot viewer"), [["correlationType","Spearman"],["showPearsonR","false"],["xColumnNames","AGE, HEIGHT, WEIGHT"],["yColumnNames","AGE, HEIGHT"]]));
      await session.step(90, "Then the \"cells\" reading of correlation plot viewer should be 6", () => readingIs(page, "cells", el("correlation plot viewer"), 6));
      await session.step(91, "And user saves the layout of the current table view", () => saveLayout(page));
      await session.step(92, "When user sets properties of correlation plot viewer:", () => setProperties(page, el("correlation plot viewer"), [["correlationType","Pearson"],["showPearsonR","true"],["xColumnNames","AGE, HEIGHT, WEIGHT, STARTED"],["yColumnNames","AGE, HEIGHT, WEIGHT, STARTED"]]));
      await session.step(97, "Then the \"cells\" reading of correlation plot viewer should be 16", () => readingIs(page, "cells", el("correlation plot viewer"), 16));
      await session.step(98, "When user loads the saved layout", () => loadLayout(page));
      await session.step(99, "Then correlation plot viewer should be visible", () => shouldBe(page, el("correlation plot viewer"), "visible"));
      await session.step(100, "And the \"cells\" reading of correlation plot viewer should be 6", () => readingIs(page, "cells", el("correlation plot viewer"), 6));
      await session.step(101, "And the \"correlation type\" reading of correlation plot viewer should be \"Spearman\"", () => readingReads(page, "correlation type", el("correlation plot viewer"), "Spearman"));
      await session.step(102, "And the \"show pearson r\" reading of correlation plot viewer should be \"false\"", () => readingReads(page, "show pearson r", el("correlation plot viewer"), "false"));
      await session.step(103, "And the \"x columns\" reading of correlation plot viewer should be \"AGE, HEIGHT, WEIGHT\"", () => readingReads(page, "x columns", el("correlation plot viewer"), "AGE, HEIGHT, WEIGHT"));
      await session.step(104, "And the \"y columns\" reading of correlation plot viewer should be \"AGE, HEIGHT\"", () => readingReads(page, "y columns", el("correlation plot viewer"), "AGE, HEIGHT"));
      await session.step(105, "And the correlation of \"WEIGHT\" and \"HEIGHT\" of correlation plot viewer should match the Spearman coefficient of the table", () => correlationMatches(page, "WEIGHT", "HEIGHT", el("correlation plot viewer"), "Spearman"));
      await session.step(106, "When user sets properties of correlation plot viewer:", () => setProperties(page, el("correlation plot viewer"), [["correlationType","Pearson"],["showPearsonR","true"],["xColumnNames","AGE, HEIGHT, WEIGHT, STARTED"],["yColumnNames","AGE, HEIGHT, WEIGHT, STARTED"]]));
      await session.step(111, "Then the \"cells\" reading of correlation plot viewer should be 16", () => readingIs(page, "cells", el("correlation plot viewer"), 16));
      await session.step(112, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A project round-trip brings it back too", async () => {
      await session.step(115, "When user sets properties of correlation plot viewer:", () => setProperties(page, el("correlation plot viewer"), [["correlationType","Spearman"],["xColumnNames","AGE, HEIGHT, WEIGHT"],["yColumnNames","AGE, HEIGHT"]]));
      await session.step(119, "Then the \"cells\" reading of correlation plot viewer should be 6", () => readingIs(page, "cells", el("correlation plot viewer"), 6));
      await session.step(120, "When user saves the current view as project \"bdd correlation matrix\"", () => saveAsProject(page, "bdd correlation matrix"));
      await session.step(121, "And user closes all views", () => closeAllViews(page));
      await session.step(122, "And user opens the \"bdd correlation matrix\" project", () => openProject(page, "bdd correlation matrix"));
      await session.step(123, "Then correlation plot viewer should be visible", () => shouldBe(page, el("correlation plot viewer"), "visible"));
      await session.step(124, "And the \"cells\" reading of correlation plot viewer should be 6", () => readingIs(page, "cells", el("correlation plot viewer"), 6));
      await session.step(125, "And the \"correlation type\" reading of correlation plot viewer should be \"Spearman\"", () => readingReads(page, "correlation type", el("correlation plot viewer"), "Spearman"));
      await session.step(126, "And the \"x columns\" reading of correlation plot viewer should be \"AGE, HEIGHT, WEIGHT\"", () => readingReads(page, "x columns", el("correlation plot viewer"), "AGE, HEIGHT, WEIGHT"));
      await session.step(127, "And the \"rows shown\" reading of correlation plot viewer should be 1000", () => readingIs(page, "rows shown", el("correlation plot viewer"), 1000));
      await session.step(128, "And the correlation of \"WEIGHT\" and \"HEIGHT\" of correlation plot viewer should match the Spearman coefficient of the table", () => correlationMatches(page, "WEIGHT", "HEIGHT", el("correlation plot viewer"), "Spearman"));
      await session.step(129, "When user sets properties of correlation plot viewer:", () => setProperties(page, el("correlation plot viewer"), [["correlationType","Pearson"],["xColumnNames","AGE, HEIGHT, WEIGHT, STARTED"],["yColumnNames","AGE, HEIGHT, WEIGHT, STARTED"]]));
      await session.step(133, "Then the \"cells\" reading of correlation plot viewer should be 16", () => readingIs(page, "cells", el("correlation plot viewer"), 16));
      await session.step(134, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A bigger cell font makes the rows taller and moves no coefficient", async () => {
      await session.step(137, "Then the \"row height\" reading of correlation plot viewer should be 20", () => readingIs(page, "row height", el("correlation plot viewer"), 20));
      await session.step(138, "When user sets \"defaultCellFont\" property of correlation plot viewer to 'normal normal 20px \"Roboto\"'", () => setProperty(page, "defaultCellFont", el("correlation plot viewer"), "normal normal 20px \"Roboto\""));
      await session.step(139, "Then the \"row height\" reading of correlation plot viewer should be between 27.9 and 28.1", () => readingBetween(page, "row height", el("correlation plot viewer"), 27.9, 28.1));
      await session.step(140, "And the \"correlation of HEIGHT and AGE\" reading of correlation plot viewer should be between -0.2349 and -0.2348", () => readingBetween(page, "correlation of HEIGHT and AGE", el("correlation plot viewer"), -0.2349, -0.2348));
      await session.step(141, "And the \"text of cell HEIGHT x AGE\" reading of correlation plot viewer should be \"-0.23\"", () => readingReads(page, "text of cell HEIGHT x AGE", el("correlation plot viewer"), "-0.23"));
      await session.step(142, "And correlation plot viewer should have repainted", () => repainted(page, el("correlation plot viewer")));
      await session.step(143, "When user sets \"defaultCellFont\" property of correlation plot viewer to 'normal normal 13px \"Roboto\"'", () => setProperty(page, "defaultCellFont", el("correlation plot viewer"), "normal normal 13px \"Roboto\""));
      await session.step(144, "Then the \"row height\" reading of correlation plot viewer should be between 18.1 and 18.3", () => readingBetween(page, "row height", el("correlation plot viewer"), 18.1, 18.3));
      await session.step(145, "And the \"correlation of HEIGHT and AGE\" reading of correlation plot viewer should be between -0.2349 and -0.2348", () => readingBetween(page, "correlation of HEIGHT and AGE", el("correlation plot viewer"), -0.2349, -0.2348));
      await session.step(146, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
