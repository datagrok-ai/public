/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/trellis-plot/trellis-plot-categories.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.trellis-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {cellsWideTall} from '../../../bindings/trellis-plot.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCategoricalFilter, filterPasses, filterPassesAll} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, hasArea, hasNoArea, noErrors, readingIs, readingLower, setProperties, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Trellis plot categories, labels and scrolling", () => {
  const session = feature(test, "features/viewers/trellis-plot/trellis-plot-categories.feature", import.meta.url);
  test("Trellis plot categories, labels and scrolling", {tag: ["@journey", "@viewers", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"],["Pack Categories","false"]]));
    await session.step(22, "Then the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
    await run.scenario("A second split column grows the grid to the clamped product", async () => {
      await session.step(25, "Then the \"x categories\" reading of trellis plot viewer should be 2", () => readingIs(page, "x categories", el("trellis plot viewer"), 2));
      await session.step(26, "When user sets \"X Column Names\" property of trellis plot viewer to \"SEX, DIS_POP\"", () => setProperty(page, "X Column Names", el("trellis plot viewer"), "SEX, DIS_POP"));
      await session.step(27, "Then the \"x categories\" reading of trellis plot viewer should be 12", () => readingIs(page, "x categories", el("trellis plot viewer"), 12));
      await session.step(28, "And the cells of trellis plot viewer should be 5 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 5, 4));
      await session.step(29, "And the \"cells\" reading of trellis plot viewer should be 20", () => readingIs(page, "cells", el("trellis plot viewer"), 20));
      await session.step(30, "And the \"cells drawn\" reading of trellis plot viewer should be 18", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 18));
      await session.step(31, "When user sets \"X Column Names\" property of trellis plot viewer to \"SEX\"", () => setProperty(page, "X Column Names", el("trellis plot viewer"), "SEX"));
      await session.step(32, "Then the \"x categories\" reading of trellis plot viewer should be 2", () => readingIs(page, "x categories", el("trellis plot viewer"), 2));
      await session.step(33, "And the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
      await session.step(34, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show X and Y Labels remove and restore the label strips", async () => {
      await session.step(37, "Then the \"x labels shown\" reading of trellis plot viewer should be 2", () => readingIs(page, "x labels shown", el("trellis plot viewer"), 2));
      await session.step(38, "And the \"y labels shown\" reading of trellis plot viewer should be 4", () => readingIs(page, "y labels shown", el("trellis plot viewer"), 4));
      await session.step(39, "And trellis plot viewer should have an \"x label F\" area", () => hasArea(page, el("trellis plot viewer"), "x label F"));
      await session.step(40, "And trellis plot viewer should have a \"y label Caucasian\" area", () => hasArea(page, el("trellis plot viewer"), "y label Caucasian"));
      await session.step(41, "When user sets \"Show X Labels\" property of trellis plot viewer to \"false\"", () => setProperty(page, "Show X Labels", el("trellis plot viewer"), "false"));
      await session.step(42, "Then the \"x labels shown\" reading of trellis plot viewer should be 0", () => readingIs(page, "x labels shown", el("trellis plot viewer"), 0));
      await session.step(43, "And trellis plot viewer should not have an \"x label F\" area", () => hasNoArea(page, el("trellis plot viewer"), "x label F"));
      await session.step(44, "And the \"y labels shown\" reading of trellis plot viewer should be 4", () => readingIs(page, "y labels shown", el("trellis plot viewer"), 4));
      await session.step(45, "When user sets \"Show Y Labels\" property of trellis plot viewer to \"false\"", () => setProperty(page, "Show Y Labels", el("trellis plot viewer"), "false"));
      await session.step(46, "Then the \"y labels shown\" reading of trellis plot viewer should be 0", () => readingIs(page, "y labels shown", el("trellis plot viewer"), 0));
      await session.step(47, "And trellis plot viewer should not have a \"y label Caucasian\" area", () => hasNoArea(page, el("trellis plot viewer"), "y label Caucasian"));
      await session.step(48, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Show X Labels","true"],["Show Y Labels","true"]]));
      await session.step(51, "Then the \"x labels shown\" reading of trellis plot viewer should be 2", () => readingIs(page, "x labels shown", el("trellis plot viewer"), 2));
      await session.step(52, "And the \"y labels shown\" reading of trellis plot viewer should be 4", () => readingIs(page, "y labels shown", el("trellis plot viewer"), 4));
      await session.step(53, "And trellis plot viewer should have an \"x label F\" area", () => hasArea(page, el("trellis plot viewer"), "x label F"));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Label orientation is horizontal, vertical, or one of each", async () => {
      await session.step(57, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["X Labels Orientation","Horz"],["Y Labels Orientation","Horz"]]));
      await session.step(60, "Then the \"x label angle\" reading of trellis plot viewer should be 0", () => readingIs(page, "x label angle", el("trellis plot viewer"), 0));
      await session.step(61, "And the \"y label angle\" reading of trellis plot viewer should be 0", () => readingIs(page, "y label angle", el("trellis plot viewer"), 0));
      await session.step(62, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["X Labels Orientation","Vert"],["Y Labels Orientation","Vert"]]));
      await session.step(65, "Then the \"x label angle\" reading of trellis plot viewer should be -90", () => readingIs(page, "x label angle", el("trellis plot viewer"), -90));
      await session.step(66, "And the \"y label angle\" reading of trellis plot viewer should be -90", () => readingIs(page, "y label angle", el("trellis plot viewer"), -90));
      await session.step(67, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["X Labels Orientation","Auto"],["Y Labels Orientation","Auto"]]));
      await session.step(70, "Then the \"x label angle\" reading of trellis plot viewer should be 0", () => readingIs(page, "x label angle", el("trellis plot viewer"), 0));
      await session.step(71, "And the \"y label angle\" reading of trellis plot viewer should be -90", () => readingIs(page, "y label angle", el("trellis plot viewer"), -90));
      await session.step(72, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Every category fits, so both scroll handles fill their track", async () => {
      await session.step(75, "Then the \"x scroll handle share\" reading of trellis plot viewer should be 1", () => readingIs(page, "x scroll handle share", el("trellis plot viewer"), 1));
      await session.step(76, "And the \"y scroll handle share\" reading of trellis plot viewer should be 1", () => readingIs(page, "y scroll handle share", el("trellis plot viewer"), 1));
      await session.step(77, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An overflowing X axis shortens its handle", async () => {
      await session.step(80, "When user sets \"X Column Names\" property of trellis plot viewer to \"SEX, DIS_POP\"", () => setProperty(page, "X Column Names", el("trellis plot viewer"), "SEX, DIS_POP"));
      await session.step(81, "Then the \"x categories\" reading of trellis plot viewer should be 12", () => readingIs(page, "x categories", el("trellis plot viewer"), 12));
      await session.step(82, "And the cells of trellis plot viewer should be 5 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 5, 4));
      await session.step(83, "And the \"x scroll handle share\" reading of trellis plot viewer should be lower than before", () => readingLower(page, "x scroll handle share", el("trellis plot viewer")));
      await session.step(84, "And the \"y scroll handle share\" reading of trellis plot viewer should be 1", () => readingIs(page, "y scroll handle share", el("trellis plot viewer"), 1));
      await session.step(85, "When user sets \"X Column Names\" property of trellis plot viewer to \"SEX\"", () => setProperty(page, "X Column Names", el("trellis plot viewer"), "SEX"));
      await session.step(86, "Then the \"x scroll handle share\" reading of trellis plot viewer should be 1", () => readingIs(page, "x scroll handle share", el("trellis plot viewer"), 1));
      await session.step(87, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An overflowing Y axis does the same", async () => {
      await session.step(90, "When user sets \"Y Column Names\" property of trellis plot viewer to \"DIS_POP, RACE\"", () => setProperty(page, "Y Column Names", el("trellis plot viewer"), "DIS_POP, RACE"));
      await session.step(91, "Then the \"y categories\" reading of trellis plot viewer should be 24", () => readingIs(page, "y categories", el("trellis plot viewer"), 24));
      await session.step(92, "And the cells of trellis plot viewer should be 2 wide and 5 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 5));
      await session.step(93, "And the \"y scroll handle share\" reading of trellis plot viewer should be lower than before", () => readingLower(page, "y scroll handle share", el("trellis plot viewer")));
      await session.step(94, "And the \"x scroll handle share\" reading of trellis plot viewer should be 1", () => readingIs(page, "x scroll handle share", el("trellis plot viewer"), 1));
      await session.step(95, "When user sets \"Y Column Names\" property of trellis plot viewer to \"RACE\"", () => setProperty(page, "Y Column Names", el("trellis plot viewer"), "RACE"));
      await session.step(96, "Then the \"y scroll handle share\" reading of trellis plot viewer should be 1", () => readingIs(page, "y scroll handle share", el("trellis plot viewer"), 1));
      await session.step(97, "And the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
      await session.step(98, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The plus and minus icons page one category row in and out", async () => {
      await session.step(101, "When user sets \"X Column Names\" property of trellis plot viewer to \"SEX, DIS_POP\"", () => setProperty(page, "X Column Names", el("trellis plot viewer"), "SEX, DIS_POP"));
      await session.step(102, "Then the cells of trellis plot viewer should be 5 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 5, 4));
      await session.step(103, "And x plus icon should be enabled", () => shouldBe(page, el("x plus icon"), "enabled"));
      await session.step(104, "When user clicks on the \"x plus\" area of trellis plot viewer", () => clickArea(page, "x plus", el("trellis plot viewer")));
      await session.step(105, "Then the cells of trellis plot viewer should be 6 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 6, 4));
      await session.step(106, "And the \"cells\" reading of trellis plot viewer should be 24", () => readingIs(page, "cells", el("trellis plot viewer"), 24));
      await session.step(107, "When user clicks on the \"x minus\" area of trellis plot viewer", () => clickArea(page, "x minus", el("trellis plot viewer")));
      await session.step(108, "Then the cells of trellis plot viewer should be 5 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 5, 4));
      await session.step(109, "And the \"cells\" reading of trellis plot viewer should be 20", () => readingIs(page, "cells", el("trellis plot viewer"), 20));
      await session.step(110, "And no errors should have been logged", () => noErrors(page));
      await session.step(111, "When user sets \"X Column Names\" property of trellis plot viewer to \"SEX\"", () => setProperty(page, "X Column Names", el("trellis plot viewer"), "SEX"));
    });
    await run.scenario("The icons go inert at the ends", async () => {
      await session.step(114, "When user sets \"X Column Names\" property of trellis plot viewer to \"SEX, DIS_POP\"", () => setProperty(page, "X Column Names", el("trellis plot viewer"), "SEX, DIS_POP"));
      await session.step(115, "Then the cells of trellis plot viewer should be 5 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 5, 4));
      await session.step(116, "When user clicks on the \"x minus\" area of trellis plot viewer", () => clickArea(page, "x minus", el("trellis plot viewer")));
      await session.step(117, "And user clicks on the \"x minus\" area of trellis plot viewer", () => clickArea(page, "x minus", el("trellis plot viewer")));
      await session.step(118, "And user clicks on the \"x minus\" area of trellis plot viewer", () => clickArea(page, "x minus", el("trellis plot viewer")));
      await session.step(119, "And user clicks on the \"x minus\" area of trellis plot viewer", () => clickArea(page, "x minus", el("trellis plot viewer")));
      await session.step(120, "Then the cells of trellis plot viewer should be 1 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 1, 4));
      await session.step(121, "And x minus icon should be disabled", () => shouldBe(page, el("x minus icon"), "disabled"));
      await session.step(122, "And x plus icon should be enabled", () => shouldBe(page, el("x plus icon"), "enabled"));
      await session.step(123, "When user clicks on the \"x minus\" area of trellis plot viewer", () => clickArea(page, "x minus", el("trellis plot viewer")));
      await session.step(124, "Then the cells of trellis plot viewer should be 1 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 1, 4));
      await session.step(125, "And y plus icon should be disabled", () => shouldBe(page, el("y plus icon"), "disabled"));
      await session.step(126, "When user clicks on the \"y plus\" area of trellis plot viewer", () => clickArea(page, "y plus", el("trellis plot viewer")));
      await session.step(127, "Then the cells of trellis plot viewer should be 1 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 1, 4));
      await session.step(128, "When user sets \"X Column Names\" property of trellis plot viewer to \"SEX\"", () => setProperty(page, "X Column Names", el("trellis plot viewer"), "SEX"));
      await session.step(129, "Then the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
      await session.step(130, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Packing drops the categories a filter leaves empty", async () => {
      await session.step(133, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["X Column Names","RACE"],["Y Column Names","SEX"],["Pack Categories","true"]]));
      await session.step(137, "Then the \"x categories packed\" reading of trellis plot viewer should be 4", () => readingIs(page, "x categories packed", el("trellis plot viewer"), 4));
      await session.step(138, "And the cells of trellis plot viewer should be 4 wide and 2 tall", () => cellsWideTall(page, el("trellis plot viewer"), 4, 2));
      await session.step(139, "When user adds a categorical filter on \"RACE\" keeping \"Caucasian\"", () => addCategoricalFilter(page, "RACE", "Caucasian"));
      await session.step(140, "Then 896 rows should pass the filter", () => filterPasses(page, 896));
      await session.step(141, "And the \"x categories\" reading of trellis plot viewer should be 4", () => readingIs(page, "x categories", el("trellis plot viewer"), 4));
      await session.step(142, "And the \"x categories packed\" reading of trellis plot viewer should be 1", () => readingIs(page, "x categories packed", el("trellis plot viewer"), 1));
      await session.step(143, "And the cells of trellis plot viewer should be 1 wide and 2 tall", () => cellsWideTall(page, el("trellis plot viewer"), 1, 2));
      await session.step(144, "And trellis plot viewer should show 896 rows", () => showsRows(page, el("trellis plot viewer"), 896));
      await session.step(145, "When user sets \"Pack Categories\" property of trellis plot viewer to \"false\"", () => setProperty(page, "Pack Categories", el("trellis plot viewer"), "false"));
      await session.step(146, "Then the \"x categories packed\" reading of trellis plot viewer should be 4", () => readingIs(page, "x categories packed", el("trellis plot viewer"), 4));
      await session.step(147, "And the cells of trellis plot viewer should be 4 wide and 2 tall", () => cellsWideTall(page, el("trellis plot viewer"), 4, 2));
      await session.step(148, "When user sets \"Pack Categories\" property of trellis plot viewer to \"true\"", () => setProperty(page, "Pack Categories", el("trellis plot viewer"), "true"));
      await session.step(149, "Then the cells of trellis plot viewer should be 1 wide and 2 tall", () => cellsWideTall(page, el("trellis plot viewer"), 1, 2));
      await session.step(150, "When user hovers over \"RACE\" filter card", () => hoverOver(page, el("\"RACE\" filter card")));
      await session.step(151, "And user clicks on close of \"RACE\" filter card", () => clickOn(page, el("close of \"RACE\" filter card")));
      await session.step(152, "Then \"RACE\" filter card should be absent", () => shouldBe(page, el("\"RACE\" filter card"), "absent"));
      await session.step(153, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(154, "And the cells of trellis plot viewer should be 4 wide and 2 tall", () => cellsWideTall(page, el("trellis plot viewer"), 4, 2));
      await session.step(155, "And no errors should have been logged", () => noErrors(page));
      await session.step(156, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["X Column Names","SEX"],["Y Column Names","RACE"],["Pack Categories","false"]]));
    });
    run.finish();
  });
});
