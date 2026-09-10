/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/line-chart/line-chart-multi-axis-and-split.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.line-chart]
--- */
import {test} from '@playwright/test';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {filterPasses} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areasSameSize, hasArea, hasNoArea, hoverArea, legendLists, noErrors, painted, pickFromAreaContextMenu, readingBetween, readingHigher, readingIs, readingReads, repainted, reportsNoError, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Line chart multi-axis layout and splitting into series", () => {
  const session = feature(test, "features/viewers/line-chart/line-chart-multi-axis-and-split.feature", import.meta.url);
  test("Line chart multi-axis layout and splitting into series", {tag: ["@journey", "@viewers", "@realizes:viewers.line-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(19, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"]]));
    await session.step(22, "Then 100 rows should pass the filter", () => filterPasses(page, 100));
    await session.step(23, "And the \"charts\" reading of line chart viewer should be 1", () => readingIs(page, "charts", el("line chart viewer"), 1));
    await session.step(24, "And the \"lines\" reading of line chart viewer should be 1", () => readingIs(page, "lines", el("line chart viewer"), 1));
    await session.step(25, "And the \"categories\" reading of line chart viewer should be 1", () => readingIs(page, "categories", el("line chart viewer"), 1));
    await session.step(26, "And the \"split columns\" reading of line chart viewer should be 0", () => readingIs(page, "split columns", el("line chart viewer"), 0));
    await session.step(27, "And the \"multi axis\" reading of line chart viewer should be \"false\"", () => readingReads(page, "multi axis", el("line chart viewer"), "false"));
    await session.step(28, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await run.scenario("Two Y columns get a chart box each, and Multi Axis folds them into one", async () => {
      await session.step(31, "When user sets \"yColumnNames\" property of line chart viewer to \"Chemical Space X, TPSA\"", () => setProperty(page, "yColumnNames", el("line chart viewer"), "Chemical Space X, TPSA"));
      await session.step(32, "Then the \"charts\" reading of line chart viewer should be 2", () => readingIs(page, "charts", el("line chart viewer"), 2));
      await session.step(33, "And the \"lines\" reading of line chart viewer should be 2", () => readingIs(page, "lines", el("line chart viewer"), 2));
      await session.step(34, "And line chart viewer should have a \"chart 2\" area", () => hasArea(page, el("line chart viewer"), "chart 2"));
      await session.step(35, "And line chart viewer should have a 'chart \"Chemical Space X\"' area", () => hasArea(page, el("line chart viewer"), "chart \"Chemical Space X\""));
      await session.step(36, "And line chart viewer should have a 'chart \"TPSA\"' area", () => hasArea(page, el("line chart viewer"), "chart \"TPSA\""));
      await session.step(37, "And the \"y axes\" reading of line chart viewer should be 2", () => readingIs(page, "y axes", el("line chart viewer"), 2));
      await session.step(38, "When user sets \"multiAxis\" property of line chart viewer to \"true\"", () => setProperty(page, "multiAxis", el("line chart viewer"), "true"));
      await session.step(39, "Then the \"charts\" reading of line chart viewer should be 1", () => readingIs(page, "charts", el("line chart viewer"), 1));
      await session.step(40, "And line chart viewer should not have a \"chart 2\" area", () => hasNoArea(page, el("line chart viewer"), "chart 2"));
      await session.step(41, "And the 'chart \"Chemical Space X\"' and 'chart \"TPSA\"' areas of line chart viewer should be the same height", () => areasSameSize(page, "chart \"Chemical Space X\"", "chart \"TPSA\"", el("line chart viewer"), "height"));
      await session.step(42, "And the \"y axes\" reading of line chart viewer should be 2", () => readingIs(page, "y axes", el("line chart viewer"), 2));
      await session.step(43, "And line chart viewer should have a \"y2 axis\" area", () => hasArea(page, el("line chart viewer"), "y2 axis"));
      await session.step(44, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(45, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["multiAxis","false"],["yColumnNames","Chemical Space X"]]));
      await session.step(48, "Then the \"charts\" reading of line chart viewer should be 1", () => readingIs(page, "charts", el("line chart viewer"), 1));
      await session.step(49, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Y Global Scale replaces the pair of scales with one that covers both columns", async () => {
      await session.step(52, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["yColumnNames","Chemical Space X, TPSA"],["multiAxis","true"]]));
      await session.step(55, "Then the \"y axes\" reading of line chart viewer should be 2", () => readingIs(page, "y axes", el("line chart viewer"), 2));
      await session.step(56, "And line chart viewer should have a \"y2 axis\" area", () => hasArea(page, el("line chart viewer"), "y2 axis"));
      await session.step(57, "And the 'y axis max of \"Chemical Space X\"' reading of line chart viewer should be between 15 and 16", () => readingBetween(page, "y axis max of \"Chemical Space X\"", el("line chart viewer"), 15, 16));
      await session.step(58, "When user sets \"yGlobalScale\" property of line chart viewer to \"true\"", () => setProperty(page, "yGlobalScale", el("line chart viewer"), "true"));
      await session.step(59, "Then the \"y axes\" reading of line chart viewer should be 1", () => readingIs(page, "y axes", el("line chart viewer"), 1));
      await session.step(60, "And line chart viewer should not have a \"y2 axis\" area", () => hasNoArea(page, el("line chart viewer"), "y2 axis"));
      await session.step(61, "And the 'y axis max of \"Chemical Space X\"' reading of line chart viewer should be higher than before", () => readingHigher(page, "y axis max of \"Chemical Space X\"", el("line chart viewer")));
      await session.step(62, "And the 'y axis max of \"Chemical Space X\"' reading of line chart viewer should be between 130 and 133", () => readingBetween(page, "y axis max of \"Chemical Space X\"", el("line chart viewer"), 130, 133));
      await session.step(63, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(64, "When user sets \"yGlobalScale\" property of line chart viewer to \"false\"", () => setProperty(page, "yGlobalScale", el("line chart viewer"), "false"));
      await session.step(65, "Then the \"y axes\" reading of line chart viewer should be 2", () => readingIs(page, "y axes", el("line chart viewer"), 2));
      await session.step(66, "And the 'y axis max of \"Chemical Space X\"' reading of line chart viewer should be between 15 and 16", () => readingBetween(page, "y axis max of \"Chemical Space X\"", el("line chart viewer"), 15, 16));
      await session.step(67, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["multiAxis","false"],["yColumnNames","Chemical Space X"]]));
      await session.step(70, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("One split column draws one series per category", async () => {
      await session.step(73, "When user sets \"splitColumnNames\" property of line chart viewer to \"Stereo Category\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), "Stereo Category"));
      await session.step(74, "Then the \"split columns\" reading of line chart viewer should be 1", () => readingIs(page, "split columns", el("line chart viewer"), 1));
      await session.step(75, "And the \"lines\" reading of line chart viewer should be 5", () => readingIs(page, "lines", el("line chart viewer"), 5));
      await session.step(76, "And the \"categories\" reading of line chart viewer should be 5", () => readingIs(page, "categories", el("line chart viewer"), 5));
      await session.step(77, "And the legend of line chart viewer should list 5 items", () => legendLists(page, el("line chart viewer"), 5));
      await session.step(78, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(79, "And line chart viewer should be painted", () => painted(page, el("line chart viewer")));
      await session.step(80, "When user sets \"splitColumnNames\" property of line chart viewer to \"\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), ""));
      await session.step(81, "Then the \"lines\" reading of line chart viewer should be 1", () => readingIs(page, "lines", el("line chart viewer"), 1));
      await session.step(82, "And the \"categories\" reading of line chart viewer should be 1", () => readingIs(page, "categories", el("line chart viewer"), 1));
      await session.step(83, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A second split column draws the combinations present, not the columns multiplied", async () => {
      await session.step(86, "When user sets \"splitColumnNames\" property of line chart viewer to \"Stereo Category\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), "Stereo Category"));
      await session.step(87, "Then the \"lines\" reading of line chart viewer should be 5", () => readingIs(page, "lines", el("line chart viewer"), 5));
      await session.step(88, "When user sets \"splitColumnNames\" property of line chart viewer to \"Stereo Category, Series\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), "Stereo Category, Series"));
      await session.step(89, "Then the \"split columns\" reading of line chart viewer should be 2", () => readingIs(page, "split columns", el("line chart viewer"), 2));
      await session.step(90, "And the \"lines\" reading of line chart viewer should be 12", () => readingIs(page, "lines", el("line chart viewer"), 12));
      await session.step(91, "And the \"categories\" reading of line chart viewer should be 12", () => readingIs(page, "categories", el("line chart viewer"), 12));
      await session.step(92, "And line chart viewer should be painted", () => painted(page, el("line chart viewer")));
      await session.step(93, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
      await session.step(94, "When user hovers over the \"plot\" area of line chart viewer", () => hoverArea(page, "plot", el("line chart viewer")));
      await session.step(95, "Then no errors should have been logged", () => noErrors(page));
      await session.step(96, "When user sets \"splitColumnNames\" property of line chart viewer to \"\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), ""));
      await session.step(97, "Then the \"lines\" reading of line chart viewer should be 1", () => readingIs(page, "lines", el("line chart viewer"), 1));
      await session.step(98, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Every further split column adds only the combinations the rows carry", async () => {
      await session.step(101, "When user sets \"splitColumnNames\" property of line chart viewer to \"Stereo Category, R1\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), "Stereo Category, R1"));
      await session.step(102, "Then the \"lines\" reading of line chart viewer should be 91", () => readingIs(page, "lines", el("line chart viewer"), 91));
      await session.step(103, "And the \"split columns\" reading of line chart viewer should be 2", () => readingIs(page, "split columns", el("line chart viewer"), 2));
      await session.step(104, "When user sets \"splitColumnNames\" property of line chart viewer to \"Stereo Category, R1, R2\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), "Stereo Category, R1, R2"));
      await session.step(105, "Then the \"lines\" reading of line chart viewer should be 94", () => readingIs(page, "lines", el("line chart viewer"), 94));
      await session.step(106, "And the \"split columns\" reading of line chart viewer should be 3", () => readingIs(page, "split columns", el("line chart viewer"), 3));
      await session.step(107, "When user sets \"splitColumnNames\" property of line chart viewer to \"Stereo Category, R1, R2, R3\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), "Stereo Category, R1, R2, R3"));
      await session.step(108, "Then the \"lines\" reading of line chart viewer should be 96", () => readingIs(page, "lines", el("line chart viewer"), 96));
      await session.step(109, "And the \"split columns\" reading of line chart viewer should be 4", () => readingIs(page, "split columns", el("line chart viewer"), 4));
      await session.step(110, "And the \"rows shown\" reading of line chart viewer should be 100", () => readingIs(page, "rows shown", el("line chart viewer"), 100));
      await session.step(111, "And line chart viewer should be painted", () => painted(page, el("line chart viewer")));
      await session.step(112, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
      await session.step(113, "When user sets \"splitColumnNames\" property of line chart viewer to \"\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), ""));
      await session.step(114, "Then the \"lines\" reading of line chart viewer should be 1", () => readingIs(page, "lines", el("line chart viewer"), 1));
      await session.step(115, "And the \"split columns\" reading of line chart viewer should be 0", () => readingIs(page, "split columns", el("line chart viewer"), 0));
      await session.step(116, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A split multiplies every Y column's series", async () => {
      await session.step(119, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["yColumnNames","Chemical Space X, TPSA"],["splitColumnNames","Stereo Category"]]));
      await session.step(122, "Then the \"charts\" reading of line chart viewer should be 2", () => readingIs(page, "charts", el("line chart viewer"), 2));
      await session.step(123, "And the \"lines\" reading of line chart viewer should be 10", () => readingIs(page, "lines", el("line chart viewer"), 10));
      await session.step(124, "And the \"categories\" reading of line chart viewer should be 5", () => readingIs(page, "categories", el("line chart viewer"), 5));
      await session.step(125, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["splitColumnNames",""],["yColumnNames","Chemical Space X"]]));
      await session.step(128, "Then the \"lines\" reading of line chart viewer should be 1", () => readingIs(page, "lines", el("line chart viewer"), 1));
      await session.step(129, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Hide other charts on the second chart leaves that column alone", async () => {
      await session.step(132, "When user sets \"yColumnNames\" property of line chart viewer to \"Chemical Space X, Chemical Space Y, TPSA\"", () => setProperty(page, "yColumnNames", el("line chart viewer"), "Chemical Space X, Chemical Space Y, TPSA"));
      await session.step(133, "Then the \"charts\" reading of line chart viewer should be 3", () => readingIs(page, "charts", el("line chart viewer"), 3));
      await session.step(134, "And line chart viewer should have a \"chart 3\" area", () => hasArea(page, el("line chart viewer"), "chart 3"));
      await session.step(135, "When user picks \"Chemical Space Y > Hide other charts\" from the context menu of the \"chart 2\" area of line chart viewer", () => pickFromAreaContextMenu(page, "Chemical Space Y > Hide other charts", "chart 2", el("line chart viewer")));
      await session.step(136, "Then the \"charts\" reading of line chart viewer should be 1", () => readingIs(page, "charts", el("line chart viewer"), 1));
      await session.step(137, "And the \"y columns\" reading of line chart viewer should be \"Chemical Space Y\"", () => readingReads(page, "y columns", el("line chart viewer"), "Chemical Space Y"));
      await session.step(138, "And line chart viewer should not have a \"chart 2\" area", () => hasNoArea(page, el("line chart viewer"), "chart 2"));
      await session.step(139, "And line chart viewer should be painted", () => painted(page, el("line chart viewer")));
      await session.step(140, "When user sets \"yColumnNames\" property of line chart viewer to \"Chemical Space X\"", () => setProperty(page, "yColumnNames", el("line chart viewer"), "Chemical Space X"));
      await session.step(141, "Then the \"y columns\" reading of line chart viewer should be \"Chemical Space X\"", () => readingReads(page, "y columns", el("line chart viewer"), "Chemical Space X"));
      await session.step(142, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
