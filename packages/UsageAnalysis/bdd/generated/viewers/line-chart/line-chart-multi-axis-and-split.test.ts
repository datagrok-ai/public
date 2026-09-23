/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/line-chart/line-chart-multi-axis-and-split.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.line-chart]
--- */
import {test} from '@playwright/test';
import '../../../bindings/grid.js';
import '../../../bindings/nx.js';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, isExpanded, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPasses} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areasSameSize, hasArea, hasNoArea, hoverArea, legendLists, noErrors, painted, pickFromAreaContextMenu, readingBetween, readingHigher, readingIs, readingReads, repainted, reportsNoError, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {checkedInColumnList, columnListStartsWith, liesWithin, toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Line chart multi-axis layout and splitting into series", () => {
  const session = feature(test, "features/viewers/line-chart/line-chart-multi-axis-and-split.feature", import.meta.url);
  test("Line chart multi-axis layout and splitting into series", {tag: ["@journey", "@viewers", "@realizes:viewers.line-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(22, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"]]), [["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"]]);
    await session.step(25, "Then 100 rows should pass the filter", () => filterPasses(page, 100));
    await session.step(26, "And the \"charts\" reading of line chart viewer should be 1", () => readingIs(page, "charts", el("line chart viewer"), 1));
    await session.step(27, "And the \"lines\" reading of line chart viewer should be 1", () => readingIs(page, "lines", el("line chart viewer"), 1));
    await session.step(28, "And the \"categories\" reading of line chart viewer should be 1", () => readingIs(page, "categories", el("line chart viewer"), 1));
    await session.step(29, "And the \"split columns\" reading of line chart viewer should be 0", () => readingIs(page, "split columns", el("line chart viewer"), 0));
    await session.step(30, "And the \"multi axis\" reading of line chart viewer should be \"false\"", () => readingReads(page, "multi axis", el("line chart viewer"), "false"));
    await session.step(31, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await run.scenario("Two Y columns get a chart box each, and Multi Axis folds them into one", async () => {
      await session.step(34, "When user sets \"yColumnNames\" property of line chart viewer to \"Chemical Space X, TPSA\"", () => setProperty(page, "yColumnNames", el("line chart viewer"), "Chemical Space X, TPSA"));
      await session.step(35, "Then the \"charts\" reading of line chart viewer should be 2", () => readingIs(page, "charts", el("line chart viewer"), 2));
      await session.step(36, "And the \"lines\" reading of line chart viewer should be 2", () => readingIs(page, "lines", el("line chart viewer"), 2));
      await session.step(37, "And line chart viewer should have a \"chart 2\" area", () => hasArea(page, el("line chart viewer"), "chart 2"));
      await session.step(38, "And line chart viewer should have a 'chart \"Chemical Space X\"' area", () => hasArea(page, el("line chart viewer"), "chart \"Chemical Space X\""));
      await session.step(39, "And line chart viewer should have a 'chart \"TPSA\"' area", () => hasArea(page, el("line chart viewer"), "chart \"TPSA\""));
      await session.step(40, "And the \"y axes\" reading of line chart viewer should be 2", () => readingIs(page, "y axes", el("line chart viewer"), 2));
      await session.step(41, "When user sets \"multiAxis\" property of line chart viewer to \"true\"", () => setProperty(page, "multiAxis", el("line chart viewer"), "true"));
      await session.step(42, "Then the \"charts\" reading of line chart viewer should be 1", () => readingIs(page, "charts", el("line chart viewer"), 1));
      await session.step(43, "And line chart viewer should not have a \"chart 2\" area", () => hasNoArea(page, el("line chart viewer"), "chart 2"));
      await session.step(44, "And the 'chart \"Chemical Space X\"' and 'chart \"TPSA\"' areas of line chart viewer should be the same height", () => areasSameSize(page, "chart \"Chemical Space X\"", "chart \"TPSA\"", el("line chart viewer"), "height"));
      await session.step(45, "And the \"y axes\" reading of line chart viewer should be 2", () => readingIs(page, "y axes", el("line chart viewer"), 2));
      await session.step(46, "And line chart viewer should have a \"y2 axis\" area", () => hasArea(page, el("line chart viewer"), "y2 axis"));
      await session.step(47, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(48, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["multiAxis","false"],["yColumnNames","Chemical Space X"]]), [["multiAxis","false"],["yColumnNames","Chemical Space X"]]);
      await session.step(51, "Then the \"charts\" reading of line chart viewer should be 1", () => readingIs(page, "charts", el("line chart viewer"), 1));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Y Global Scale replaces the pair of scales with one that covers both columns", async () => {
      await session.step(55, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["yColumnNames","Chemical Space X, TPSA"],["multiAxis","true"]]), [["yColumnNames","Chemical Space X, TPSA"],["multiAxis","true"]]);
      await session.step(58, "Then the \"y axes\" reading of line chart viewer should be 2", () => readingIs(page, "y axes", el("line chart viewer"), 2));
      await session.step(59, "And line chart viewer should have a \"y2 axis\" area", () => hasArea(page, el("line chart viewer"), "y2 axis"));
      await session.step(60, "And the 'y axis max of \"Chemical Space X\"' reading of line chart viewer should be between 15 and 16", () => readingBetween(page, "y axis max of \"Chemical Space X\"", el("line chart viewer"), 15, 16));
      await session.step(61, "When user sets \"yGlobalScale\" property of line chart viewer to \"true\"", () => setProperty(page, "yGlobalScale", el("line chart viewer"), "true"));
      await session.step(62, "Then the \"y axes\" reading of line chart viewer should be 1", () => readingIs(page, "y axes", el("line chart viewer"), 1));
      await session.step(63, "And line chart viewer should not have a \"y2 axis\" area", () => hasNoArea(page, el("line chart viewer"), "y2 axis"));
      await session.step(64, "And the 'y axis max of \"Chemical Space X\"' reading of line chart viewer should be higher than before", () => readingHigher(page, "y axis max of \"Chemical Space X\"", el("line chart viewer")));
      await session.step(65, "And the 'y axis max of \"Chemical Space X\"' reading of line chart viewer should be between 130 and 133", () => readingBetween(page, "y axis max of \"Chemical Space X\"", el("line chart viewer"), 130, 133));
      await session.step(66, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(67, "When user sets \"yGlobalScale\" property of line chart viewer to \"false\"", () => setProperty(page, "yGlobalScale", el("line chart viewer"), "false"));
      await session.step(68, "Then the \"y axes\" reading of line chart viewer should be 2", () => readingIs(page, "y axes", el("line chart viewer"), 2));
      await session.step(69, "And the 'y axis max of \"Chemical Space X\"' reading of line chart viewer should be between 15 and 16", () => readingBetween(page, "y axis max of \"Chemical Space X\"", el("line chart viewer"), 15, 16));
      await session.step(70, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["multiAxis","false"],["yColumnNames","Chemical Space X"]]), [["multiAxis","false"],["yColumnNames","Chemical Space X"]]);
      await session.step(73, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("One split column draws one series per category", async () => {
      await session.step(76, "When user sets \"splitColumnNames\" property of line chart viewer to \"Stereo Category\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), "Stereo Category"));
      await session.step(77, "Then the \"split columns\" reading of line chart viewer should be 1", () => readingIs(page, "split columns", el("line chart viewer"), 1));
      await session.step(78, "And the \"lines\" reading of line chart viewer should be 5", () => readingIs(page, "lines", el("line chart viewer"), 5));
      await session.step(79, "And the \"categories\" reading of line chart viewer should be 5", () => readingIs(page, "categories", el("line chart viewer"), 5));
      await session.step(80, "And the legend of line chart viewer should list 5 items", () => legendLists(page, el("line chart viewer"), 5));
      await session.step(81, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(82, "And line chart viewer should be painted", () => painted(page, el("line chart viewer")));
      await session.step(83, "When user sets \"splitColumnNames\" property of line chart viewer to \"\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), ""));
      await session.step(84, "Then the \"lines\" reading of line chart viewer should be 1", () => readingIs(page, "lines", el("line chart viewer"), 1));
      await session.step(85, "And the \"categories\" reading of line chart viewer should be 1", () => readingIs(page, "categories", el("line chart viewer"), 1));
      await session.step(86, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A second split column draws the combinations present, not the columns multiplied", async () => {
      await session.step(89, "When user sets \"splitColumnNames\" property of line chart viewer to \"Stereo Category\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), "Stereo Category"));
      await session.step(90, "Then the \"lines\" reading of line chart viewer should be 5", () => readingIs(page, "lines", el("line chart viewer"), 5));
      await session.step(91, "When user sets \"splitColumnNames\" property of line chart viewer to \"Stereo Category, Series\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), "Stereo Category, Series"));
      await session.step(92, "Then the \"split columns\" reading of line chart viewer should be 2", () => readingIs(page, "split columns", el("line chart viewer"), 2));
      await session.step(93, "And the \"lines\" reading of line chart viewer should be 12", () => readingIs(page, "lines", el("line chart viewer"), 12));
      await session.step(94, "And the \"categories\" reading of line chart viewer should be 12", () => readingIs(page, "categories", el("line chart viewer"), 12));
      await session.step(95, "And line chart viewer should be painted", () => painted(page, el("line chart viewer")));
      await session.step(96, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
      await session.step(97, "When user hovers over the \"plot\" area of line chart viewer", () => hoverArea(page, "plot", el("line chart viewer")));
      await session.step(98, "Then no errors should have been logged", () => noErrors(page));
      await session.step(99, "When user sets \"splitColumnNames\" property of line chart viewer to \"\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), ""));
      await session.step(100, "Then the \"lines\" reading of line chart viewer should be 1", () => readingIs(page, "lines", el("line chart viewer"), 1));
      await session.step(101, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Multi Axis with two split columns keeps drawing, and a hover over it logs nothing", async () => {
      await session.step(104, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["yColumnNames","Chemical Space X, TPSA"],["multiAxis","true"]]), [["yColumnNames","Chemical Space X, TPSA"],["multiAxis","true"]]);
      await session.step(107, "Then the \"charts\" reading of line chart viewer should be 1", () => readingIs(page, "charts", el("line chart viewer"), 1));
      await session.step(108, "And the \"y axes\" reading of line chart viewer should be 2", () => readingIs(page, "y axes", el("line chart viewer"), 2));
      await session.step(109, "And no errors should have been logged", () => noErrors(page));
      await session.step(110, "When user sets \"splitColumnNames\" property of line chart viewer to \"Stereo Category\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), "Stereo Category"));
      await session.step(111, "Then the \"split columns\" reading of line chart viewer should be 1", () => readingIs(page, "split columns", el("line chart viewer"), 1));
      await session.step(112, "And the \"categories\" reading of line chart viewer should be 5", () => readingIs(page, "categories", el("line chart viewer"), 5));
      await session.step(113, "And line chart viewer should be painted", () => painted(page, el("line chart viewer")));
      await session.step(114, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
      await session.step(115, "When user sets \"splitColumnNames\" property of line chart viewer to \"Stereo Category, Series\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), "Stereo Category, Series"));
      await session.step(116, "Then the \"split columns\" reading of line chart viewer should be 2", () => readingIs(page, "split columns", el("line chart viewer"), 2));
      await session.step(117, "And the \"categories\" reading of line chart viewer should be 12", () => readingIs(page, "categories", el("line chart viewer"), 12));
      await session.step(118, "And the \"lines\" reading of line chart viewer should be 24", () => readingIs(page, "lines", el("line chart viewer"), 24));
      await session.step(119, "And the \"charts\" reading of line chart viewer should be 1", () => readingIs(page, "charts", el("line chart viewer"), 1));
      await session.step(120, "And line chart viewer should be painted", () => painted(page, el("line chart viewer")));
      await session.step(121, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
      await session.step(122, "When user hovers over the \"plot\" area of line chart viewer", () => hoverArea(page, "plot", el("line chart viewer")));
      await session.step(123, "Then no errors should have been logged", () => noErrors(page));
      await session.step(124, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["splitColumnNames",""],["multiAxis","false"],["yColumnNames","Chemical Space X"]]), [["splitColumnNames",""],["multiAxis","false"],["yColumnNames","Chemical Space X"]]);
      await session.step(128, "Then the \"charts\" reading of line chart viewer should be 1", () => readingIs(page, "charts", el("line chart viewer"), 1));
      await session.step(129, "And the \"lines\" reading of line chart viewer should be 1", () => readingIs(page, "lines", el("line chart viewer"), 1));
      await session.step(130, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An edit in the Y column list keeps the three Y columns, and its search box sits inside the list", async () => {
      await session.step(133, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["yColumnNames","Chemical Space X, Chemical Space Y, TPSA"],["multiAxis","true"]]), [["yColumnNames","Chemical Space X, Chemical Space Y, TPSA"],["multiAxis","true"]]);
      await session.step(136, "Then the \"y columns\" reading of line chart viewer should be \"Chemical Space X, Chemical Space Y, TPSA\"", () => readingReads(page, "y columns", el("line chart viewer"), "Chemical Space X, Chemical Space Y, TPSA"));
      await session.step(137, "When user clicks on settings icon of line chart viewer", () => clickOn(page, el("settings icon of line chart viewer")));
      await session.step(138, "Given \"Y Axis\" category in context panel is expanded", () => isExpanded(page, el("\"Y Axis\" category in context panel")));
      await session.step(139, "When user clicks on \"...\" button in \"Y\" property", () => clickOn(page, el("\"...\" button in \"Y\" property")));
      await session.step(140, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
      await session.step(141, "And the \"Chemical Space Y\" column should be checked in the column list of \"Select columns...\" dialog", () => checkedInColumnList(page, "Chemical Space Y", el("\"Select columns...\" dialog")));
      await session.step(142, "And \"Search\" input in \"Select columns...\" dialog should lie within \"Select columns...\" dialog", () => liesWithin(page, el("\"Search\" input in \"Select columns...\" dialog"), el("\"Select columns...\" dialog")));
      await session.step(143, "When user types \"TPSA\" into \"Search\" input in \"Select columns...\" dialog", () => typeInto(page, "TPSA", el("\"Search\" input in \"Select columns...\" dialog")));
      await session.step(144, "Then the column list of \"Select columns...\" dialog should start with \"TPSA\"", () => columnListStartsWith(page, el("\"Select columns...\" dialog"), "TPSA"));
      await session.step(145, "And the \"TPSA\" column should be checked in the column list of \"Select columns...\" dialog", () => checkedInColumnList(page, "TPSA", el("\"Select columns...\" dialog")));
      await session.step(146, "When user toggles the \"TPSA\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "TPSA", el("\"Select columns...\" dialog")));
      await session.step(147, "And user toggles the \"TPSA\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "TPSA", el("\"Select columns...\" dialog")));
      await session.step(148, "Then the \"TPSA\" column should be checked in the column list of \"Select columns...\" dialog", () => checkedInColumnList(page, "TPSA", el("\"Select columns...\" dialog")));
      await session.step(149, "And \"Search\" input in \"Select columns...\" dialog should lie within \"Select columns...\" dialog", () => liesWithin(page, el("\"Search\" input in \"Select columns...\" dialog"), el("\"Select columns...\" dialog")));
      await session.step(150, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(151, "Then \"Select columns...\" dialog should be absent", () => shouldBe(page, el("\"Select columns...\" dialog"), "absent"));
      await session.step(152, "And the \"y columns\" reading of line chart viewer should be \"Chemical Space X, Chemical Space Y, TPSA\"", () => readingReads(page, "y columns", el("line chart viewer"), "Chemical Space X, Chemical Space Y, TPSA"));
      await session.step(153, "And the \"charts\" reading of line chart viewer should be 1", () => readingIs(page, "charts", el("line chart viewer"), 1));
      await session.step(154, "And the \"multi axis\" reading of line chart viewer should be \"true\"", () => readingReads(page, "multi axis", el("line chart viewer"), "true"));
      await session.step(155, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["multiAxis","false"],["yColumnNames","Chemical Space X"]]), [["multiAxis","false"],["yColumnNames","Chemical Space X"]]);
      await session.step(158, "Then the \"y columns\" reading of line chart viewer should be \"Chemical Space X\"", () => readingReads(page, "y columns", el("line chart viewer"), "Chemical Space X"));
      await session.step(159, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Every further split column adds only the combinations the rows carry", async () => {
      await session.step(162, "When user sets \"splitColumnNames\" property of line chart viewer to \"Stereo Category, R1\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), "Stereo Category, R1"));
      await session.step(163, "Then the \"lines\" reading of line chart viewer should be 91", () => readingIs(page, "lines", el("line chart viewer"), 91));
      await session.step(164, "And the \"split columns\" reading of line chart viewer should be 2", () => readingIs(page, "split columns", el("line chart viewer"), 2));
      await session.step(165, "When user sets \"splitColumnNames\" property of line chart viewer to \"Stereo Category, R1, R2\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), "Stereo Category, R1, R2"));
      await session.step(166, "Then the \"lines\" reading of line chart viewer should be 94", () => readingIs(page, "lines", el("line chart viewer"), 94));
      await session.step(167, "And the \"split columns\" reading of line chart viewer should be 3", () => readingIs(page, "split columns", el("line chart viewer"), 3));
      await session.step(168, "When user sets \"splitColumnNames\" property of line chart viewer to \"Stereo Category, R1, R2, R3\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), "Stereo Category, R1, R2, R3"));
      await session.step(169, "Then the \"lines\" reading of line chart viewer should be 96", () => readingIs(page, "lines", el("line chart viewer"), 96));
      await session.step(170, "And the \"split columns\" reading of line chart viewer should be 4", () => readingIs(page, "split columns", el("line chart viewer"), 4));
      await session.step(171, "And the \"rows shown\" reading of line chart viewer should be 100", () => readingIs(page, "rows shown", el("line chart viewer"), 100));
      await session.step(172, "And line chart viewer should be painted", () => painted(page, el("line chart viewer")));
      await session.step(173, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
      await session.step(174, "When user sets \"splitColumnNames\" property of line chart viewer to \"\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), ""));
      await session.step(175, "Then the \"lines\" reading of line chart viewer should be 1", () => readingIs(page, "lines", el("line chart viewer"), 1));
      await session.step(176, "And the \"split columns\" reading of line chart viewer should be 0", () => readingIs(page, "split columns", el("line chart viewer"), 0));
      await session.step(177, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A split multiplies every Y column's series", async () => {
      await session.step(180, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["yColumnNames","Chemical Space X, TPSA"],["splitColumnNames","Stereo Category"]]), [["yColumnNames","Chemical Space X, TPSA"],["splitColumnNames","Stereo Category"]]);
      await session.step(183, "Then the \"charts\" reading of line chart viewer should be 2", () => readingIs(page, "charts", el("line chart viewer"), 2));
      await session.step(184, "And the \"lines\" reading of line chart viewer should be 10", () => readingIs(page, "lines", el("line chart viewer"), 10));
      await session.step(185, "And the \"categories\" reading of line chart viewer should be 5", () => readingIs(page, "categories", el("line chart viewer"), 5));
      await session.step(186, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["splitColumnNames",""],["yColumnNames","Chemical Space X"]]), [["splitColumnNames",""],["yColumnNames","Chemical Space X"]]);
      await session.step(189, "Then the \"lines\" reading of line chart viewer should be 1", () => readingIs(page, "lines", el("line chart viewer"), 1));
      await session.step(190, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Hide other charts on the second chart leaves that column alone", async () => {
      await session.step(193, "When user sets \"yColumnNames\" property of line chart viewer to \"Chemical Space X, Chemical Space Y, TPSA\"", () => setProperty(page, "yColumnNames", el("line chart viewer"), "Chemical Space X, Chemical Space Y, TPSA"));
      await session.step(194, "Then the \"charts\" reading of line chart viewer should be 3", () => readingIs(page, "charts", el("line chart viewer"), 3));
      await session.step(195, "And line chart viewer should have a \"chart 3\" area", () => hasArea(page, el("line chart viewer"), "chart 3"));
      await session.step(196, "When user picks \"Chemical Space Y > Hide other charts\" from the context menu of the \"chart 2\" area of line chart viewer", () => pickFromAreaContextMenu(page, "Chemical Space Y > Hide other charts", "chart 2", el("line chart viewer")));
      await session.step(197, "Then the \"charts\" reading of line chart viewer should be 1", () => readingIs(page, "charts", el("line chart viewer"), 1));
      await session.step(198, "And the \"y columns\" reading of line chart viewer should be \"Chemical Space Y\"", () => readingReads(page, "y columns", el("line chart viewer"), "Chemical Space Y"));
      await session.step(199, "And line chart viewer should not have a \"chart 2\" area", () => hasNoArea(page, el("line chart viewer"), "chart 2"));
      await session.step(200, "And line chart viewer should be painted", () => painted(page, el("line chart viewer")));
      await session.step(201, "When user sets \"yColumnNames\" property of line chart viewer to \"Chemical Space X\"", () => setProperty(page, "yColumnNames", el("line chart viewer"), "Chemical Space X"));
      await session.step(202, "Then the \"y columns\" reading of line chart viewer should be \"Chemical Space X\"", () => readingReads(page, "y columns", el("line chart viewer"), "Chemical Space X"));
      await session.step(203, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
