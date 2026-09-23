/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/formula-lines/formula-lines-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot, viewers.line-chart]
--- */
import {test} from '@playwright/test';
import '../../../bindings/grid.js';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {renameColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, hasArea, hoverArea, noBalloons, noErrors, oneTooltip, pickFromAreaContextMenu, pointerAway, propertyShouldBe, propertyShouldContain, readingIs, repainted, resizeTo, setProperties, setProperty, tooltipSomeColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Formula lines regression checks", () => {
  const session = feature(test, "features/viewers/formula-lines/formula-lines-regressions.feature", import.meta.url);
  test("A column on both axes is renamed while a line joins them", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.line-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "Given user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE"],["yColumnName","AGE"],["formulaLines","[{\"type\":\"line\",\"formula\":\"${AGE} = ${AGE}\"}]"]]), [["xColumnName","AGE"],["yColumnName","AGE"],["formulaLines","[{\"type\":\"line\",\"formula\":\"${AGE} = ${AGE}\"}]"]]);
    await session.step(22, "Then the \"formula lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "formula lines", el("scatter plot viewer"), 1));
    await session.step(23, "When user renames \"AGE\" column to \"AGE_RENAMED\"", () => renameColumn(page, "AGE", "AGE_RENAMED"));
    await session.step(24, "Then \"xColumnName\" property of scatter plot viewer should be \"AGE_RENAMED\"", () => propertyShouldBe(page, "xColumnName", el("scatter plot viewer"), "AGE_RENAMED"));
    await session.step(25, "And \"yColumnName\" property of scatter plot viewer should be \"AGE_RENAMED\"", () => propertyShouldBe(page, "yColumnName", el("scatter plot viewer"), "AGE_RENAMED"));
    await session.step(26, "And \"formulaLines\" property of scatter plot viewer should contain \"${AGE_RENAMED} = ${AGE_RENAMED}\"", () => propertyShouldContain(page, "formulaLines", el("scatter plot viewer"), "${AGE_RENAMED} = ${AGE_RENAMED}"));
    await session.step(27, "And the \"formula lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "formula lines", el("scatter plot viewer"), 1));
    await session.step(28, "And no errors should have been logged", () => noErrors(page));
    await session.step(29, "When user renames \"AGE_RENAMED\" column to \"AGE\"", () => renameColumn(page, "AGE_RENAMED", "AGE"));
    await session.step(30, "Then \"xColumnName\" property of scatter plot viewer should be \"AGE\"", () => propertyShouldBe(page, "xColumnName", el("scatter plot viewer"), "AGE"));
    await session.step(31, "And \"formulaLines\" property of scatter plot viewer should contain \"${AGE} = ${AGE}\"", () => propertyShouldContain(page, "formulaLines", el("scatter plot viewer"), "${AGE} = ${AGE}"));
    await session.step(32, "And the \"formula lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "formula lines", el("scatter plot viewer"), 1));
    await session.step(33, "And no errors should have been logged", () => noErrors(page));
  });
  test("Renaming a column rewrites the scatter plot line that uses it", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.line-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(36, "Given user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE"],["yColumnName","HEIGHT"],["formulaLines","[{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${AGE}\"}]"]]), [["xColumnName","AGE"],["yColumnName","HEIGHT"],["formulaLines","[{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${AGE}\"}]"]]);
    await session.step(40, "Then the \"formula lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "formula lines", el("scatter plot viewer"), 1));
    await session.step(41, "When user renames \"AGE\" column to \"AGE_R\"", () => renameColumn(page, "AGE", "AGE_R"));
    await session.step(42, "Then \"formulaLines\" property of scatter plot viewer should contain \"${HEIGHT} = ${AGE_R}\"", () => propertyShouldContain(page, "formulaLines", el("scatter plot viewer"), "${HEIGHT} = ${AGE_R}"));
    await session.step(43, "And the \"formula lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "formula lines", el("scatter plot viewer"), 1));
    await session.step(44, "When user renames \"AGE_R\" column to \"AGE\"", () => renameColumn(page, "AGE_R", "AGE"));
    await session.step(45, "Then \"formulaLines\" property of scatter plot viewer should contain \"${HEIGHT} = ${AGE}\"", () => propertyShouldContain(page, "formulaLines", el("scatter plot viewer"), "${HEIGHT} = ${AGE}"));
    await session.step(46, "And no errors should have been logged", () => noErrors(page));
  });
  test("Renaming a column rewrites the line chart line added from its axis", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.line-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(49, "Given user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","AGE"],["yColumnNames","HEIGHT"]]), [["xColumnName","AGE"],["yColumnNames","HEIGHT"]]);
    await session.step(52, "And user resizes line chart viewer to 800 by 500", () => resizeTo(page, el("line chart viewer"), 800, 500));
    await session.step(53, "When user picks \"Annotations > Add Line\" from the context menu of the \"bottom edge of y axis\" area of line chart viewer", () => pickFromAreaContextMenu(page, "Annotations > Add Line", "bottom edge of y axis", el("line chart viewer")));
    await session.step(54, "And user clicks OK button in \"Formula Lines\" dialog", () => clickOn(page, el("OK button in \"Formula Lines\" dialog")));
    await session.step(55, "Then \"formulaLines\" property of line chart viewer should contain \"${avg(HEIGHT)} = 168.8\"", () => propertyShouldContain(page, "formulaLines", el("line chart viewer"), "${avg(HEIGHT)} = 168.8"));
    await session.step(56, "And line chart viewer should have a \"formula line avg(HEIGHT) = 168.8\" area", () => hasArea(page, el("line chart viewer"), "formula line avg(HEIGHT) = 168.8"));
    await session.step(57, "When user renames \"HEIGHT\" column to \"HEIGHT_R\"", () => renameColumn(page, "HEIGHT", "HEIGHT_R"));
    await session.step(58, "Then \"formulaLines\" property of line chart viewer should contain \"${avg(HEIGHT_R)} = 168.8\"", () => propertyShouldContain(page, "formulaLines", el("line chart viewer"), "${avg(HEIGHT_R)} = 168.8"));
    await session.step(59, "And line chart viewer should have a \"formula line avg(HEIGHT_R) = 168.8\" area", () => hasArea(page, el("line chart viewer"), "formula line avg(HEIGHT_R) = 168.8"));
    await session.step(60, "When user renames \"HEIGHT_R\" column to \"HEIGHT\"", () => renameColumn(page, "HEIGHT_R", "HEIGHT"));
    await session.step(61, "Then \"formulaLines\" property of line chart viewer should contain \"${avg(HEIGHT)} = 168.8\"", () => propertyShouldContain(page, "formulaLines", el("line chart viewer"), "${avg(HEIGHT)} = 168.8"));
    await session.step(62, "And no errors should have been logged", () => noErrors(page));
  });
  test("A line survives a change of the axis columns untouched", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.line-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(65, "Given user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"],["formulaLines","[{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${WEIGHT}\"}]"]]), [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"],["formulaLines","[{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${WEIGHT}\"}]"]]);
    await session.step(69, "Then the \"formula lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "formula lines", el("scatter plot viewer"), 1));
    await session.step(70, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["xColumnName","AGE"],["yColumnName","WEIGHT"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"]]);
    await session.step(73, "Then \"formulaLines\" property of scatter plot viewer should be '[{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${WEIGHT}\"}]'", () => propertyShouldBe(page, "formulaLines", el("scatter plot viewer"), "[{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${WEIGHT}\"}]"));
    await session.step(74, "And the \"formula lines\" reading of scatter plot viewer should be 0", () => readingIs(page, "formula lines", el("scatter plot viewer"), 0));
    await session.step(75, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"]]), [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"]]);
    await session.step(78, "Then the \"formula lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "formula lines", el("scatter plot viewer"), 1));
    await session.step(79, "And no errors should have been logged", () => noErrors(page));
  });
  test("A horizontal band survives a logarithmic Y axis on the scatter plot", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.line-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(82, "Given user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"],["formulaLines","[{\"type\":\"band\",\"formula\":\"${HEIGHT} in (160.9, 177.6)\",\"orientation\":\"Horizontal\",\"column2\":\"WEIGHT\"}]"]]), [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"],["formulaLines","[{\"type\":\"band\",\"formula\":\"${HEIGHT} in (160.9, 177.6)\",\"orientation\":\"Horizontal\",\"column2\":\"WEIGHT\"}]"]]);
    await session.step(86, "Then the \"formula lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "formula lines", el("scatter plot viewer"), 1));
    await session.step(87, "When user sets \"yAxisType\" property of scatter plot viewer to \"logarithmic\"", () => setProperty(page, "yAxisType", el("scatter plot viewer"), "logarithmic"));
    await session.step(88, "Then the \"formula lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "formula lines", el("scatter plot viewer"), 1));
    await session.step(89, "And scatter plot viewer should have a \"formula band 1\" area", () => hasArea(page, el("scatter plot viewer"), "formula band 1"));
    await session.step(90, "And scatter plot viewer should have repainted", () => repainted(page, el("scatter plot viewer")));
    await session.step(91, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(92, "When user sets \"yAxisType\" property of scatter plot viewer to \"linear\"", () => setProperty(page, "yAxisType", el("scatter plot viewer"), "linear"));
    await session.step(93, "Then the \"formula lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "formula lines", el("scatter plot viewer"), 1));
    await session.step(94, "And no errors should have been logged", () => noErrors(page));
  });
  test("A line survives a logarithmic Y axis on the line chart", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.line-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(97, "Given user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","AGE"],["yColumnNames","HEIGHT"]]), [["xColumnName","AGE"],["yColumnNames","HEIGHT"]]);
    await session.step(100, "And user resizes line chart viewer to 800 by 500", () => resizeTo(page, el("line chart viewer"), 800, 500));
    await session.step(101, "When user picks \"Annotations > Add Line\" from the context menu of the \"bottom edge of y axis\" area of line chart viewer", () => pickFromAreaContextMenu(page, "Annotations > Add Line", "bottom edge of y axis", el("line chart viewer")));
    await session.step(102, "And user clicks OK button in \"Formula Lines\" dialog", () => clickOn(page, el("OK button in \"Formula Lines\" dialog")));
    await session.step(103, "Then line chart viewer should have a \"formula line avg(HEIGHT) = 168.8\" area", () => hasArea(page, el("line chart viewer"), "formula line avg(HEIGHT) = 168.8"));
    await session.step(104, "When user sets \"yAxisType\" property of line chart viewer to \"logarithmic\"", () => setProperty(page, "yAxisType", el("line chart viewer"), "logarithmic"));
    await session.step(105, "Then line chart viewer should have a \"formula line avg(HEIGHT) = 168.8\" area", () => hasArea(page, el("line chart viewer"), "formula line avg(HEIGHT) = 168.8"));
    await session.step(106, "And the \"formula lines\" reading of line chart viewer should be 1", () => readingIs(page, "formula lines", el("line chart viewer"), 1));
    await session.step(107, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(108, "When user sets \"yAxisType\" property of line chart viewer to \"linear\"", () => setProperty(page, "yAxisType", el("line chart viewer"), "linear"));
    await session.step(109, "Then line chart viewer should have a \"formula line avg(HEIGHT) = 168.8\" area", () => hasArea(page, el("line chart viewer"), "formula line avg(HEIGHT) = 168.8"));
    await session.step(110, "And no errors should have been logged", () => noErrors(page));
  });
  test("Hovering markers next to a formula line raises no errors", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.line-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(113, "Given user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"],["formulaLines","[{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${WEIGHT}\"}]"]]), [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"],["formulaLines","[{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${WEIGHT}\"}]"]]);
    await session.step(117, "And user resizes scatter plot viewer to 800 by 500", () => resizeTo(page, el("scatter plot viewer"), 800, 500));
    await session.step(118, "Then the \"formula lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "formula lines", el("scatter plot viewer"), 1));
    await session.step(119, "When user hovers over the \"marker of row 1\" area of scatter plot viewer", () => hoverArea(page, "marker of row 1", el("scatter plot viewer")));
    await session.step(120, "Then the tooltip should show some columns", () => tooltipSomeColumns(page));
    await session.step(121, "When user hovers over the \"marker of row 2\" area of scatter plot viewer", () => hoverArea(page, "marker of row 2", el("scatter plot viewer")));
    await session.step(122, "And user hovers over the \"marker of row 3\" area of scatter plot viewer", () => hoverArea(page, "marker of row 3", el("scatter plot viewer")));
    await session.step(123, "And user hovers over the \"marker of row 4\" area of scatter plot viewer", () => hoverArea(page, "marker of row 4", el("scatter plot viewer")));
    await session.step(124, "And user hovers over the \"marker of row 5\" area of scatter plot viewer", () => hoverArea(page, "marker of row 5", el("scatter plot viewer")));
    await session.step(125, "And user hovers over the \"marker of row 6\" area of scatter plot viewer", () => hoverArea(page, "marker of row 6", el("scatter plot viewer")));
    await session.step(126, "And user hovers over the \"marker of row 7\" area of scatter plot viewer", () => hoverArea(page, "marker of row 7", el("scatter plot viewer")));
    await session.step(127, "And user hovers over the \"marker of row 8\" area of scatter plot viewer", () => hoverArea(page, "marker of row 8", el("scatter plot viewer")));
    await session.step(128, "Then the tooltip should show some columns", () => tooltipSomeColumns(page));
    await session.step(129, "And exactly one tooltip should be shown", () => oneTooltip(page));
    await session.step(130, "When user moves the pointer away from scatter plot viewer", () => pointerAway(page, el("scatter plot viewer")));
    await session.step(131, "Then no errors should have been logged", () => noErrors(page));
  });
});
