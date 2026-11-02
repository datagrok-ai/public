/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/formula-lines/formula-lines-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot, viewers.line-chart, viewers.trellis-plot, powerpack.dialogs.formula-lines, public#3177, public#542, public#699, public#3390, public#3389, public#425]
--- */
import {test} from '@playwright/test';
import '../../bindings/add-new-column.js';
import '../../bindings/enrichment.js';
import '../../bindings/home.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {renameColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {drawsFormulaLines, lineDrawnAtValue} from '@datagrok-libraries/bdd/bindings/tiers/viewers/formula-lines';
import {addViewerWith, clickArea, hasArea, hoverArea, loadLayout, noBalloons, noErrors, oneTooltip, pickFromAreaContextMenu, pickFromContextMenu, pickFromOpenMenu, pointerAway, propertyShouldBe, propertyShouldContain, readingBetween, readingIs, readingReads, repainted, resizeTo, saveLayout, setProperties, setProperty, tooltipSomeColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {setInnerProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Formula lines regression checks", () => {
  const session = feature(test, "features/formula-lines/formula-lines-regressions.feature", import.meta.url);
  test("A column on both axes is renamed while a line joins them", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.line-chart", "@realizes:viewers.trellis-plot", "@realizes:powerpack.dialogs.formula-lines"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(24, "Given user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE"],["yColumnName","AGE"],["formulaLines","[{\"type\":\"line\",\"formula\":\"${AGE} = ${AGE}\"}]"]]), [["xColumnName","AGE"],["yColumnName","AGE"],["formulaLines","[{\"type\":\"line\",\"formula\":\"${AGE} = ${AGE}\"}]"]]);
    await session.step(28, "Then the \"formula lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "formula lines", el("scatter plot viewer"), 1));
    await session.step(29, "When user renames \"AGE\" column to \"AGE_RENAMED\"", () => renameColumn(page, "AGE", "AGE_RENAMED"));
    await session.step(30, "Then \"xColumnName\" property of scatter plot viewer should be \"AGE_RENAMED\"", () => propertyShouldBe(page, "xColumnName", el("scatter plot viewer"), "AGE_RENAMED"));
    await session.step(31, "And \"yColumnName\" property of scatter plot viewer should be \"AGE_RENAMED\"", () => propertyShouldBe(page, "yColumnName", el("scatter plot viewer"), "AGE_RENAMED"));
    await session.step(32, "And \"formulaLines\" property of scatter plot viewer should contain \"${AGE_RENAMED} = ${AGE_RENAMED}\"", () => propertyShouldContain(page, "formulaLines", el("scatter plot viewer"), "${AGE_RENAMED} = ${AGE_RENAMED}"));
    await session.step(33, "And the \"formula lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "formula lines", el("scatter plot viewer"), 1));
    await session.step(34, "And no errors should have been logged", () => noErrors(page));
    await session.step(35, "When user renames \"AGE_RENAMED\" column to \"AGE\"", () => renameColumn(page, "AGE_RENAMED", "AGE"));
    await session.step(36, "Then \"xColumnName\" property of scatter plot viewer should be \"AGE\"", () => propertyShouldBe(page, "xColumnName", el("scatter plot viewer"), "AGE"));
    await session.step(37, "And \"formulaLines\" property of scatter plot viewer should contain \"${AGE} = ${AGE}\"", () => propertyShouldContain(page, "formulaLines", el("scatter plot viewer"), "${AGE} = ${AGE}"));
    await session.step(38, "And the \"formula lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "formula lines", el("scatter plot viewer"), 1));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
  });
  test("Renaming a column rewrites the scatter plot line that uses it", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.line-chart", "@realizes:viewers.trellis-plot", "@realizes:powerpack.dialogs.formula-lines"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(42, "Given user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE"],["yColumnName","HEIGHT"],["formulaLines","[{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${AGE}\"}]"]]), [["xColumnName","AGE"],["yColumnName","HEIGHT"],["formulaLines","[{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${AGE}\"}]"]]);
    await session.step(46, "Then the \"formula lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "formula lines", el("scatter plot viewer"), 1));
    await session.step(47, "When user renames \"AGE\" column to \"AGE_R\"", () => renameColumn(page, "AGE", "AGE_R"));
    await session.step(48, "Then \"formulaLines\" property of scatter plot viewer should contain \"${HEIGHT} = ${AGE_R}\"", () => propertyShouldContain(page, "formulaLines", el("scatter plot viewer"), "${HEIGHT} = ${AGE_R}"));
    await session.step(49, "And the \"formula lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "formula lines", el("scatter plot viewer"), 1));
    await session.step(50, "When user renames \"AGE_R\" column to \"AGE\"", () => renameColumn(page, "AGE_R", "AGE"));
    await session.step(51, "Then \"formulaLines\" property of scatter plot viewer should contain \"${HEIGHT} = ${AGE}\"", () => propertyShouldContain(page, "formulaLines", el("scatter plot viewer"), "${HEIGHT} = ${AGE}"));
    await session.step(52, "And no errors should have been logged", () => noErrors(page));
  });
  test("Renaming a column rewrites the line chart line added from its axis", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.line-chart", "@realizes:viewers.trellis-plot", "@realizes:powerpack.dialogs.formula-lines"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(55, "Given user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","AGE"],["yColumnNames","HEIGHT"]]), [["xColumnName","AGE"],["yColumnNames","HEIGHT"]]);
    await session.step(58, "And user resizes line chart viewer to 800 by 500", () => resizeTo(page, el("line chart viewer"), 800, 500));
    await session.step(59, "When user picks \"Annotations > Add Line\" from the context menu of the \"bottom edge of y axis\" area of line chart viewer", () => pickFromAreaContextMenu(page, "Annotations > Add Line", "bottom edge of y axis", el("line chart viewer")));
    await session.step(60, "And user clicks OK button in \"Formula Lines\" dialog", () => clickOn(page, el("OK button in \"Formula Lines\" dialog")));
    await session.step(61, "Then \"formulaLines\" property of line chart viewer should contain \"${avg(HEIGHT)} = 168.8\"", () => propertyShouldContain(page, "formulaLines", el("line chart viewer"), "${avg(HEIGHT)} = 168.8"));
    await session.step(62, "And line chart viewer should have a \"formula line avg(HEIGHT) = 168.8\" area", () => hasArea(page, el("line chart viewer"), "formula line avg(HEIGHT) = 168.8"));
    await session.step(63, "When user renames \"HEIGHT\" column to \"HEIGHT_R\"", () => renameColumn(page, "HEIGHT", "HEIGHT_R"));
    await session.step(64, "Then \"formulaLines\" property of line chart viewer should contain \"${avg(HEIGHT_R)} = 168.8\"", () => propertyShouldContain(page, "formulaLines", el("line chart viewer"), "${avg(HEIGHT_R)} = 168.8"));
    await session.step(65, "And line chart viewer should have a \"formula line avg(HEIGHT_R) = 168.8\" area", () => hasArea(page, el("line chart viewer"), "formula line avg(HEIGHT_R) = 168.8"));
    await session.step(66, "When user renames \"HEIGHT_R\" column to \"HEIGHT\"", () => renameColumn(page, "HEIGHT_R", "HEIGHT"));
    await session.step(67, "Then \"formulaLines\" property of line chart viewer should contain \"${avg(HEIGHT)} = 168.8\"", () => propertyShouldContain(page, "formulaLines", el("line chart viewer"), "${avg(HEIGHT)} = 168.8"));
    await session.step(68, "And no errors should have been logged", () => noErrors(page));
  });
  test("A line on the columns survives a change of the axis columns untouched, a line on the axis tokens follows the axes (public#3177)", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.line-chart", "@realizes:viewers.trellis-plot", "@realizes:powerpack.dialogs.formula-lines", "@realizes:public#3177"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(72, "Given user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"],["formulaLines","[{\"type\":\"line\",\"formula\":\"${Y} = ${X}\"},{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${WEIGHT} + 90\"}]"]]), [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"],["formulaLines","[{\"type\":\"line\",\"formula\":\"${Y} = ${X}\"},{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${WEIGHT} + 90\"}]"]]);
    await session.step(76, "And user resizes scatter plot viewer to 800 by 500", () => resizeTo(page, el("scatter plot viewer"), 800, 500));
    await session.step(77, "Then scatter plot viewer should draw 2 formula lines", () => drawsFormulaLines(page, el("scatter plot viewer"), 2));
    await session.step(78, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["xColumnName","AGE"],["yColumnName","WEIGHT"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"]]);
    await session.step(81, "Then \"formulaLines\" property of scatter plot viewer should be '[{\"type\":\"line\",\"formula\":\"${Y} = ${X}\"},{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${WEIGHT} + 90\"}]'", () => propertyShouldBe(page, "formulaLines", el("scatter plot viewer"), "[{\"type\":\"line\",\"formula\":\"${Y} = ${X}\"},{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${WEIGHT} + 90\"}]"));
    await session.step(82, "And scatter plot viewer should draw 1 formula line", () => drawsFormulaLines(page, el("scatter plot viewer"), 1));
    await session.step(83, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"]]), [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"]]);
    await session.step(86, "Then scatter plot viewer should draw 2 formula lines", () => drawsFormulaLines(page, el("scatter plot viewer"), 2));
    await session.step(87, "And no errors should have been logged", () => noErrors(page));
  });
  test("A horizontal band survives a logarithmic Y axis on the scatter plot", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.line-chart", "@realizes:viewers.trellis-plot", "@realizes:powerpack.dialogs.formula-lines"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(90, "Given user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"],["formulaLines","[{\"type\":\"band\",\"formula\":\"${HEIGHT} in (160.9, 177.6)\",\"orientation\":\"Horizontal\",\"column2\":\"WEIGHT\"}]"]]), [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"],["formulaLines","[{\"type\":\"band\",\"formula\":\"${HEIGHT} in (160.9, 177.6)\",\"orientation\":\"Horizontal\",\"column2\":\"WEIGHT\"}]"]]);
    await session.step(94, "Then the \"formula lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "formula lines", el("scatter plot viewer"), 1));
    await session.step(95, "When user sets \"yAxisType\" property of scatter plot viewer to \"logarithmic\"", () => setProperty(page, "yAxisType", el("scatter plot viewer"), "logarithmic"));
    await session.step(96, "Then the \"formula lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "formula lines", el("scatter plot viewer"), 1));
    await session.step(97, "And scatter plot viewer should have a \"formula band 1\" area", () => hasArea(page, el("scatter plot viewer"), "formula band 1"));
    await session.step(98, "And scatter plot viewer should have repainted", () => repainted(page, el("scatter plot viewer")));
    await session.step(99, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(100, "When user sets \"yAxisType\" property of scatter plot viewer to \"linear\"", () => setProperty(page, "yAxisType", el("scatter plot viewer"), "linear"));
    await session.step(101, "Then the \"formula lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "formula lines", el("scatter plot viewer"), 1));
    await session.step(102, "And no errors should have been logged", () => noErrors(page));
  });
  test("A line survives a logarithmic Y axis on the line chart", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.line-chart", "@realizes:viewers.trellis-plot", "@realizes:powerpack.dialogs.formula-lines"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(105, "Given user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","AGE"],["yColumnNames","HEIGHT"]]), [["xColumnName","AGE"],["yColumnNames","HEIGHT"]]);
    await session.step(108, "And user resizes line chart viewer to 800 by 500", () => resizeTo(page, el("line chart viewer"), 800, 500));
    await session.step(109, "When user picks \"Annotations > Add Line\" from the context menu of the \"bottom edge of y axis\" area of line chart viewer", () => pickFromAreaContextMenu(page, "Annotations > Add Line", "bottom edge of y axis", el("line chart viewer")));
    await session.step(110, "And user clicks OK button in \"Formula Lines\" dialog", () => clickOn(page, el("OK button in \"Formula Lines\" dialog")));
    await session.step(111, "Then line chart viewer should have a \"formula line avg(HEIGHT) = 168.8\" area", () => hasArea(page, el("line chart viewer"), "formula line avg(HEIGHT) = 168.8"));
    await session.step(112, "When user sets \"yAxisType\" property of line chart viewer to \"logarithmic\"", () => setProperty(page, "yAxisType", el("line chart viewer"), "logarithmic"));
    await session.step(113, "Then line chart viewer should have a \"formula line avg(HEIGHT) = 168.8\" area", () => hasArea(page, el("line chart viewer"), "formula line avg(HEIGHT) = 168.8"));
    await session.step(114, "And the \"formula lines\" reading of line chart viewer should be 1", () => readingIs(page, "formula lines", el("line chart viewer"), 1));
    await session.step(115, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(116, "When user sets \"yAxisType\" property of line chart viewer to \"linear\"", () => setProperty(page, "yAxisType", el("line chart viewer"), "linear"));
    await session.step(117, "Then line chart viewer should have a \"formula line avg(HEIGHT) = 168.8\" area", () => hasArea(page, el("line chart viewer"), "formula line avg(HEIGHT) = 168.8"));
    await session.step(118, "And no errors should have been logged", () => noErrors(page));
  });
  test("Hovering markers next to a formula line raises no errors", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.line-chart", "@realizes:viewers.trellis-plot", "@realizes:powerpack.dialogs.formula-lines"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(121, "Given user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"],["formulaLines","[{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${WEIGHT} + 100\"}]"]]), [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"],["formulaLines","[{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${WEIGHT} + 100\"}]"]]);
    await session.step(125, "And user resizes scatter plot viewer to 800 by 500", () => resizeTo(page, el("scatter plot viewer"), 800, 500));
    await session.step(126, "Then the \"formula lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "formula lines", el("scatter plot viewer"), 1));
    await session.step(127, "And scatter plot viewer should draw 1 formula line", () => drawsFormulaLines(page, el("scatter plot viewer"), 1));
    await session.step(128, "When user hovers over the \"marker of row 1\" area of scatter plot viewer", () => hoverArea(page, "marker of row 1", el("scatter plot viewer")));
    await session.step(129, "Then the tooltip should show some columns", () => tooltipSomeColumns(page));
    await session.step(130, "When user hovers over the \"marker of row 2\" area of scatter plot viewer", () => hoverArea(page, "marker of row 2", el("scatter plot viewer")));
    await session.step(131, "And user hovers over the \"marker of row 3\" area of scatter plot viewer", () => hoverArea(page, "marker of row 3", el("scatter plot viewer")));
    await session.step(132, "And user hovers over the \"marker of row 4\" area of scatter plot viewer", () => hoverArea(page, "marker of row 4", el("scatter plot viewer")));
    await session.step(133, "And user hovers over the \"marker of row 5\" area of scatter plot viewer", () => hoverArea(page, "marker of row 5", el("scatter plot viewer")));
    await session.step(134, "And user hovers over the \"marker of row 6\" area of scatter plot viewer", () => hoverArea(page, "marker of row 6", el("scatter plot viewer")));
    await session.step(135, "And user hovers over the \"marker of row 7\" area of scatter plot viewer", () => hoverArea(page, "marker of row 7", el("scatter plot viewer")));
    await session.step(136, "And user hovers over the \"marker of row 8\" area of scatter plot viewer", () => hoverArea(page, "marker of row 8", el("scatter plot viewer")));
    await session.step(137, "Then the tooltip should show some columns", () => tooltipSomeColumns(page));
    await session.step(138, "And exactly one tooltip should be shown", () => oneTooltip(page));
    await session.step(139, "When user moves the pointer away from scatter plot viewer", () => pointerAway(page, el("scatter plot viewer")));
    await session.step(140, "Then no errors should have been logged", () => noErrors(page));
  });
  test("A band on a column with parentheses in its name is drawn (public#542)", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.line-chart", "@realizes:viewers.trellis-plot", "@realizes:powerpack.dialogs.formula-lines", "@realizes:public#542"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(144, "Given user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"]]), [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"]]);
    await session.step(147, "And user resizes scatter plot viewer to 800 by 500", () => resizeTo(page, el("scatter plot viewer"), 800, 500));
    await session.step(148, "When user renames \"WEIGHT\" column to \"WEIGHT (kg)\"", () => renameColumn(page, "WEIGHT", "WEIGHT (kg)"));
    await session.step(149, "And user picks \"Annotations > Add Band\" from the context menu of the \"left edge of x axis\" area of scatter plot viewer", () => pickFromAreaContextMenu(page, "Annotations > Add Band", "left edge of x axis", el("scatter plot viewer")));
    await session.step(150, "And user clicks OK button in \"Formula Lines\" dialog", () => clickOn(page, el("OK button in \"Formula Lines\" dialog")));
    await session.step(151, "Then \"formulaLines\" property of scatter plot viewer should contain \"${WEIGHT (kg)} in (64.2, 91.0)\"", () => propertyShouldContain(page, "formulaLines", el("scatter plot viewer"), "${WEIGHT (kg)} in (64.2, 91.0)"));
    await session.step(152, "And scatter plot viewer should draw 1 formula line", () => drawsFormulaLines(page, el("scatter plot viewer"), 1));
    await session.step(153, "And no errors should have been logged", () => noErrors(page));
  });
  test("A valid line is drawn next to an invalid one, also after a layout round trip (public#699)", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.line-chart", "@realizes:viewers.trellis-plot", "@realizes:powerpack.dialogs.formula-lines", "@realizes:public#699"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(157, "Given user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"],["formulaLines","[{\"type\":\"line\",\"formula\":\"${HEIGHT} = 170\"},{\"type\":\"line\",\"formula\":\"${WEIGHT} + 1\"}]"]]), [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"],["formulaLines","[{\"type\":\"line\",\"formula\":\"${HEIGHT} = 170\"},{\"type\":\"line\",\"formula\":\"${WEIGHT} + 1\"}]"]]);
    await session.step(161, "And user resizes scatter plot viewer to 800 by 500", () => resizeTo(page, el("scatter plot viewer"), 800, 500));
    await session.step(162, "Then scatter plot viewer should draw 1 formula line", () => drawsFormulaLines(page, el("scatter plot viewer"), 1));
    await session.step(163, "And the \"HEIGHT\" line of scatter plot viewer should be drawn at its value on the y axis", () => lineDrawnAtValue(page, "HEIGHT", el("scatter plot viewer"), "y"));
    await session.step(164, "When user saves the layout of the current table view", () => saveLayout(page));
    await session.step(165, "And user loads the saved layout", () => loadLayout(page));
    await session.step(166, "Then scatter plot viewer should draw 1 formula line", () => drawsFormulaLines(page, el("scatter plot viewer"), 1));
    await session.step(167, "And the \"HEIGHT\" line of scatter plot viewer should be drawn at its value on the y axis", () => lineDrawnAtValue(page, "HEIGHT", el("scatter plot viewer"), "y"));
    await session.step(168, "When user picks \"Tools > Formula Lines...\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Tools > Formula Lines...", el("scatter plot viewer")));
    await session.step(169, "And user clicks on the \"cell 1 of title\" area of grid in \"Formula Lines\" dialog", () => clickArea(page, "cell 1 of title", el("grid in \"Formula Lines\" dialog")));
    await session.step(170, "Then the \"current item\" reading of scatter plot viewer in \"Formula Lines\" dialog should be \"${HEIGHT} = 170\"", () => readingReads(page, "current item", el("scatter plot viewer in \"Formula Lines\" dialog"), "${HEIGHT} = 170"));
    await session.step(171, "And scatter plot viewer in \"Formula Lines\" dialog should draw 1 formula line", () => drawsFormulaLines(page, el("scatter plot viewer in \"Formula Lines\" dialog"), 1));
    await session.step(172, "When user clicks CANCEL button in \"Formula Lines\" dialog", () => clickOn(page, el("CANCEL button in \"Formula Lines\" dialog")));
    await session.step(173, "Then no errors should have been logged", () => noErrors(page));
  });
  test("The Formula Lines dialog opened from a trellis cell previews the inner scatter plot's X and Y (public#3390, public#3389)", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.line-chart", "@realizes:viewers.trellis-plot", "@realizes:powerpack.dialogs.formula-lines", "@realizes:public#3390", "@realizes:public#3389"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(177, "Given user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Viewer Type","Scatter plot"]]), [["X Column Names","SEX"],["Viewer Type","Scatter plot"]]);
    await session.step(180, "And user sets \"xColumnName\" inner property of trellis plot viewer to \"WEIGHT\"", () => setInnerProperty(page, "xColumnName", el("trellis plot viewer"), "WEIGHT"));
    await session.step(181, "And user sets \"yColumnName\" inner property of trellis plot viewer to \"HEIGHT\"", () => setInnerProperty(page, "yColumnName", el("trellis plot viewer"), "HEIGHT"));
    await session.step(182, "When user picks \"Scatter plot > Tools > Formula Lines...\" from the context menu of the \"cell body F | Critical\" area of trellis plot viewer", () => pickFromAreaContextMenu(page, "Scatter plot > Tools > Formula Lines...", "cell body F | Critical", el("trellis plot viewer")));
    await session.step(183, "Then \"Formula Lines\" dialog should be visible", () => shouldBe(page, el("\"Formula Lines\" dialog"), "visible"));
    await session.step(184, "When user clicks on \"ADD NEW\" button in \"Formula Lines\" dialog", () => clickOn(page, el("\"ADD NEW\" button in \"Formula Lines\" dialog")));
    await session.step(185, "And user picks \"Line - Vertical\" from the open menu", () => pickFromOpenMenu(page, "Line - Vertical"));
    await session.step(186, "Then the \"current item\" reading of scatter plot viewer in \"Formula Lines\" dialog should be \"${WEIGHT} = 77.4\"", () => readingReads(page, "current item", el("scatter plot viewer in \"Formula Lines\" dialog"), "${WEIGHT} = 77.4"));
    await session.step(187, "And \"xColumnName\" property of scatter plot viewer in \"Formula Lines\" dialog should be \"WEIGHT\"", () => propertyShouldBe(page, "xColumnName", el("scatter plot viewer in \"Formula Lines\" dialog"), "WEIGHT"));
    await session.step(188, "And \"yColumnName\" property of scatter plot viewer in \"Formula Lines\" dialog should be \"HEIGHT\"", () => propertyShouldBe(page, "yColumnName", el("scatter plot viewer in \"Formula Lines\" dialog"), "HEIGHT"));
    await session.step(189, "And the \"y axis min\" reading of scatter plot viewer in \"Formula Lines\" dialog should be between 134.24 and 134.25", () => readingBetween(page, "y axis min", el("scatter plot viewer in \"Formula Lines\" dialog"), 134.24, 134.25));
    await session.step(190, "When user clicks OK button in \"Formula Lines\" dialog", () => clickOn(page, el("OK button in \"Formula Lines\" dialog")));
    await session.step(191, "And user picks \"Scatter plot > Tools > Formula Lines...\" from the context menu of the \"cell body F | Critical\" area of trellis plot viewer", () => pickFromAreaContextMenu(page, "Scatter plot > Tools > Formula Lines...", "cell body F | Critical", el("trellis plot viewer")));
    await session.step(192, "Then the \"current item\" reading of scatter plot viewer in \"Formula Lines\" dialog should be \"${WEIGHT} = 77.4\"", () => readingReads(page, "current item", el("scatter plot viewer in \"Formula Lines\" dialog"), "${WEIGHT} = 77.4"));
    await session.step(193, "When user clicks on \"ADD NEW\" button in \"Formula Lines\" dialog", () => clickOn(page, el("\"ADD NEW\" button in \"Formula Lines\" dialog")));
    await session.step(194, "And user picks \"Band - Vertical\" from the open menu", () => pickFromOpenMenu(page, "Band - Vertical"));
    await session.step(195, "Then the \"current item\" reading of scatter plot viewer in \"Formula Lines\" dialog should be \"${WEIGHT} in(64.2, 91.0)\"", () => readingReads(page, "current item", el("scatter plot viewer in \"Formula Lines\" dialog"), "${WEIGHT} in(64.2, 91.0)"));
    await session.step(196, "And \"xColumnName\" property of scatter plot viewer in \"Formula Lines\" dialog should be \"WEIGHT\"", () => propertyShouldBe(page, "xColumnName", el("scatter plot viewer in \"Formula Lines\" dialog"), "WEIGHT"));
    await session.step(197, "And \"yColumnName\" property of scatter plot viewer in \"Formula Lines\" dialog should be \"HEIGHT\"", () => propertyShouldBe(page, "yColumnName", el("scatter plot viewer in \"Formula Lines\" dialog"), "HEIGHT"));
    await session.step(198, "When user clicks CANCEL button in \"Formula Lines\" dialog", () => clickOn(page, el("CANCEL button in \"Formula Lines\" dialog")));
    await session.step(199, "Then no errors should have been logged", () => noErrors(page));
  });
  test("Formula Lines stays in the context menu after the axis type is switched from the axis menu (public#425)", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.line-chart", "@realizes:viewers.trellis-plot", "@realizes:powerpack.dialogs.formula-lines", "@realizes:public#425"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(203, "Given user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"]]), [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"]]);
    await session.step(206, "And user resizes scatter plot viewer to 800 by 500", () => resizeTo(page, el("scatter plot viewer"), 800, 500));
    await session.step(207, "When user picks \"X Axis Type > logarithmic\" from the context menu of the \"left edge of x axis\" area of scatter plot viewer", () => pickFromAreaContextMenu(page, "X Axis Type > logarithmic", "left edge of x axis", el("scatter plot viewer")));
    await session.step(208, "Then \"xAxisType\" property of scatter plot viewer should be \"logarithmic\"", () => propertyShouldBe(page, "xAxisType", el("scatter plot viewer"), "logarithmic"));
    await session.step(209, "When user picks \"Tools > Formula Lines...\" from the context menu of the \"view\" area of scatter plot viewer", () => pickFromAreaContextMenu(page, "Tools > Formula Lines...", "view", el("scatter plot viewer")));
    await session.step(210, "Then \"Formula Lines\" dialog should be visible", () => shouldBe(page, el("\"Formula Lines\" dialog"), "visible"));
    await session.step(211, "When user clicks CANCEL button in \"Formula Lines\" dialog", () => clickOn(page, el("CANCEL button in \"Formula Lines\" dialog")));
    await session.step(212, "Then no errors should have been logged", () => noErrors(page));
  });
});
