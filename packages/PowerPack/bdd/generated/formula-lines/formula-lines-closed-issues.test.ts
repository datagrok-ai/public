/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/formula-lines/formula-lines-closed-issues.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot, viewers.trellis-plot, powerpack.dialogs.formula-lines, public#3177, public#3019, public#542, public#699, public#3390, public#3389, public#425]
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
import {closeAllViews, openDataset, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {drawsFormulaLines, lineDrawnAtValue} from '@datagrok-libraries/bdd/bindings/tiers/viewers/formula-lines';
import {addViewerWith, clickArea, loadLayout, noErrors, pickFromAreaContextMenu, pickFromContextMenu, pickFromOpenMenu, propertyShouldBe, propertyShouldContain, readingAtLeast, readingReads, resizeTo, saveLayout, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {setInnerProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Formula lines — closed GitHub issues guarded on the scatter plot", () => {
  const session = feature(test, "features/formula-lines/formula-lines-closed-issues.feature", import.meta.url);
  test("A line on the axis tokens follows the columns on the axes, a line on column names does not", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.trellis-plot", "@realizes:powerpack.dialogs.formula-lines", "@realizes:public#3177"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(24, "Given user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"],["formulaLines","[{\"type\":\"line\",\"formula\":\"${Y} = ${X}\"},{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${WEIGHT} + 90\"}]"]]), [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"],["formulaLines","[{\"type\":\"line\",\"formula\":\"${Y} = ${X}\"},{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${WEIGHT} + 90\"}]"]]);
    await session.step(28, "And user resizes scatter plot viewer to 800 by 500", () => resizeTo(page, el("scatter plot viewer"), 800, 500));
    await session.step(29, "Then scatter plot viewer should draw 2 formula lines", () => drawsFormulaLines(page, el("scatter plot viewer"), 2));
    await session.step(30, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["xColumnName","AGE"],["yColumnName","WEIGHT"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"]]);
    await session.step(33, "Then scatter plot viewer should draw 1 formula line", () => drawsFormulaLines(page, el("scatter plot viewer"), 1));
    await session.step(34, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"]]), [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"]]);
    await session.step(37, "Then scatter plot viewer should draw 2 formula lines", () => drawsFormulaLines(page, el("scatter plot viewer"), 2));
    await session.step(38, "And no errors should have been logged", () => noErrors(page));
  });
  test("Viewer formula lines on two scatter plots come back with the project", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.trellis-plot", "@realizes:powerpack.dialogs.formula-lines", "@realizes:public#3019"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(43, "Given user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"]]), [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"]]);
    await session.step(46, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"]]);
    await session.step(49, "When user sets \"formulaLines\" property of first scatter plot viewer to '[{\"type\":\"line\",\"formula\":\"${HEIGHT} = 170\"}]'", () => setProperty(page, "formulaLines", el("first scatter plot viewer"), "[{\"type\":\"line\",\"formula\":\"${HEIGHT} = 170\"}]"));
    await session.step(50, "And user sets \"formulaLines\" property of second scatter plot viewer to '[{\"type\":\"line\",\"formula\":\"${WEIGHT} = 100\"},{\"type\":\"line\",\"formula\":\"${AGE} = 50\"}]'", () => setProperty(page, "formulaLines", el("second scatter plot viewer"), "[{\"type\":\"line\",\"formula\":\"${WEIGHT} = 100\"},{\"type\":\"line\",\"formula\":\"${AGE} = 50\"}]"));
    await session.step(51, "Then first scatter plot viewer should draw 1 formula line", () => drawsFormulaLines(page, el("first scatter plot viewer"), 1));
    await session.step(52, "And second scatter plot viewer should draw 2 formula lines", () => drawsFormulaLines(page, el("second scatter plot viewer"), 2));
    await session.step(53, "When user saves the current view as project \"bdd-fl-two-plots-{time}\"", () => saveAsProject(page, session.text("bdd-fl-two-plots-{time}")));
    await session.step(54, "And user closes all views", () => closeAllViews(page));
    await session.step(55, "And user opens the \"bdd-fl-two-plots-{time}\" project", () => openProject(page, session.text("bdd-fl-two-plots-{time}")));
    await session.step(56, "Then first scatter plot viewer should draw 1 formula line", () => drawsFormulaLines(page, el("first scatter plot viewer"), 1));
    await session.step(57, "And second scatter plot viewer should draw 2 formula lines", () => drawsFormulaLines(page, el("second scatter plot viewer"), 2));
    await session.step(58, "And the \"HEIGHT\" line of first scatter plot viewer should be drawn at its value on the y axis", () => lineDrawnAtValue(page, "HEIGHT", el("first scatter plot viewer"), "y"));
    await session.step(59, "And no errors should have been logged", () => noErrors(page));
  });
  test("Bands on columns with parentheses in their names are drawn", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.trellis-plot", "@realizes:powerpack.dialogs.formula-lines", "@realizes:public#542"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(64, "Given user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"]]), [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"]]);
    await session.step(67, "And user resizes scatter plot viewer to 800 by 500", () => resizeTo(page, el("scatter plot viewer"), 800, 500));
    await session.step(68, "When user renames \"WEIGHT\" column to \"WEIGHT (kg)\"", () => renameColumn(page, "WEIGHT", "WEIGHT (kg)"));
    await session.step(69, "And user renames \"HEIGHT\" column to \"HEIGHT (cm)\"", () => renameColumn(page, "HEIGHT", "HEIGHT (cm)"));
    await session.step(70, "And user picks \"Annotations > Add Band\" from the context menu of the \"left edge of x axis\" area of scatter plot viewer", () => pickFromAreaContextMenu(page, "Annotations > Add Band", "left edge of x axis", el("scatter plot viewer")));
    await session.step(71, "And user clicks OK button in \"Formula Lines\" dialog", () => clickOn(page, el("OK button in \"Formula Lines\" dialog")));
    await session.step(72, "Then \"formulaLines\" property of scatter plot viewer should contain \"${WEIGHT (kg)} in (64.2, 91.0)\"", () => propertyShouldContain(page, "formulaLines", el("scatter plot viewer"), "${WEIGHT (kg)} in (64.2, 91.0)"));
    await session.step(73, "And scatter plot viewer should draw 1 formula line", () => drawsFormulaLines(page, el("scatter plot viewer"), 1));
    await session.step(74, "When user picks \"Annotations > Add Band\" from the context menu of the \"bottom edge of y axis\" area of scatter plot viewer", () => pickFromAreaContextMenu(page, "Annotations > Add Band", "bottom edge of y axis", el("scatter plot viewer")));
    await session.step(75, "And user clicks OK button in \"Formula Lines\" dialog", () => clickOn(page, el("OK button in \"Formula Lines\" dialog")));
    await session.step(76, "Then \"formulaLines\" property of scatter plot viewer should contain \"${HEIGHT (cm)} in (160.9, 177.6)\"", () => propertyShouldContain(page, "formulaLines", el("scatter plot viewer"), "${HEIGHT (cm)} in (160.9, 177.6)"));
    await session.step(77, "And scatter plot viewer should draw 2 formula lines", () => drawsFormulaLines(page, el("scatter plot viewer"), 2));
    await session.step(78, "When user sets \"formulaLines\" property of scatter plot viewer to \"\"", () => setProperty(page, "formulaLines", el("scatter plot viewer"), ""));
    await session.step(79, "And user renames \"WEIGHT (kg)\" column to \"WEIGHT\"", () => renameColumn(page, "WEIGHT (kg)", "WEIGHT"));
    await session.step(80, "And user renames \"HEIGHT (cm)\" column to \"HEIGHT\"", () => renameColumn(page, "HEIGHT (cm)", "HEIGHT"));
    await session.step(81, "Then no errors should have been logged", () => noErrors(page));
  });
  test("A valid line is drawn next to an invalid one, also after a layout round trip", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.trellis-plot", "@realizes:powerpack.dialogs.formula-lines", "@realizes:public#699"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(86, "Given user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"],["formulaLines","[{\"type\":\"line\",\"formula\":\"${HEIGHT} = 170\"},{\"type\":\"line\",\"formula\":\"${WEIGHT} + 1\"}]"]]), [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"],["formulaLines","[{\"type\":\"line\",\"formula\":\"${HEIGHT} = 170\"},{\"type\":\"line\",\"formula\":\"${WEIGHT} + 1\"}]"]]);
    await session.step(90, "And user resizes scatter plot viewer to 800 by 500", () => resizeTo(page, el("scatter plot viewer"), 800, 500));
    await session.step(91, "Then scatter plot viewer should draw 1 formula line", () => drawsFormulaLines(page, el("scatter plot viewer"), 1));
    await session.step(92, "And the \"HEIGHT\" line of scatter plot viewer should be drawn at its value on the y axis", () => lineDrawnAtValue(page, "HEIGHT", el("scatter plot viewer"), "y"));
    await session.step(93, "When user saves the layout of the current table view", () => saveLayout(page));
    await session.step(94, "And user loads the saved layout", () => loadLayout(page));
    await session.step(95, "Then scatter plot viewer should draw 1 formula line", () => drawsFormulaLines(page, el("scatter plot viewer"), 1));
    await session.step(96, "And the \"HEIGHT\" line of scatter plot viewer should be drawn at its value on the y axis", () => lineDrawnAtValue(page, "HEIGHT", el("scatter plot viewer"), "y"));
    await session.step(97, "When user picks \"Tools > Formula Lines...\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Tools > Formula Lines...", el("scatter plot viewer")));
    await session.step(98, "And user clicks on the \"cell 1 of title\" area of grid in \"Formula Lines\" dialog", () => clickArea(page, "cell 1 of title", el("grid in \"Formula Lines\" dialog")));
    await session.step(99, "Then the \"current item\" reading of scatter plot viewer in \"Formula Lines\" dialog should be \"${HEIGHT} = 170\"", () => readingReads(page, "current item", el("scatter plot viewer in \"Formula Lines\" dialog"), "${HEIGHT} = 170"));
    await session.step(100, "And scatter plot viewer in \"Formula Lines\" dialog should draw 1 formula line", () => drawsFormulaLines(page, el("scatter plot viewer in \"Formula Lines\" dialog"), 1));
    await session.step(101, "When user clicks CANCEL button in \"Formula Lines\" dialog", () => clickOn(page, el("CANCEL button in \"Formula Lines\" dialog")));
    await session.step(102, "Then no errors should have been logged", () => noErrors(page));
  });
  test("The Formula Lines dialog opened from a trellis cell previews the inner scatter plot's X and Y", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.trellis-plot", "@realizes:powerpack.dialogs.formula-lines", "@realizes:public#3390", "@realizes:public#3389"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(107, "Given user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Viewer Type","Scatter plot"]]), [["X Column Names","SEX"],["Viewer Type","Scatter plot"]]);
    await session.step(110, "And user sets \"xColumnName\" inner property of trellis plot viewer to \"WEIGHT\"", () => setInnerProperty(page, "xColumnName", el("trellis plot viewer"), "WEIGHT"));
    await session.step(111, "And user sets \"yColumnName\" inner property of trellis plot viewer to \"HEIGHT\"", () => setInnerProperty(page, "yColumnName", el("trellis plot viewer"), "HEIGHT"));
    await session.step(112, "When user picks \"Scatter plot > Tools > Formula Lines...\" from the context menu of the \"cell body F | Critical\" area of trellis plot viewer", () => pickFromAreaContextMenu(page, "Scatter plot > Tools > Formula Lines...", "cell body F | Critical", el("trellis plot viewer")));
    await session.step(113, "Then \"Formula Lines\" dialog should be visible", () => shouldBe(page, el("\"Formula Lines\" dialog"), "visible"));
    await session.step(114, "When user clicks on \"ADD NEW\" button in \"Formula Lines\" dialog", () => clickOn(page, el("\"ADD NEW\" button in \"Formula Lines\" dialog")));
    await session.step(115, "And user picks \"Line - Vertical\" from the open menu", () => pickFromOpenMenu(page, "Line - Vertical"));
    await session.step(116, "Then the \"current item\" reading of scatter plot viewer in \"Formula Lines\" dialog should be \"${WEIGHT} = 77.4\"", () => readingReads(page, "current item", el("scatter plot viewer in \"Formula Lines\" dialog"), "${WEIGHT} = 77.4"));
    await session.step(117, "And \"xColumnName\" property of scatter plot viewer in \"Formula Lines\" dialog should be \"WEIGHT\"", () => propertyShouldBe(page, "xColumnName", el("scatter plot viewer in \"Formula Lines\" dialog"), "WEIGHT"));
    await session.step(118, "And \"yColumnName\" property of scatter plot viewer in \"Formula Lines\" dialog should be \"HEIGHT\"", () => propertyShouldBe(page, "yColumnName", el("scatter plot viewer in \"Formula Lines\" dialog"), "HEIGHT"));
    await session.step(119, "And the \"y axis min\" reading of scatter plot viewer in \"Formula Lines\" dialog should be at least 100", () => readingAtLeast(page, "y axis min", el("scatter plot viewer in \"Formula Lines\" dialog"), 100));
    await session.step(120, "When user clicks OK button in \"Formula Lines\" dialog", () => clickOn(page, el("OK button in \"Formula Lines\" dialog")));
    await session.step(121, "And user picks \"Scatter plot > Tools > Formula Lines...\" from the context menu of the \"cell body F | Critical\" area of trellis plot viewer", () => pickFromAreaContextMenu(page, "Scatter plot > Tools > Formula Lines...", "cell body F | Critical", el("trellis plot viewer")));
    await session.step(122, "Then the \"current item\" reading of scatter plot viewer in \"Formula Lines\" dialog should be \"${WEIGHT} = 77.4\"", () => readingReads(page, "current item", el("scatter plot viewer in \"Formula Lines\" dialog"), "${WEIGHT} = 77.4"));
    await session.step(123, "When user clicks on \"ADD NEW\" button in \"Formula Lines\" dialog", () => clickOn(page, el("\"ADD NEW\" button in \"Formula Lines\" dialog")));
    await session.step(124, "And user picks \"Band - Vertical\" from the open menu", () => pickFromOpenMenu(page, "Band - Vertical"));
    await session.step(125, "Then the \"current item\" reading of scatter plot viewer in \"Formula Lines\" dialog should be \"${WEIGHT} in(64.2, 91.0)\"", () => readingReads(page, "current item", el("scatter plot viewer in \"Formula Lines\" dialog"), "${WEIGHT} in(64.2, 91.0)"));
    await session.step(126, "And \"xColumnName\" property of scatter plot viewer in \"Formula Lines\" dialog should be \"WEIGHT\"", () => propertyShouldBe(page, "xColumnName", el("scatter plot viewer in \"Formula Lines\" dialog"), "WEIGHT"));
    await session.step(127, "And \"yColumnName\" property of scatter plot viewer in \"Formula Lines\" dialog should be \"HEIGHT\"", () => propertyShouldBe(page, "yColumnName", el("scatter plot viewer in \"Formula Lines\" dialog"), "HEIGHT"));
    await session.step(128, "When user clicks CANCEL button in \"Formula Lines\" dialog", () => clickOn(page, el("CANCEL button in \"Formula Lines\" dialog")));
    await session.step(129, "And user sets \"formulaLines\" inner property of trellis plot viewer to \"\"", () => setInnerProperty(page, "formulaLines", el("trellis plot viewer"), ""));
    await session.step(130, "Then no errors should have been logged", () => noErrors(page));
  });
  test("Formula Lines stays in the context menu after the axis type is switched from the axis menu", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.trellis-plot", "@realizes:powerpack.dialogs.formula-lines", "@realizes:public#425"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(135, "Given user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"]]), [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"]]);
    await session.step(138, "And user resizes scatter plot viewer to 800 by 500", () => resizeTo(page, el("scatter plot viewer"), 800, 500));
    await session.step(139, "When user picks \"X Axis Type > logarithmic\" from the context menu of the \"left edge of x axis\" area of scatter plot viewer", () => pickFromAreaContextMenu(page, "X Axis Type > logarithmic", "left edge of x axis", el("scatter plot viewer")));
    await session.step(140, "Then \"xAxisType\" property of scatter plot viewer should be \"logarithmic\"", () => propertyShouldBe(page, "xAxisType", el("scatter plot viewer"), "logarithmic"));
    await session.step(141, "When user picks \"Tools > Formula Lines...\" from the context menu of the \"view\" area of scatter plot viewer", () => pickFromAreaContextMenu(page, "Tools > Formula Lines...", "view", el("scatter plot viewer")));
    await session.step(142, "Then \"Formula Lines\" dialog should be visible", () => shouldBe(page, el("\"Formula Lines\" dialog"), "visible"));
    await session.step(143, "When user clicks CANCEL button in \"Formula Lines\" dialog", () => clickOn(page, el("CANCEL button in \"Formula Lines\" dialog")));
    await session.step(144, "And user picks \"X Axis Type > linear\" from the context menu of the \"left edge of x axis\" area of scatter plot viewer", () => pickFromAreaContextMenu(page, "X Axis Type > linear", "left edge of x axis", el("scatter plot viewer")));
    await session.step(145, "Then \"xAxisType\" property of scatter plot viewer should be \"linear\"", () => propertyShouldBe(page, "xAxisType", el("scatter plot viewer"), "linear"));
    await session.step(146, "When user picks \"Tools > Formula Lines...\" from the context menu of the \"view\" area of scatter plot viewer", () => pickFromAreaContextMenu(page, "Tools > Formula Lines...", "view", el("scatter plot viewer")));
    await session.step(147, "Then \"Formula Lines\" dialog should be visible", () => shouldBe(page, el("\"Formula Lines\" dialog"), "visible"));
    await session.step(148, "When user clicks CANCEL button in \"Formula Lines\" dialog", () => clickOn(page, el("CANCEL button in \"Formula Lines\" dialog")));
    await session.step(149, "Then no errors should have been logged", () => noErrors(page));
  });
});
