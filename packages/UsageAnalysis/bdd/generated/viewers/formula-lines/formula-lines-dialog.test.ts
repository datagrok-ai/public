/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/formula-lines/formula-lines-dialog.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot, viewers.line-chart, powerpack.dialogs.formula-lines]
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
import {clickOn, enterInto, shouldBe, shouldHaveText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {setTableTag, setTableTagText} from '@datagrok-libraries/bdd/bindings/platform/data';
import {dialogCloses, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, hasArea, lessInk, loadLayout, moreInk, noErrors, pickFromContextMenu, pickFromOpenMenu, propertyShouldContain, readingIs, readingReads, resizeTo, saveLayout, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Formula Lines dialog and the look it writes", () => {
  const session = feature(test, "features/viewers/formula-lines/formula-lines-dialog.feature", import.meta.url);
  test("A horizontal line added in the dialog lands in the look and leaves with its trash button", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.line-chart", "@realizes:powerpack.dialogs.formula-lines"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "Given user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"]]), [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"]]);
    await session.step(20, "And user resizes scatter plot viewer to 800 by 500", () => resizeTo(page, el("scatter plot viewer"), 800, 500));
    await session.step(21, "Then the \"formula lines\" reading of scatter plot viewer should be 0", () => readingIs(page, "formula lines", el("scatter plot viewer"), 0));
    await session.step(22, "When user picks \"Tools > Formula Lines...\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Tools > Formula Lines...", el("scatter plot viewer")));
    await session.step(23, "Then \"Formula Lines\" dialog should be visible", () => shouldBe(page, el("\"Formula Lines\" dialog"), "visible"));
    await session.step(24, "When user clicks on \"ADD NEW\" button in \"Formula Lines\" dialog", () => clickOn(page, el("\"ADD NEW\" button in \"Formula Lines\" dialog")));
    await session.step(25, "And user picks \"Line - Horizontal\" from the open menu", () => pickFromOpenMenu(page, "Line - Horizontal"));
    await session.step(26, "Then editor of Column input in \"Formula Lines\" dialog should have text \"HEIGHT\"", () => shouldHaveText(page, el("editor of Column input in \"Formula Lines\" dialog"), "HEIGHT"));
    await session.step(27, "And Value input in \"Formula Lines\" dialog should have value \"168.5\"", () => shouldHaveValue(page, el("Value input in \"Formula Lines\" dialog"), "168.5"));
    await session.step(28, "When user enters \"Median height\" into Title input in \"Formula Lines\" dialog", () => enterInto(page, "Median height", el("Title input in \"Formula Lines\" dialog")));
    await session.step(29, "And user clicks OK button in \"Formula Lines\" dialog", () => clickOn(page, el("OK button in \"Formula Lines\" dialog")));
    await session.step(30, "Then the \"Formula Lines\" dialog should close", () => dialogCloses(page, "Formula Lines"));
    await session.step(31, "And \"formulaLines\" property of scatter plot viewer should contain \"${HEIGHT} = 168.5\"", () => propertyShouldContain(page, "formulaLines", el("scatter plot viewer"), "${HEIGHT} = 168.5"));
    await session.step(32, "And \"formulaLines\" property of scatter plot viewer should contain \"\\\"title\\\":\\\"Median height\\\"\"", () => propertyShouldContain(page, "formulaLines", el("scatter plot viewer"), "\"title\":\"Median height\""));
    await session.step(33, "And the \"formula lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "formula lines", el("scatter plot viewer"), 1));
    await session.step(34, "When user picks \"Tools > Formula Lines...\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Tools > Formula Lines...", el("scatter plot viewer")));
    await session.step(35, "And user clicks on Delete button in \"Formula Lines\" dialog", () => clickOn(page, el("Delete button in \"Formula Lines\" dialog")));
    await session.step(36, "And user clicks OK button in \"Formula Lines\" dialog", () => clickOn(page, el("OK button in \"Formula Lines\" dialog")));
    await session.step(37, "Then the \"formula lines\" reading of scatter plot viewer should be 0", () => readingIs(page, "formula lines", el("scatter plot viewer"), 0));
    await session.step(38, "And no errors should have been logged", () => noErrors(page));
  });
  test("Two lines with the same formula and different ranges are both drawn", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.line-chart", "@realizes:powerpack.dialogs.formula-lines"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(41, "Given user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"]]), [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"]]);
    await session.step(44, "And user resizes scatter plot viewer to 800 by 500", () => resizeTo(page, el("scatter plot viewer"), 800, 500));
    await session.step(45, "When user sets \"formulaLines\" property of scatter plot viewer to '[{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${WEIGHT} + 100\",\"min\":60,\"max\":90,\"title\":\"Light\"},{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${WEIGHT} + 100\",\"min\":100,\"max\":150,\"title\":\"Heavy\"}]'", () => setProperty(page, "formulaLines", el("scatter plot viewer"), "[{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${WEIGHT} + 100\",\"min\":60,\"max\":90,\"title\":\"Light\"},{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${WEIGHT} + 100\",\"min\":100,\"max\":150,\"title\":\"Heavy\"}]"));
    await session.step(46, "Then the \"formula lines\" reading of scatter plot viewer should be 2", () => readingIs(page, "formula lines", el("scatter plot viewer"), 2));
    await session.step(47, "And scatter plot viewer should have more ink than before", () => moreInk(page, el("scatter plot viewer")));
    await session.step(48, "When user sets \"formulaLines\" property of scatter plot viewer to \"\"", () => setProperty(page, "formulaLines", el("scatter plot viewer"), ""));
    await session.step(49, "Then the \"formula lines\" reading of scatter plot viewer should be 0", () => readingIs(page, "formula lines", el("scatter plot viewer"), 0));
    await session.step(50, "And no errors should have been logged", () => noErrors(page));
  });
  test("Unchecked items are not drawn and come back when checked again", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.line-chart", "@realizes:powerpack.dialogs.formula-lines"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(53, "Given user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"]]), [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"]]);
    await session.step(56, "And user resizes scatter plot viewer to 800 by 500", () => resizeTo(page, el("scatter plot viewer"), 800, 500));
    await session.step(57, "When user sets \"formulaLines\" property of scatter plot viewer to '[{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${WEIGHT} + 100\",\"title\":\"Light\"},{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${WEIGHT} + 80\",\"title\":\"Heavy\"},{\"type\":\"band\",\"formula\":\"${HEIGHT} in (160.9, 177.6)\",\"orientation\":\"Horizontal\",\"column2\":\"WEIGHT\",\"title\":\"Band\"}]'", () => setProperty(page, "formulaLines", el("scatter plot viewer"), "[{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${WEIGHT} + 100\",\"title\":\"Light\"},{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${WEIGHT} + 80\",\"title\":\"Heavy\"},{\"type\":\"band\",\"formula\":\"${HEIGHT} in (160.9, 177.6)\",\"orientation\":\"Horizontal\",\"column2\":\"WEIGHT\",\"title\":\"Band\"}]"));
    await session.step(58, "Then the \"formula lines\" reading of scatter plot viewer should be 3", () => readingIs(page, "formula lines", el("scatter plot viewer"), 3));
    await session.step(59, "When user sets \"formulaLines\" property of scatter plot viewer to '[{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${WEIGHT} + 100\",\"title\":\"Light\",\"visible\":false},{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${WEIGHT} + 80\",\"title\":\"Heavy\"},{\"type\":\"band\",\"formula\":\"${HEIGHT} in (160.9, 177.6)\",\"orientation\":\"Horizontal\",\"column2\":\"WEIGHT\",\"title\":\"Band\",\"visible\":false}]'", () => setProperty(page, "formulaLines", el("scatter plot viewer"), "[{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${WEIGHT} + 100\",\"title\":\"Light\",\"visible\":false},{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${WEIGHT} + 80\",\"title\":\"Heavy\"},{\"type\":\"band\",\"formula\":\"${HEIGHT} in (160.9, 177.6)\",\"orientation\":\"Horizontal\",\"column2\":\"WEIGHT\",\"title\":\"Band\",\"visible\":false}]"));
    await session.step(60, "Then the \"formula lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "formula lines", el("scatter plot viewer"), 1));
    await session.step(61, "And scatter plot viewer should have less ink than before", () => lessInk(page, el("scatter plot viewer")));
    await session.step(62, "When user sets \"formulaLines\" property of scatter plot viewer to '[{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${WEIGHT} + 100\",\"title\":\"Light\"},{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${WEIGHT} + 80\",\"title\":\"Heavy\"},{\"type\":\"band\",\"formula\":\"${HEIGHT} in (160.9, 177.6)\",\"orientation\":\"Horizontal\",\"column2\":\"WEIGHT\",\"title\":\"Band\"}]'", () => setProperty(page, "formulaLines", el("scatter plot viewer"), "[{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${WEIGHT} + 100\",\"title\":\"Light\"},{\"type\":\"line\",\"formula\":\"${HEIGHT} = ${WEIGHT} + 80\",\"title\":\"Heavy\"},{\"type\":\"band\",\"formula\":\"${HEIGHT} in (160.9, 177.6)\",\"orientation\":\"Horizontal\",\"column2\":\"WEIGHT\",\"title\":\"Band\"}]"));
    await session.step(63, "Then the \"formula lines\" reading of scatter plot viewer should be 3", () => readingIs(page, "formula lines", el("scatter plot viewer"), 3));
    await session.step(64, "And scatter plot viewer should have more ink than before", () => moreInk(page, el("scatter plot viewer")));
    await session.step(65, "When user sets \"formulaLines\" property of scatter plot viewer to \"\"", () => setProperty(page, "formulaLines", el("scatter plot viewer"), ""));
    await session.step(66, "Then no errors should have been logged", () => noErrors(page));
  });
  test("The color and style of a line are kept in the look", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.line-chart", "@realizes:powerpack.dialogs.formula-lines"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(69, "Given user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"]]), [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"]]);
    await session.step(72, "When user sets \"formulaLines\" property of scatter plot viewer to '[{\"type\":\"line\",\"formula\":\"${HEIGHT} = 168.5\",\"color\":\"#ff0000\",\"style\":\"dashed\"}]'", () => setProperty(page, "formulaLines", el("scatter plot viewer"), "[{\"type\":\"line\",\"formula\":\"${HEIGHT} = 168.5\",\"color\":\"#ff0000\",\"style\":\"dashed\"}]"));
    await session.step(73, "Then the \"formula lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "formula lines", el("scatter plot viewer"), 1));
    await session.step(74, "And \"formulaLines\" property of scatter plot viewer should contain \"\\\"color\\\":\\\"#ff0000\\\"\"", () => propertyShouldContain(page, "formulaLines", el("scatter plot viewer"), "\"color\":\"#ff0000\""));
    await session.step(75, "And \"formulaLines\" property of scatter plot viewer should contain \"\\\"style\\\":\\\"dashed\\\"\"", () => propertyShouldContain(page, "formulaLines", el("scatter plot viewer"), "\"style\":\"dashed\""));
    await session.step(76, "When user saves the layout of the current table view", () => saveLayout(page));
    await session.step(77, "And user sets \"formulaLines\" property of scatter plot viewer to \"\"", () => setProperty(page, "formulaLines", el("scatter plot viewer"), ""));
    await session.step(78, "Then the \"formula lines\" reading of scatter plot viewer should be 0", () => readingIs(page, "formula lines", el("scatter plot viewer"), 0));
    await session.step(79, "When user loads the saved layout", () => loadLayout(page));
    await session.step(80, "Then the \"formula lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "formula lines", el("scatter plot viewer"), 1));
    await session.step(81, "And \"formulaLines\" property of scatter plot viewer should contain \"\\\"style\\\":\\\"dashed\\\"\"", () => propertyShouldContain(page, "formulaLines", el("scatter plot viewer"), "\"style\":\"dashed\""));
    await session.step(82, "When user sets \"formulaLines\" property of scatter plot viewer to \"\"", () => setProperty(page, "formulaLines", el("scatter plot viewer"), ""));
    await session.step(83, "Then no errors should have been logged", () => noErrors(page));
  });
  test("A dataframe line is drawn wherever an axis carries its column", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.line-chart", "@realizes:powerpack.dialogs.formula-lines"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(86, "Given user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"]]);
    await session.step(89, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","SEX"],["yColumnNames","WEIGHT"]]), [["xColumnName","SEX"],["yColumnNames","WEIGHT"]]);
    await session.step(92, "Then the \"formula lines\" reading of scatter plot viewer should be 0", () => readingIs(page, "formula lines", el("scatter plot viewer"), 0));
    await session.step(93, "And the \"formula lines\" reading of line chart viewer should be 0", () => readingIs(page, "formula lines", el("line chart viewer"), 0));
    await session.step(94, "When user sets the \".formula-lines\" tag of the table to:", () => setTableTagText(page, ".formula-lines", "[{\"type\":\"line\",\"formula\":\"${WEIGHT} = 100\",\"title\":\"Reference weight\"}]"));
    await session.step(98, "Then the \"formula lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "formula lines", el("scatter plot viewer"), 1));
    await session.step(99, "And the \"formula lines\" reading of line chart viewer should be 0", () => readingIs(page, "formula lines", el("line chart viewer"), 0));
    await session.step(100, "When user sets \"xColumnName\" property of line chart viewer to \"USUBJID\"", () => setProperty(page, "xColumnName", el("line chart viewer"), "USUBJID"));
    await session.step(101, "Then the \"aggregated\" reading of line chart viewer should be \"false\"", () => readingReads(page, "aggregated", el("line chart viewer"), "false"));
    await session.step(102, "And the \"formula lines\" reading of line chart viewer should be 1", () => readingIs(page, "formula lines", el("line chart viewer"), 1));
    await session.step(103, "And line chart viewer should have a \"formula line Reference weight\" area", () => hasArea(page, el("line chart viewer"), "formula line Reference weight"));
    await session.step(104, "When user sets \"xColumnName\" property of line chart viewer to \"SEX\"", () => setProperty(page, "xColumnName", el("line chart viewer"), "SEX"));
    await session.step(105, "And user sets the \".formula-lines\" tag of the table to \"[]\"", () => setTableTag(page, ".formula-lines", "[]"));
    await session.step(106, "Then the \"formula lines\" reading of scatter plot viewer should be 0", () => readingIs(page, "formula lines", el("scatter plot viewer"), 0));
    await session.step(107, "And the \"formula lines\" reading of line chart viewer should be 0", () => readingIs(page, "formula lines", el("line chart viewer"), 0));
    await session.step(108, "And no errors should have been logged", () => noErrors(page));
  });
});
