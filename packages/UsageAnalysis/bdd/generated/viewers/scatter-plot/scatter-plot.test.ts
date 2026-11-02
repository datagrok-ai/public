/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/scatter-plot/scatter-plot.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot, public#3041, public#1156, public#344, core#1141, public#2285]
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
import {clickOn, shouldBe, shouldHaveText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {distinctValues} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {addCalculated, filterPassesAll, removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaPainted, boundTable, closeContextMenu, hasArea, hasNoArea, lessInk, moreInk, noErrors, notRepainted, openContextMenu, paintedInColors, propertyShouldNotBe, readingAtLeast, readingIs, repainted, repaintedBy, setProperties, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Scatter plot property surface", () => {
  const session = feature(test, "features/viewers/scatter-plot/scatter-plot.feature", import.meta.url);
  test("Scatter plot property surface", {tag: ["@journey", "@viewers", "@realizes:viewers.scatter-plot", "@realizes:public#3041", "@realizes:public#1156", "@realizes:public#344", "@realizes:core#1141", "@realizes:public#2285"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"]]), [["X","WEIGHT"],["Y","HEIGHT"]]);
    await session.step(23, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
    await run.scenario("The axis histograms appear beside the plot and re-bin", async () => {
      await session.step(26, "Then scatter plot viewer should not have an \"x histogram\" area", () => hasNoArea(page, el("scatter plot viewer"), "x histogram"));
      await session.step(27, "And scatter plot viewer should not have a \"y histogram\" area", () => hasNoArea(page, el("scatter plot viewer"), "y histogram"));
      await session.step(28, "When user sets \"Show X Histogram\" property of scatter plot viewer to \"true\"", () => setProperty(page, "Show X Histogram", el("scatter plot viewer"), "true"));
      await session.step(29, "Then scatter plot viewer should have an \"x histogram\" area", () => hasArea(page, el("scatter plot viewer"), "x histogram"));
      await session.step(30, "And the \"x histogram\" area of scatter plot viewer should be painted", () => areaPainted(page, "x histogram", el("scatter plot viewer")));
      await session.step(31, "And scatter plot viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("scatter plot viewer"), 500));
      await session.step(32, "When user sets \"Show Y Histogram\" property of scatter plot viewer to \"true\"", () => setProperty(page, "Show Y Histogram", el("scatter plot viewer"), "true"));
      await session.step(33, "Then scatter plot viewer should have a \"y histogram\" area", () => hasArea(page, el("scatter plot viewer"), "y histogram"));
      await session.step(34, "And the \"y histogram\" area of scatter plot viewer should be painted", () => areaPainted(page, "y histogram", el("scatter plot viewer")));
      await session.step(35, "When user sets \"Histogram Bins\" property of scatter plot viewer to \"20\"", () => setProperty(page, "Histogram Bins", el("scatter plot viewer"), "20"));
      await session.step(36, "Then scatter plot viewer should have repainted", () => repainted(page, el("scatter plot viewer")));
      await session.step(37, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Show X Histogram","false"],["Show Y Histogram","false"],["Histogram Bins","10"]]), [["Show X Histogram","false"],["Show Y Histogram","false"],["Histogram Bins","10"]]);
      await session.step(41, "Then scatter plot viewer should not have an \"x histogram\" area", () => hasNoArea(page, el("scatter plot viewer"), "x histogram"));
      await session.step(42, "And scatter plot viewer should not have a \"y histogram\" area", () => hasNoArea(page, el("scatter plot viewer"), "y histogram"));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Grid lines, the axes and the on-viewer selectors can be turned off", async () => {
      await session.step(46, "Then X column input in scatter plot viewer should be visible", () => shouldBe(page, el("X column input in scatter plot viewer"), "visible"));
      await session.step(47, "And Y column input in scatter plot viewer should be visible", () => shouldBe(page, el("Y column input in scatter plot viewer"), "visible"));
      await session.step(48, "And scatter plot viewer should have an \"x axis\" area", () => hasArea(page, el("scatter plot viewer"), "x axis"));
      await session.step(49, "And scatter plot viewer should have a \"y axis\" area", () => hasArea(page, el("scatter plot viewer"), "y axis"));
      await session.step(50, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Show Vertical Grid Lines","false"],["Show Horizontal Grid Lines","false"]]), [["Show Vertical Grid Lines","false"],["Show Horizontal Grid Lines","false"]]);
      await session.step(53, "Then scatter plot viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("scatter plot viewer"), 500));
      await session.step(54, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Show X Axis","false"],["Show Y Axis","false"]]), [["Show X Axis","false"],["Show Y Axis","false"]]);
      await session.step(57, "Then scatter plot viewer should not have an \"x axis\" area", () => hasNoArea(page, el("scatter plot viewer"), "x axis"));
      await session.step(58, "And scatter plot viewer should not have a \"y axis\" area", () => hasNoArea(page, el("scatter plot viewer"), "y axis"));
      await session.step(59, "And scatter plot viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("scatter plot viewer"), 500));
      await session.step(60, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Show X Selector","false"],["Show Y Selector","false"]]), [["Show X Selector","false"],["Show Y Selector","false"]]);
      await session.step(63, "Then X column input in scatter plot viewer should be hidden", () => shouldBe(page, el("X column input in scatter plot viewer"), "hidden"));
      await session.step(64, "And X column input in scatter plot viewer should be present", () => shouldBe(page, el("X column input in scatter plot viewer"), "present"));
      await session.step(65, "And Y column input in scatter plot viewer should be hidden", () => shouldBe(page, el("Y column input in scatter plot viewer"), "hidden"));
      await session.step(66, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Show X Axis","true"],["Show Y Axis","true"],["Show X Selector","true"],["Show Y Selector","true"],["Show Vertical Grid Lines","true"],["Show Horizontal Grid Lines","true"]]), [["Show X Axis","true"],["Show Y Axis","true"],["Show X Selector","true"],["Show Y Selector","true"],["Show Vertical Grid Lines","true"],["Show Horizontal Grid Lines","true"]]);
      await session.step(73, "Then scatter plot viewer should have an \"x axis\" area", () => hasArea(page, el("scatter plot viewer"), "x axis"));
      await session.step(74, "And scatter plot viewer should have a \"y axis\" area", () => hasArea(page, el("scatter plot viewer"), "y axis"));
      await session.step(75, "And X column input in scatter plot viewer should be visible", () => shouldBe(page, el("X column input in scatter plot viewer"), "visible"));
      await session.step(76, "And Y column input in scatter plot viewer should be visible", () => shouldBe(page, el("Y column input in scatter plot viewer"), "visible"));
      await session.step(77, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Whisker columns draw error bars around the markers", async () => {
      await session.step(80, "Then scatter plot viewer should not have a \"whiskers\" area", () => hasNoArea(page, el("scatter plot viewer"), "whiskers"));
      await session.step(81, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["X Whisker Min Column","AGE"],["X Whisker Max Column","WEIGHT"],["Y Whisker Min Column","HEIGHT"],["Y Whisker Max Column","WEIGHT"]]), [["X Whisker Min Column","AGE"],["X Whisker Max Column","WEIGHT"],["Y Whisker Min Column","HEIGHT"],["Y Whisker Max Column","WEIGHT"]]);
      await session.step(86, "Then scatter plot viewer should have a \"whiskers\" area", () => hasArea(page, el("scatter plot viewer"), "whiskers"));
      await session.step(87, "And scatter plot viewer should have more ink than before", () => moreInk(page, el("scatter plot viewer")));
      await session.step(88, "And no errors should have been logged", () => noErrors(page));
      await session.step(89, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["X Whisker Min Column",""],["X Whisker Max Column",""],["Y Whisker Min Column",""],["Y Whisker Max Column",""]]), [["X Whisker Min Column",""],["X Whisker Max Column",""],["Y Whisker Min Column",""],["Y Whisker Max Column",""]]);
      await session.step(94, "Then scatter plot viewer should not have a \"whiskers\" area", () => hasNoArea(page, el("scatter plot viewer"), "whiskers"));
      await session.step(95, "And scatter plot viewer should have less ink than before", () => lessInk(page, el("scatter plot viewer")));
      await session.step(96, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The context menu carries Reset View, the Lasso Tool and the tool groups", async () => {
      await session.step(99, "When user opens the context menu of scatter plot viewer", () => openContextMenu(page, el("scatter plot viewer")));
      await session.step(100, "Then \"Reset View\" menu item in context menu should be visible", () => shouldBe(page, el("\"Reset View\" menu item in context menu"), "visible"));
      await session.step(101, "And \"Lasso Tool\" menu item in context menu should be visible", () => shouldBe(page, el("\"Lasso Tool\" menu item in context menu"), "visible"));
      await session.step(102, "And \"Tools\" menu item in context menu should be visible", () => shouldBe(page, el("\"Tools\" menu item in context menu"), "visible"));
      await session.step(103, "And \"Markers\" menu item in context menu should be visible", () => shouldBe(page, el("\"Markers\" menu item in context menu"), "visible"));
      await session.step(104, "And \"Labels\" menu item in context menu should be visible", () => shouldBe(page, el("\"Labels\" menu item in context menu"), "visible"));
      await session.step(105, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(106, "Then context menu should be hidden", () => shouldBe(page, el("context menu"), "hidden"));
      await session.step(107, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(108, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Title and description show on the viewer and go away again", async () => {
      await session.step(111, "Then title of scatter plot viewer should not contain the text \"Test Plot\"", () => shouldNotContainText(page, el("title of scatter plot viewer"), "Test Plot"));
      await session.step(112, "When user sets \"Title\" property of scatter plot viewer to \"Test Plot\"", () => setProperty(page, "Title", el("scatter plot viewer"), "Test Plot"));
      await session.step(113, "Then title of scatter plot viewer should have text \"Test Plot\"", () => shouldHaveText(page, el("title of scatter plot viewer"), "Test Plot"));
      await session.step(114, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Description","Test description"],["Description Visibility Mode","Always"]]), [["Description","Test description"],["Description Visibility Mode","Always"]]);
      await session.step(117, "Then description of scatter plot viewer should have text \"Test description\"", () => shouldHaveText(page, el("description of scatter plot viewer"), "Test description"));
      await session.step(118, "And title of scatter plot viewer should have text \"Test Plot\"", () => shouldHaveText(page, el("title of scatter plot viewer"), "Test Plot"));
      await session.step(119, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Title",""],["Description",""]]), [["Title",""],["Description",""]]);
      await session.step(122, "Then title of scatter plot viewer should not contain the text \"Test Plot\"", () => shouldNotContainText(page, el("title of scatter plot viewer"), "Test Plot"));
      await session.step(123, "And description of scatter plot viewer should be absent", () => shouldBe(page, el("description of scatter plot viewer"), "absent"));
      await session.step(124, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Lines Order draws connecting lines and Lines By defaults to the color column", async () => {
      await session.step(127, "When user clicks on settings icon of scatter plot viewer", () => clickOn(page, el("settings icon of scatter plot viewer")));
      await session.step(128, "Then \"Lines By\" property should be disabled", () => shouldBe(page, el("\"Lines By\" property"), "disabled"));
      await session.step(129, "And scatter plot viewer should not have a \"lines\" area", () => hasNoArea(page, el("scatter plot viewer"), "lines"));
      await session.step(130, "And \"RACE\" column should have at least 3 distinct values", () => distinctValues(page, "RACE", 3));
      await session.step(131, "When user sets \"Color\" property of scatter plot viewer to \"RACE\"", () => setProperty(page, "Color", el("scatter plot viewer"), "RACE"));
      await session.step(132, "Then scatter plot viewer should be painted in at least 3 colors", () => paintedInColors(page, el("scatter plot viewer"), 3));
      await session.step(133, "When user sets \"Lines Order\" property of scatter plot viewer to \"AGE\"", () => setProperty(page, "Lines Order", el("scatter plot viewer"), "AGE"));
      await session.step(134, "Then scatter plot viewer should have a \"lines\" area", () => hasArea(page, el("scatter plot viewer"), "lines"));
      await session.step(135, "And \"Lines By\" property should be enabled", () => shouldBe(page, el("\"Lines By\" property"), "enabled"));
      await session.step(136, "And scatter plot viewer should have more ink than before", () => moreInk(page, el("scatter plot viewer")));
      await session.step(137, "When user sets \"Lines By\" property of scatter plot viewer to \"RACE\"", () => setProperty(page, "Lines By", el("scatter plot viewer"), "RACE"));
      await session.step(138, "Then scatter plot viewer should not have repainted", () => notRepainted(page, el("scatter plot viewer")));
      await session.step(139, "When user sets \"Lines By\" property of scatter plot viewer to \"SEX\"", () => setProperty(page, "Lines By", el("scatter plot viewer"), "SEX"));
      await session.step(140, "Then scatter plot viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("scatter plot viewer"), 500));
      await session.step(141, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Lines By",""],["Lines Order",""],["Color",""]]), [["Lines By",""],["Lines Order",""],["Color",""]]);
      await session.step(145, "Then scatter plot viewer should not have a \"lines\" area", () => hasNoArea(page, el("scatter plot viewer"), "lines"));
      await session.step(146, "And \"Lines By\" property should be disabled", () => shouldBe(page, el("\"Lines By\" property"), "disabled"));
      await session.step(147, "And scatter plot viewer should have less ink than before", () => lessInk(page, el("scatter plot viewer")));
      await session.step(148, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Data > Table rebinds the plot to another table and back, and the Filter formula narrows it, survives a new column and an invalid one is ignored", async () => {
      await session.step(152, "Given user opens spgi dataset", () => openDataset(page, ds("spgi")));
      await session.step(153, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
      await session.step(154, "When user sets \"table\" property of scatter plot viewer to \"spgi-100\"", () => setProperty(page, "table", el("scatter plot viewer"), "spgi-100"));
      await session.step(155, "Then scatter plot viewer should be bound to table \"spgi-100\"", () => boundTable(page, el("scatter plot viewer"), "spgi-100"));
      await session.step(156, "And \"X\" property of scatter plot viewer should not be \"WEIGHT\"", () => propertyShouldNotBe(page, "X", el("scatter plot viewer"), "WEIGHT"));
      await session.step(157, "And \"X\" property of scatter plot viewer should not be \"\"", () => propertyShouldNotBe(page, "X", el("scatter plot viewer"), ""));
      await session.step(158, "And \"Y\" property of scatter plot viewer should not be \"HEIGHT\"", () => propertyShouldNotBe(page, "Y", el("scatter plot viewer"), "HEIGHT"));
      await session.step(159, "And the \"rows shown\" reading of scatter plot viewer should be at least 1", () => readingAtLeast(page, "rows shown", el("scatter plot viewer"), 1));
      await session.step(160, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["X","Stereo Category"],["Y","Chemist 521"],["Label Columns","Core"]]), [["X","Stereo Category"],["Y","Chemist 521"],["Label Columns","Core"]]);
      await session.step(164, "Then the \"labels shown\" reading of scatter plot viewer should be 43", () => readingIs(page, "labels shown", el("scatter plot viewer"), 43));
      await session.step(165, "When user sets \"Label Columns\" property of scatter plot viewer to \"\"", () => setProperty(page, "Label Columns", el("scatter plot viewer"), ""));
      await session.step(166, "When user sets \"table\" property of scatter plot viewer to \"demog-1000\"", () => setProperty(page, "table", el("scatter plot viewer"), "demog-1000"));
      await session.step(167, "Then scatter plot viewer should be bound to table \"demog-1000\"", () => boundTable(page, el("scatter plot viewer"), "demog-1000"));
      await session.step(168, "And \"X\" property of scatter plot viewer should not be \"\"", () => propertyShouldNotBe(page, "X", el("scatter plot viewer"), ""));
      await session.step(169, "And the \"rows shown\" reading of scatter plot viewer should be at least 1", () => readingAtLeast(page, "rows shown", el("scatter plot viewer"), 1));
      await session.step(170, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["X","WEIGHT"],["Y","HEIGHT"]]), [["X","WEIGHT"],["Y","HEIGHT"]]);
      await session.step(173, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(174, "When user sets \"Filter\" property of scatter plot viewer to \"${AGE} > 44\"", () => setProperty(page, "Filter", el("scatter plot viewer"), "${AGE} > 44"));
      await session.step(175, "Then scatter plot viewer should show 466 rows", () => showsRows(page, el("scatter plot viewer"), 466));
      await session.step(176, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(177, "When user adds a calculated column \"AGE2\" with formula \"${AGE} * 2\"", () => addCalculated(page, "AGE2", "${AGE} * 2"));
      await session.step(178, "And user sets \"Color\" property of scatter plot viewer to \"AGE2\"", () => setProperty(page, "Color", el("scatter plot viewer"), "AGE2"));
      await session.step(179, "Then scatter plot viewer should show 466 rows", () => showsRows(page, el("scatter plot viewer"), 466));
      await session.step(180, "When user sets \"Color\" property of scatter plot viewer to \"\"", () => setProperty(page, "Color", el("scatter plot viewer"), ""));
      await session.step(181, "And user removes \"AGE2\" column", () => removeColumn(page, "AGE2"));
      await session.step(182, "And user sets \"Filter\" property of scatter plot viewer to '${RACE} = \"Black\"'", () => setProperty(page, "Filter", el("scatter plot viewer"), "${RACE} = \"Black\""));
      await session.step(183, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(184, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(185, "And no errors should have been logged", () => noErrors(page));
      await session.step(186, "When user sets \"Filter\" property of scatter plot viewer to \"\"", () => setProperty(page, "Filter", el("scatter plot viewer"), ""));
      await session.step(187, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(188, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Marker Default Size changes the markers (public#2285)", async () => {
      await session.step(192, "When user sets \"markerDefaultSize\" property of scatter plot viewer to \"12\"", () => setProperty(page, "markerDefaultSize", el("scatter plot viewer"), "12"));
      await session.step(193, "Then scatter plot viewer should have more ink than before", () => moreInk(page, el("scatter plot viewer")));
      await session.step(194, "When user sets \"markerDefaultSize\" property of scatter plot viewer to \"2\"", () => setProperty(page, "markerDefaultSize", el("scatter plot viewer"), "2"));
      await session.step(195, "Then scatter plot viewer should have less ink than before", () => lessInk(page, el("scatter plot viewer")));
      await session.step(196, "When user sets \"markerDefaultSize\" property of scatter plot viewer to \"\"", () => setProperty(page, "markerDefaultSize", el("scatter plot viewer"), ""));
      await session.step(197, "Then scatter plot viewer should have more ink than before", () => moreInk(page, el("scatter plot viewer")));
      await session.step(198, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
