/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/legend/legend-structures.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.legend]
--- */
import {test} from '@playwright/test';
import '../../../bindings/connections.js';
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
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {closeAllViews, openDataset, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, legendItemsAsStructures, legendItemsAsText, legendLists, loadLayout, noErrors, propertyShouldBe, saveLayoutToServer, setProperties, setProperty, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {setInnerProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Molecules in the legend", () => {
  const session = feature(test, "features/viewers/legend/legend-structures.feature", import.meta.url);
  test("Molecules in the legend", {tag: ["@journey", "@viewers", "@realizes:viewers.legend"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 10, page);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(25, "Then \"Core\" column should have semantic type \"Molecule\"", () => columnSemType(page, "Core", "Molecule"));
    await session.step(26, "When user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["colorColumnName","Core"],["Legend Visibility","Always"],["Legend Position","Right"]]), [["colorColumnName","Core"],["Legend Visibility","Always"],["Legend Position","Right"]]);
    await session.step(30, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["splitColumnName","Core"],["Legend Visibility","Always"],["Legend Position","Right"]]), [["splitColumnName","Core"],["Legend Visibility","Always"],["Legend Position","Right"]]);
    await session.step(34, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["splitColumnNames","Core"],["Legend Visibility","Always"],["Legend Position","Right"]]), [["splitColumnNames","Core"],["Legend Visibility","Always"],["Legend Position","Right"]]);
    await session.step(38, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["splitColumnName","Series"],["stackColumnName","Core"],["Legend Visibility","Always"],["Legend Position","Right"]]), [["splitColumnName","Series"],["stackColumnName","Core"],["Legend Visibility","Always"],["Legend Position","Right"]]);
    await session.step(43, "And user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["categoryColumnName","Core"],["Legend Visibility","Always"],["Legend Position","Right"]]), [["categoryColumnName","Core"],["Legend Visibility","Always"],["Legend Position","Right"]]);
    await session.step(47, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["xColumnNames","Core"],["Viewer Type","Scatter plot"],["Legend Visibility","Always"],["Legend Position","Right"]]), [["xColumnNames","Core"],["Viewer Type","Scatter plot"],["Legend Visibility","Always"],["Legend Position","Right"]]);
    await session.step(52, "And user sets \"colorColumnName\" inner property of trellis plot viewer to \"Core\"", () => setInnerProperty(page, "colorColumnName", el("trellis plot viewer"), "Core"));
    await session.step(53, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["categoryColumnNames","Core"],["markerColorColumnName","Core"],["Legend Visibility","Always"],["Legend Position","Right"]]), [["categoryColumnNames","Core"],["markerColorColumnName","Core"],["Legend Visibility","Always"],["Legend Position","Right"]]);
    await run.scenario("The scatter plot legend draws the seven cores as structures [viewer=scatter plot]", async () => {
      await session.step(60, "Then the legend of scatter plot viewer should list 7 items", () => legendLists(page, el("scatter plot viewer"), 7));
      await session.step(61, "And every item in the legend of scatter plot viewer should be drawn as a structure", () => legendItemsAsStructures(page, el("scatter plot viewer")));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The histogram legend draws the seven cores as structures [viewer=histogram]", async () => {
      await session.step(60, "Then the legend of histogram viewer should list 7 items", () => legendLists(page, el("histogram viewer"), 7));
      await session.step(61, "And every item in the legend of histogram viewer should be drawn as a structure", () => legendItemsAsStructures(page, el("histogram viewer")));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The line chart legend draws the seven cores as structures [viewer=line chart]", async () => {
      await session.step(60, "Then the legend of line chart viewer should list 7 items", () => legendLists(page, el("line chart viewer"), 7));
      await session.step(61, "And every item in the legend of line chart viewer should be drawn as a structure", () => legendItemsAsStructures(page, el("line chart viewer")));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The bar chart legend draws the seven cores as structures [viewer=bar chart]", async () => {
      await session.step(60, "Then the legend of bar chart viewer should list 7 items", () => legendLists(page, el("bar chart viewer"), 7));
      await session.step(61, "And every item in the legend of bar chart viewer should be drawn as a structure", () => legendItemsAsStructures(page, el("bar chart viewer")));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The pie chart legend draws the seven cores as structures [viewer=pie chart]", async () => {
      await session.step(60, "Then the legend of pie chart viewer should list 7 items", () => legendLists(page, el("pie chart viewer"), 7));
      await session.step(61, "And every item in the legend of pie chart viewer should be drawn as a structure", () => legendItemsAsStructures(page, el("pie chart viewer")));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The trellis plot legend draws the seven cores as structures [viewer=trellis plot]", async () => {
      await session.step(60, "Then the legend of trellis plot viewer should list 7 items", () => legendLists(page, el("trellis plot viewer"), 7));
      await session.step(61, "And every item in the legend of trellis plot viewer should be drawn as a structure", () => legendItemsAsStructures(page, el("trellis plot viewer")));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The box plot legend draws the seven cores as structures [viewer=box plot]", async () => {
      await session.step(60, "Then the legend of box plot viewer should list 7 items", () => legendLists(page, el("box plot viewer"), 7));
      await session.step(61, "And every item in the legend of box plot viewer should be drawn as a structure", () => legendItemsAsStructures(page, el("box plot viewer")));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Markers on Core add structure items next to the color items of Series", async () => {
      await session.step(75, "When user sets \"Color\" property of scatter plot viewer to \"Series\"", () => setProperty(page, "Color", el("scatter plot viewer"), "Series"));
      await session.step(76, "Then the legend of scatter plot viewer should list 5 items", () => legendLists(page, el("scatter plot viewer"), 5));
      await session.step(77, "And every item in the legend of scatter plot viewer should be drawn as text", () => legendItemsAsText(page, el("scatter plot viewer")));
      await session.step(78, "When user sets \"Markers\" property of scatter plot viewer to \"Core\"", () => setProperty(page, "Markers", el("scatter plot viewer"), "Core"));
      await session.step(79, "Then the legend of scatter plot viewer should list 12 items", () => legendLists(page, el("scatter plot viewer"), 12));
      await session.step(80, "And thumbnail of last legend item in legend of scatter plot viewer should be visible", () => shouldBe(page, el("thumbnail of last legend item in legend of scatter plot viewer"), "visible"));
      await session.step(81, "And thumbnail of first legend item in legend of scatter plot viewer should be absent", () => shouldBe(page, el("thumbnail of first legend item in legend of scatter plot viewer"), "absent"));
      await session.step(82, "When user sets \"Color\" property of scatter plot viewer to \"Id\"", () => setProperty(page, "Color", el("scatter plot viewer"), "Id"));
      await session.step(83, "Then the legend of scatter plot viewer should list 107 items", () => legendLists(page, el("scatter plot viewer"), 107));
      await session.step(84, "And thumbnail of last legend item in legend of scatter plot viewer should be visible", () => shouldBe(page, el("thumbnail of last legend item in legend of scatter plot viewer"), "visible"));
      await session.step(85, "And thumbnail of first legend item in legend of scatter plot viewer should be absent", () => shouldBe(page, el("thumbnail of first legend item in legend of scatter plot viewer"), "absent"));
      await session.step(86, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Color","Core"],["Markers",""]]), [["Color","Core"],["Markers",""]]);
      await session.step(89, "Then the legend of scatter plot viewer should list 7 items", () => legendLists(page, el("scatter plot viewer"), 7));
      await session.step(90, "And every item in the legend of scatter plot viewer should be drawn as a structure", () => legendItemsAsStructures(page, el("scatter plot viewer")));
      await session.step(91, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The structures come back from a saved layout", async () => {
      await session.step(94, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(95, "And user sets \"Color\" property of scatter plot viewer to \"Series\"", () => setProperty(page, "Color", el("scatter plot viewer"), "Series"));
      await session.step(96, "And user sets \"splitColumnName\" property of histogram viewer to \"Series\"", () => setProperty(page, "splitColumnName", el("histogram viewer"), "Series"));
      await session.step(97, "And user sets \"splitColumnNames\" property of line chart viewer to \"Series\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), "Series"));
      await session.step(98, "And user sets \"stackColumnName\" property of bar chart viewer to \"Series\"", () => setProperty(page, "stackColumnName", el("bar chart viewer"), "Series"));
      await session.step(99, "And user sets \"Category\" property of pie chart viewer to \"Series\"", () => setProperty(page, "Category", el("pie chart viewer"), "Series"));
      await session.step(100, "And user sets \"colorColumnName\" inner property of trellis plot viewer to \"Series\"", () => setInnerProperty(page, "colorColumnName", el("trellis plot viewer"), "Series"));
      await session.step(101, "And user sets \"markerColorColumnName\" property of box plot viewer to \"Series\"", () => setProperty(page, "markerColorColumnName", el("box plot viewer"), "Series"));
      await session.step(102, "Then every item in the legend of scatter plot viewer should be drawn as text", () => legendItemsAsText(page, el("scatter plot viewer")));
      await session.step(103, "And every item in the legend of histogram viewer should be drawn as text", () => legendItemsAsText(page, el("histogram viewer")));
      await session.step(104, "And every item in the legend of line chart viewer should be drawn as text", () => legendItemsAsText(page, el("line chart viewer")));
      await session.step(105, "And every item in the legend of bar chart viewer should be drawn as text", () => legendItemsAsText(page, el("bar chart viewer")));
      await session.step(106, "And every item in the legend of pie chart viewer should be drawn as text", () => legendItemsAsText(page, el("pie chart viewer")));
      await session.step(107, "And every item in the legend of trellis plot viewer should be drawn as text", () => legendItemsAsText(page, el("trellis plot viewer")));
      await session.step(108, "And every item in the legend of box plot viewer should be drawn as text", () => legendItemsAsText(page, el("box plot viewer")));
      await session.step(109, "When user loads the saved layout", () => loadLayout(page));
      await session.step(110, "Then \"Category\" property of pie chart viewer should be \"Core\"", () => propertyShouldBe(page, "Category", el("pie chart viewer"), "Core"));
      await session.step(111, "And the legend of pie chart viewer should list 7 items", () => legendLists(page, el("pie chart viewer"), 7));
      await session.step(112, "And every item in the legend of pie chart viewer should be drawn as a structure", () => legendItemsAsStructures(page, el("pie chart viewer")));
      await session.step(113, "And every item in the legend of scatter plot viewer should be drawn as a structure", () => legendItemsAsStructures(page, el("scatter plot viewer")));
      await session.step(114, "And every item in the legend of histogram viewer should be drawn as a structure", () => legendItemsAsStructures(page, el("histogram viewer")));
      await session.step(115, "And every item in the legend of line chart viewer should be drawn as a structure", () => legendItemsAsStructures(page, el("line chart viewer")));
      await session.step(116, "And every item in the legend of bar chart viewer should be drawn as a structure", () => legendItemsAsStructures(page, el("bar chart viewer")));
      await session.step(117, "And every item in the legend of trellis plot viewer should be drawn as a structure", () => legendItemsAsStructures(page, el("trellis plot viewer")));
      await session.step(118, "And every item in the legend of box plot viewer should be drawn as a structure", () => legendItemsAsStructures(page, el("box plot viewer")));
      await session.step(119, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The structures come back from a saved project", async () => {
      await session.step(122, "When user saves the current view as project \"bdd-legend-structures\"", () => saveAsProject(page, "bdd-legend-structures"));
      await session.step(123, "And user closes all views", () => closeAllViews(page));
      await session.step(124, "And user opens the \"bdd-legend-structures\" project", () => openProject(page, "bdd-legend-structures"));
      await session.step(125, "Then \"Core\" column should have semantic type \"Molecule\"", () => columnSemType(page, "Core", "Molecule"));
      await session.step(126, "And the open tableview should have 1 box plot viewer", () => viewerCount(page, 1, "box plot"));
      await session.step(127, "And the legend of scatter plot viewer should list 7 items", () => legendLists(page, el("scatter plot viewer"), 7));
      await session.step(128, "And every item in the legend of scatter plot viewer should be drawn as a structure", () => legendItemsAsStructures(page, el("scatter plot viewer")));
      await session.step(129, "And every item in the legend of histogram viewer should be drawn as a structure", () => legendItemsAsStructures(page, el("histogram viewer")));
      await session.step(130, "And every item in the legend of line chart viewer should be drawn as a structure", () => legendItemsAsStructures(page, el("line chart viewer")));
      await session.step(131, "And every item in the legend of bar chart viewer should be drawn as a structure", () => legendItemsAsStructures(page, el("bar chart viewer")));
      await session.step(132, "And every item in the legend of pie chart viewer should be drawn as a structure", () => legendItemsAsStructures(page, el("pie chart viewer")));
      await session.step(133, "And every item in the legend of trellis plot viewer should be drawn as a structure", () => legendItemsAsStructures(page, el("trellis plot viewer")));
      await session.step(134, "And every item in the legend of box plot viewer should be drawn as a structure", () => legendItemsAsStructures(page, el("box plot viewer")));
      await session.step(135, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
