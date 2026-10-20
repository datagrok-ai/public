/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/molecule-categories.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/file-formats.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {openDataset, packageInstalled} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, hasArea, noErrors, readingAtLeast, readingHigher, readingIs, readingLower, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {dragRangeHandle} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Molecules as the categories of five viewers", () => {
  const session = feature(test, "features/viewers/molecule-categories.feature", import.meta.url);
  test("Molecules as the categories of five viewers", {tag: ["@journey", "@viewers"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(20, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(21, "Then \"Core\" column should have semantic type \"Molecule\"", () => columnSemType(page, "Core", "Molecule"));
    await session.step(22, "And \"Structure\" column should have semantic type \"Molecule\"", () => columnSemType(page, "Structure", "Molecule"));
    await run.scenario("The bar chart splits into a bar per core", async () => {
      await session.step(25, "When user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["splitColumnName","Core"]]), [["splitColumnName","Core"]]);
      await session.step(27, "Then the \"bars\" reading of bar chart viewer should be 7", () => readingIs(page, "bars", el("bar chart viewer"), 7));
      await session.step(28, "And bar chart viewer should have a \"bar n1nc([*:3])n([*:2])c1[*:1]\" area", () => hasArea(page, el("bar chart viewer"), "bar n1nc([*:3])n([*:2])c1[*:1]"));
      await session.step(29, "And bar chart viewer should have a \"y axis\" area", () => hasArea(page, el("bar chart viewer"), "y axis"));
      await session.step(30, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Dragging the bar chart's category slider over the molecules shows fewer bars", async () => {
      await session.step(33, "When user sets \"splitColumnName\" property of bar chart viewer to \"Structure\"", () => setProperty(page, "splitColumnName", el("bar chart viewer"), "Structure"));
      await session.step(34, "Then the \"bars\" reading of bar chart viewer should be 19", () => readingIs(page, "bars", el("bar chart viewer"), 19));
      await session.step(35, "When user drags the max handle of the \"y\" range slider of bar chart viewer by -120 pixels", () => dragRangeHandle(page, "max", "y", el("bar chart viewer"), -120));
      await session.step(36, "Then the \"bars\" reading of bar chart viewer should be lower than before", () => readingLower(page, "bars", el("bar chart viewer")));
      await session.step(37, "And the \"bars\" reading of bar chart viewer should be at least 1", () => readingAtLeast(page, "bars", el("bar chart viewer"), 1));
      await session.step(38, "And bar chart viewer should have a \"y axis\" area", () => hasArea(page, el("bar chart viewer"), "y axis"));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
      await session.step(40, "When user clicks on close icon of bar chart viewer", () => clickOn(page, el("close icon of bar chart viewer")));
      await session.step(41, "Then bar chart viewer should be absent", () => shouldBe(page, el("bar chart viewer"), "absent"));
    });
    await run.scenario("The box plot lays out a column per core", async () => {
      await session.step(44, "When user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["categoryColumnNames","Core"],["valueColumnName","Chemical Space X"]]), [["categoryColumnNames","Core"],["valueColumnName","Chemical Space X"]]);
      await session.step(47, "Then box plot viewer should have a \"category n1nc([*:3])n([*:2])c1[*:1]\" area", () => hasArea(page, el("box plot viewer"), "category n1nc([*:3])n([*:2])c1[*:1]"));
      await session.step(48, "And box plot viewer should have a \"category CC1(N[*:1])CCN([*:2])C1\" area", () => hasArea(page, el("box plot viewer"), "category CC1(N[*:1])CCN([*:2])C1"));
      await session.step(49, "And no errors should have been logged", () => noErrors(page));
      await session.step(50, "When user clicks on close icon of box plot viewer", () => clickOn(page, el("close icon of box plot viewer")));
      await session.step(51, "Then box plot viewer should be absent", () => shouldBe(page, el("box plot viewer"), "absent"));
    });
    await run.scenario("The pie chart cuts a slice per core", async () => {
      await session.step(54, "When user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["categoryColumnName","Core"]]), [["categoryColumnName","Core"]]);
      await session.step(56, "Then the \"slices\" reading of pie chart viewer should be 7", () => readingIs(page, "slices", el("pie chart viewer"), 7));
      await session.step(57, "And pie chart viewer should have a 'slice \"n1nc([*:3])n([*:2])c1[*:1]\"' area", () => hasArea(page, el("pie chart viewer"), "slice \"n1nc([*:3])n([*:2])c1[*:1]\""));
      await session.step(58, "And no errors should have been logged", () => noErrors(page));
      await session.step(59, "When user clicks on close icon of pie chart viewer", () => clickOn(page, el("close icon of pie chart viewer")));
      await session.step(60, "Then pie chart viewer should be absent", () => shouldBe(page, el("pie chart viewer"), "absent"));
    });
    await run.scenario("The trellis plot labels a column per core", async () => {
      await session.step(63, "When user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["xColumnNames","Core"]]), [["xColumnNames","Core"]]);
      await session.step(65, "Then the \"x categories\" reading of trellis plot viewer should be 7", () => readingIs(page, "x categories", el("trellis plot viewer"), 7));
      await session.step(66, "And the \"x labels shown\" reading of trellis plot viewer should be 7", () => readingIs(page, "x labels shown", el("trellis plot viewer"), 7));
      await session.step(67, "And no errors should have been logged", () => noErrors(page));
      await session.step(68, "When user clicks on close icon of trellis plot viewer", () => clickOn(page, el("close icon of trellis plot viewer")));
      await session.step(69, "Then trellis plot viewer should be absent", () => shouldBe(page, el("trellis plot viewer"), "absent"));
    });
    await run.scenario("The scatter plot spreads its x axis over the cores, and its x slider moves over the molecules", async () => {
      await session.step(72, "When user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","Core"],["yColumnName","Chemical Space Y"]]), [["xColumnName","Core"],["yColumnName","Chemical Space Y"]]);
      await session.step(75, "Then the \"x axis span\" reading of scatter plot viewer should be 6.6", () => readingIs(page, "x axis span", el("scatter plot viewer"), 6.6));
      await session.step(76, "When user sets \"xColumnName\" property of scatter plot viewer to \"Structure\"", () => setProperty(page, "xColumnName", el("scatter plot viewer"), "Structure"));
      await session.step(77, "Then the \"x axis span\" reading of scatter plot viewer should be at least 90", () => readingAtLeast(page, "x axis span", el("scatter plot viewer"), 90));
      await session.step(78, "When user drags the min handle of the \"x\" range slider of scatter plot viewer by 150 pixels", () => dragRangeHandle(page, "min", "x", el("scatter plot viewer"), 150));
      await session.step(79, "Then the \"x axis min\" reading of scatter plot viewer should be higher than before", () => readingHigher(page, "x axis min", el("scatter plot viewer")));
      await session.step(80, "And scatter plot viewer should have a \"x axis\" area", () => hasArea(page, el("scatter plot viewer"), "x axis"));
      await session.step(81, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
