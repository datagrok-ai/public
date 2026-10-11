/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/trellis-plot/trellis-plot-use-in-trellis.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.trellis-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, addViewerWith, noErrors, painted, pickFromAreaContextMenu, pickFromContextMenu, readingIs, readingReads, reportsNoError, setProperty, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {innerPropertyShouldBe} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Use in Trellis from another viewer", () => {
  const session = feature(test, "features/viewers/trellis-plot/trellis-plot-use-in-trellis.feature", import.meta.url);
  test("A scatter plot in a trellis keeps its X, Y and color columns", {tag: ["@viewers", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "Given user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","AGE"],["Y","HEIGHT"],["Color","SEX"]]), [["X","AGE"],["Y","HEIGHT"],["Color","SEX"]]);
    await session.step(21, "When user picks \"General > Use in Trellis\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "General > Use in Trellis", el("scatter plot viewer")));
    await session.step(22, "Then the open tableview should have 1 trellis plot viewer", () => viewerCount(page, 1, "trellis plot"));
    await session.step(23, "And the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
    await session.step(24, "And the \"inner viewer type\" reading of trellis plot viewer should be \"Scatter plot\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Scatter plot"));
    await session.step(25, "And \"xColumnName\" inner property of trellis plot viewer should be \"AGE\"", () => innerPropertyShouldBe(page, "xColumnName", el("trellis plot viewer"), "AGE"));
    await session.step(26, "And \"yColumnName\" inner property of trellis plot viewer should be \"HEIGHT\"", () => innerPropertyShouldBe(page, "yColumnName", el("trellis plot viewer"), "HEIGHT"));
    await session.step(27, "And \"colorColumnName\" inner property of trellis plot viewer should be \"SEX\"", () => innerPropertyShouldBe(page, "colorColumnName", el("trellis plot viewer"), "SEX"));
    await session.step(28, "And the \"cells drawn\" reading of trellis plot viewer should be 23", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 23));
    await session.step(29, "And trellis plot viewer should report no error", () => reportsNoError(page, el("trellis plot viewer")));
    await session.step(30, "And no errors should have been logged", () => noErrors(page));
  });
  test("A bar chart in a trellis becomes its inner viewer, its Value carried over [viewer=bar chart, type=Bar chart, area=view, caption=Value, value=WEIGHT, inner=valueColumnName, cells=28]", {tag: ["@viewers", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(33, "Given user adds a bar chart viewer", () => addViewer(page, "bar chart"));
    await session.step(34, "And user sets \"Value\" property of bar chart viewer to \"WEIGHT\"", () => setProperty(page, "Value", el("bar chart viewer"), "WEIGHT"));
    await session.step(35, "Then bar chart viewer should be painted", () => painted(page, el("bar chart viewer")));
    await session.step(36, "When user picks \"General > Use in Trellis\" from the context menu of the \"view\" area of bar chart viewer", () => pickFromAreaContextMenu(page, "General > Use in Trellis", "view", el("bar chart viewer")));
    await session.step(37, "Then the open tableview should have 1 trellis plot viewer", () => viewerCount(page, 1, "trellis plot"));
    await session.step(38, "And the open tableview should have 1 bar chart viewer", () => viewerCount(page, 1, "bar chart"));
    await session.step(39, "And the \"inner viewer type\" reading of trellis plot viewer should be \"Bar chart\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Bar chart"));
    await session.step(40, "And \"valueColumnName\" inner property of trellis plot viewer should be \"WEIGHT\"", () => innerPropertyShouldBe(page, "valueColumnName", el("trellis plot viewer"), "WEIGHT"));
    await session.step(41, "And the \"cells drawn\" reading of trellis plot viewer should be 28", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 28));
    await session.step(42, "And trellis plot viewer should report no error", () => reportsNoError(page, el("trellis plot viewer")));
    await session.step(43, "And no errors should have been logged", () => noErrors(page));
  });
  test("A histogram in a trellis becomes its inner viewer, its Value carried over [viewer=histogram, type=Histogram, area=view, caption=Value, value=HEIGHT, inner=valueColumnName, cells=23]", {tag: ["@viewers", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(33, "Given user adds a histogram viewer", () => addViewer(page, "histogram"));
    await session.step(34, "And user sets \"Value\" property of histogram viewer to \"HEIGHT\"", () => setProperty(page, "Value", el("histogram viewer"), "HEIGHT"));
    await session.step(35, "Then histogram viewer should be painted", () => painted(page, el("histogram viewer")));
    await session.step(36, "When user picks \"General > Use in Trellis\" from the context menu of the \"view\" area of histogram viewer", () => pickFromAreaContextMenu(page, "General > Use in Trellis", "view", el("histogram viewer")));
    await session.step(37, "Then the open tableview should have 1 trellis plot viewer", () => viewerCount(page, 1, "trellis plot"));
    await session.step(38, "And the open tableview should have 1 histogram viewer", () => viewerCount(page, 1, "histogram"));
    await session.step(39, "And the \"inner viewer type\" reading of trellis plot viewer should be \"Histogram\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Histogram"));
    await session.step(40, "And \"valueColumnName\" inner property of trellis plot viewer should be \"HEIGHT\"", () => innerPropertyShouldBe(page, "valueColumnName", el("trellis plot viewer"), "HEIGHT"));
    await session.step(41, "And the \"cells drawn\" reading of trellis plot viewer should be 23", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 23));
    await session.step(42, "And trellis plot viewer should report no error", () => reportsNoError(page, el("trellis plot viewer")));
    await session.step(43, "And no errors should have been logged", () => noErrors(page));
  });
  test("A line chart in a trellis becomes its inner viewer, its X carried over [viewer=line chart, type=Line chart, area=plot, caption=X, value=STARTED, inner=xColumnName, cells=28]", {tag: ["@viewers", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(33, "Given user adds a line chart viewer", () => addViewer(page, "line chart"));
    await session.step(34, "And user sets \"X\" property of line chart viewer to \"STARTED\"", () => setProperty(page, "X", el("line chart viewer"), "STARTED"));
    await session.step(35, "Then line chart viewer should be painted", () => painted(page, el("line chart viewer")));
    await session.step(36, "When user picks \"General > Use in Trellis\" from the context menu of the \"plot\" area of line chart viewer", () => pickFromAreaContextMenu(page, "General > Use in Trellis", "plot", el("line chart viewer")));
    await session.step(37, "Then the open tableview should have 1 trellis plot viewer", () => viewerCount(page, 1, "trellis plot"));
    await session.step(38, "And the open tableview should have 1 line chart viewer", () => viewerCount(page, 1, "line chart"));
    await session.step(39, "And the \"inner viewer type\" reading of trellis plot viewer should be \"Line chart\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Line chart"));
    await session.step(40, "And \"xColumnName\" inner property of trellis plot viewer should be \"STARTED\"", () => innerPropertyShouldBe(page, "xColumnName", el("trellis plot viewer"), "STARTED"));
    await session.step(41, "And the \"cells drawn\" reading of trellis plot viewer should be 28", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 28));
    await session.step(42, "And trellis plot viewer should report no error", () => reportsNoError(page, el("trellis plot viewer")));
    await session.step(43, "And no errors should have been logged", () => noErrors(page));
  });
  test("A box plot in a trellis becomes its inner viewer, its Value carried over [viewer=box plot, type=Box plot, area=view, caption=Value, value=WEIGHT, inner=valueColumnName, cells=28]", {tag: ["@viewers", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(33, "Given user adds a box plot viewer", () => addViewer(page, "box plot"));
    await session.step(34, "And user sets \"Value\" property of box plot viewer to \"WEIGHT\"", () => setProperty(page, "Value", el("box plot viewer"), "WEIGHT"));
    await session.step(35, "Then box plot viewer should be painted", () => painted(page, el("box plot viewer")));
    await session.step(36, "When user picks \"General > Use in Trellis\" from the context menu of the \"view\" area of box plot viewer", () => pickFromAreaContextMenu(page, "General > Use in Trellis", "view", el("box plot viewer")));
    await session.step(37, "Then the open tableview should have 1 trellis plot viewer", () => viewerCount(page, 1, "trellis plot"));
    await session.step(38, "And the open tableview should have 1 box plot viewer", () => viewerCount(page, 1, "box plot"));
    await session.step(39, "And the \"inner viewer type\" reading of trellis plot viewer should be \"Box plot\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Box plot"));
    await session.step(40, "And \"valueColumnName\" inner property of trellis plot viewer should be \"WEIGHT\"", () => innerPropertyShouldBe(page, "valueColumnName", el("trellis plot viewer"), "WEIGHT"));
    await session.step(41, "And the \"cells drawn\" reading of trellis plot viewer should be 28", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 28));
    await session.step(42, "And trellis plot viewer should report no error", () => reportsNoError(page, el("trellis plot viewer")));
    await session.step(43, "And no errors should have been logged", () => noErrors(page));
  });
});
