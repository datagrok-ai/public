/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/row-source/row-source-mouse-over-group.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot, viewers.line-chart, viewers.histogram, viewers.bar-chart, viewers.pie-chart, viewers.box-plot, viewers.pc-plot]
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
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, hoverArea, noErrors, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Row Source MouseOverGroup on seven viewers", () => {
  const session = feature(test, "features/viewers/row-source/row-source-mouse-over-group.feature", import.meta.url);
  test("scatter plot is empty until a group of pie chart is hovered, then shows that group within its filter [viewer=scatter plot, source=pie chart, small group=slice Asian, small rows=5, large group=slice Caucasian, large rows=463]", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.line-chart", "@realizes:viewers.histogram", "@realizes:viewers.bar-chart", "@realizes:viewers.pie-chart", "@realizes:viewers.box-plot", "@realizes:viewers.pc-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"],["colorColumnName","RACE"],["filter","${AGE} > 44"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"],["colorColumnName","RACE"],["filter","${AGE} > 44"]]);
    await session.step(21, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","AGE"],["yColumnNames","WEIGHT"],["filter","${AGE} > 44"]]), [["xColumnName","AGE"],["yColumnNames","WEIGHT"],["filter","${AGE} > 44"]]);
    await session.step(25, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["valueColumnName","AGE"],["filter","${AGE} > 44"]]), [["valueColumnName","AGE"],["filter","${AGE} > 44"]]);
    await session.step(28, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["valueColumnName","AGE"],["splitColumnName","RACE"],["filter","${AGE} > 44"]]), [["valueColumnName","AGE"],["splitColumnName","RACE"],["filter","${AGE} > 44"]]);
    await session.step(32, "And user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["categoryColumnName","RACE"],["filter","${AGE} > 44"]]), [["categoryColumnName","RACE"],["filter","${AGE} > 44"]]);
    await session.step(35, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["categoryColumnNames","RACE"],["valueColumnName","AGE"],["filter","${AGE} > 44"]]), [["categoryColumnNames","RACE"],["valueColumnName","AGE"],["filter","${AGE} > 44"]]);
    await session.step(39, "And user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["columnNames","AGE, WEIGHT"],["filter","${AGE} > 44"]]), [["columnNames","AGE, WEIGHT"],["filter","${AGE} > 44"]]);
    await session.step(44, "When user sets \"filter\" property of pie chart viewer to \"\"", () => setProperty(page, "filter", el("pie chart viewer"), ""));
    await session.step(45, "And user sets \"rowSource\" property of scatter plot viewer to \"MouseOverGroup\"", () => setProperty(page, "rowSource", el("scatter plot viewer"), "MouseOverGroup"));
    await session.step(46, "Then scatter plot viewer should show 0 rows", () => showsRows(page, el("scatter plot viewer"), 0));
    await session.step(47, "When user hovers over the \"slice Asian\" area of pie chart viewer", () => hoverArea(page, "slice Asian", el("pie chart viewer")));
    await session.step(48, "Then scatter plot viewer should show 5 rows", () => showsRows(page, el("scatter plot viewer"), 5));
    await session.step(49, "When user hovers over the \"slice Caucasian\" area of pie chart viewer", () => hoverArea(page, "slice Caucasian", el("pie chart viewer")));
    await session.step(50, "Then scatter plot viewer should show 463 rows", () => showsRows(page, el("scatter plot viewer"), 463));
    await session.step(51, "And no errors should have been logged", () => noErrors(page));
  });
  test("line chart is empty until a group of pie chart is hovered, then shows that group within its filter [viewer=line chart, source=pie chart, small group=slice Asian, small rows=5, large group=slice Caucasian, large rows=463]", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.line-chart", "@realizes:viewers.histogram", "@realizes:viewers.bar-chart", "@realizes:viewers.pie-chart", "@realizes:viewers.box-plot", "@realizes:viewers.pc-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"],["colorColumnName","RACE"],["filter","${AGE} > 44"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"],["colorColumnName","RACE"],["filter","${AGE} > 44"]]);
    await session.step(21, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","AGE"],["yColumnNames","WEIGHT"],["filter","${AGE} > 44"]]), [["xColumnName","AGE"],["yColumnNames","WEIGHT"],["filter","${AGE} > 44"]]);
    await session.step(25, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["valueColumnName","AGE"],["filter","${AGE} > 44"]]), [["valueColumnName","AGE"],["filter","${AGE} > 44"]]);
    await session.step(28, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["valueColumnName","AGE"],["splitColumnName","RACE"],["filter","${AGE} > 44"]]), [["valueColumnName","AGE"],["splitColumnName","RACE"],["filter","${AGE} > 44"]]);
    await session.step(32, "And user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["categoryColumnName","RACE"],["filter","${AGE} > 44"]]), [["categoryColumnName","RACE"],["filter","${AGE} > 44"]]);
    await session.step(35, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["categoryColumnNames","RACE"],["valueColumnName","AGE"],["filter","${AGE} > 44"]]), [["categoryColumnNames","RACE"],["valueColumnName","AGE"],["filter","${AGE} > 44"]]);
    await session.step(39, "And user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["columnNames","AGE, WEIGHT"],["filter","${AGE} > 44"]]), [["columnNames","AGE, WEIGHT"],["filter","${AGE} > 44"]]);
    await session.step(44, "When user sets \"filter\" property of pie chart viewer to \"\"", () => setProperty(page, "filter", el("pie chart viewer"), ""));
    await session.step(45, "And user sets \"rowSource\" property of line chart viewer to \"MouseOverGroup\"", () => setProperty(page, "rowSource", el("line chart viewer"), "MouseOverGroup"));
    await session.step(46, "Then line chart viewer should show 0 rows", () => showsRows(page, el("line chart viewer"), 0));
    await session.step(47, "When user hovers over the \"slice Asian\" area of pie chart viewer", () => hoverArea(page, "slice Asian", el("pie chart viewer")));
    await session.step(48, "Then line chart viewer should show 5 rows", () => showsRows(page, el("line chart viewer"), 5));
    await session.step(49, "When user hovers over the \"slice Caucasian\" area of pie chart viewer", () => hoverArea(page, "slice Caucasian", el("pie chart viewer")));
    await session.step(50, "Then line chart viewer should show 463 rows", () => showsRows(page, el("line chart viewer"), 463));
    await session.step(51, "And no errors should have been logged", () => noErrors(page));
  });
  test("histogram is empty until a group of pie chart is hovered, then shows that group within its filter [viewer=histogram, source=pie chart, small group=slice Asian, small rows=5, large group=slice Caucasian, large rows=463]", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.line-chart", "@realizes:viewers.histogram", "@realizes:viewers.bar-chart", "@realizes:viewers.pie-chart", "@realizes:viewers.box-plot", "@realizes:viewers.pc-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"],["colorColumnName","RACE"],["filter","${AGE} > 44"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"],["colorColumnName","RACE"],["filter","${AGE} > 44"]]);
    await session.step(21, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","AGE"],["yColumnNames","WEIGHT"],["filter","${AGE} > 44"]]), [["xColumnName","AGE"],["yColumnNames","WEIGHT"],["filter","${AGE} > 44"]]);
    await session.step(25, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["valueColumnName","AGE"],["filter","${AGE} > 44"]]), [["valueColumnName","AGE"],["filter","${AGE} > 44"]]);
    await session.step(28, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["valueColumnName","AGE"],["splitColumnName","RACE"],["filter","${AGE} > 44"]]), [["valueColumnName","AGE"],["splitColumnName","RACE"],["filter","${AGE} > 44"]]);
    await session.step(32, "And user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["categoryColumnName","RACE"],["filter","${AGE} > 44"]]), [["categoryColumnName","RACE"],["filter","${AGE} > 44"]]);
    await session.step(35, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["categoryColumnNames","RACE"],["valueColumnName","AGE"],["filter","${AGE} > 44"]]), [["categoryColumnNames","RACE"],["valueColumnName","AGE"],["filter","${AGE} > 44"]]);
    await session.step(39, "And user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["columnNames","AGE, WEIGHT"],["filter","${AGE} > 44"]]), [["columnNames","AGE, WEIGHT"],["filter","${AGE} > 44"]]);
    await session.step(44, "When user sets \"filter\" property of pie chart viewer to \"\"", () => setProperty(page, "filter", el("pie chart viewer"), ""));
    await session.step(45, "And user sets \"rowSource\" property of histogram viewer to \"MouseOverGroup\"", () => setProperty(page, "rowSource", el("histogram viewer"), "MouseOverGroup"));
    await session.step(46, "Then histogram viewer should show 0 rows", () => showsRows(page, el("histogram viewer"), 0));
    await session.step(47, "When user hovers over the \"slice Asian\" area of pie chart viewer", () => hoverArea(page, "slice Asian", el("pie chart viewer")));
    await session.step(48, "Then histogram viewer should show 5 rows", () => showsRows(page, el("histogram viewer"), 5));
    await session.step(49, "When user hovers over the \"slice Caucasian\" area of pie chart viewer", () => hoverArea(page, "slice Caucasian", el("pie chart viewer")));
    await session.step(50, "Then histogram viewer should show 463 rows", () => showsRows(page, el("histogram viewer"), 463));
    await session.step(51, "And no errors should have been logged", () => noErrors(page));
  });
  test("bar chart is empty until a group of pie chart is hovered, then shows that group within its filter [viewer=bar chart, source=pie chart, small group=slice Asian, small rows=5, large group=slice Caucasian, large rows=463]", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.line-chart", "@realizes:viewers.histogram", "@realizes:viewers.bar-chart", "@realizes:viewers.pie-chart", "@realizes:viewers.box-plot", "@realizes:viewers.pc-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"],["colorColumnName","RACE"],["filter","${AGE} > 44"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"],["colorColumnName","RACE"],["filter","${AGE} > 44"]]);
    await session.step(21, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","AGE"],["yColumnNames","WEIGHT"],["filter","${AGE} > 44"]]), [["xColumnName","AGE"],["yColumnNames","WEIGHT"],["filter","${AGE} > 44"]]);
    await session.step(25, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["valueColumnName","AGE"],["filter","${AGE} > 44"]]), [["valueColumnName","AGE"],["filter","${AGE} > 44"]]);
    await session.step(28, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["valueColumnName","AGE"],["splitColumnName","RACE"],["filter","${AGE} > 44"]]), [["valueColumnName","AGE"],["splitColumnName","RACE"],["filter","${AGE} > 44"]]);
    await session.step(32, "And user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["categoryColumnName","RACE"],["filter","${AGE} > 44"]]), [["categoryColumnName","RACE"],["filter","${AGE} > 44"]]);
    await session.step(35, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["categoryColumnNames","RACE"],["valueColumnName","AGE"],["filter","${AGE} > 44"]]), [["categoryColumnNames","RACE"],["valueColumnName","AGE"],["filter","${AGE} > 44"]]);
    await session.step(39, "And user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["columnNames","AGE, WEIGHT"],["filter","${AGE} > 44"]]), [["columnNames","AGE, WEIGHT"],["filter","${AGE} > 44"]]);
    await session.step(44, "When user sets \"filter\" property of pie chart viewer to \"\"", () => setProperty(page, "filter", el("pie chart viewer"), ""));
    await session.step(45, "And user sets \"rowSource\" property of bar chart viewer to \"MouseOverGroup\"", () => setProperty(page, "rowSource", el("bar chart viewer"), "MouseOverGroup"));
    await session.step(46, "Then bar chart viewer should show 0 rows", () => showsRows(page, el("bar chart viewer"), 0));
    await session.step(47, "When user hovers over the \"slice Asian\" area of pie chart viewer", () => hoverArea(page, "slice Asian", el("pie chart viewer")));
    await session.step(48, "Then bar chart viewer should show 5 rows", () => showsRows(page, el("bar chart viewer"), 5));
    await session.step(49, "When user hovers over the \"slice Caucasian\" area of pie chart viewer", () => hoverArea(page, "slice Caucasian", el("pie chart viewer")));
    await session.step(50, "Then bar chart viewer should show 463 rows", () => showsRows(page, el("bar chart viewer"), 463));
    await session.step(51, "And no errors should have been logged", () => noErrors(page));
  });
  test("box plot is empty until a group of pie chart is hovered, then shows that group within its filter [viewer=box plot, source=pie chart, small group=slice Asian, small rows=5, large group=slice Caucasian, large rows=463]", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.line-chart", "@realizes:viewers.histogram", "@realizes:viewers.bar-chart", "@realizes:viewers.pie-chart", "@realizes:viewers.box-plot", "@realizes:viewers.pc-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"],["colorColumnName","RACE"],["filter","${AGE} > 44"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"],["colorColumnName","RACE"],["filter","${AGE} > 44"]]);
    await session.step(21, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","AGE"],["yColumnNames","WEIGHT"],["filter","${AGE} > 44"]]), [["xColumnName","AGE"],["yColumnNames","WEIGHT"],["filter","${AGE} > 44"]]);
    await session.step(25, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["valueColumnName","AGE"],["filter","${AGE} > 44"]]), [["valueColumnName","AGE"],["filter","${AGE} > 44"]]);
    await session.step(28, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["valueColumnName","AGE"],["splitColumnName","RACE"],["filter","${AGE} > 44"]]), [["valueColumnName","AGE"],["splitColumnName","RACE"],["filter","${AGE} > 44"]]);
    await session.step(32, "And user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["categoryColumnName","RACE"],["filter","${AGE} > 44"]]), [["categoryColumnName","RACE"],["filter","${AGE} > 44"]]);
    await session.step(35, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["categoryColumnNames","RACE"],["valueColumnName","AGE"],["filter","${AGE} > 44"]]), [["categoryColumnNames","RACE"],["valueColumnName","AGE"],["filter","${AGE} > 44"]]);
    await session.step(39, "And user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["columnNames","AGE, WEIGHT"],["filter","${AGE} > 44"]]), [["columnNames","AGE, WEIGHT"],["filter","${AGE} > 44"]]);
    await session.step(44, "When user sets \"filter\" property of pie chart viewer to \"\"", () => setProperty(page, "filter", el("pie chart viewer"), ""));
    await session.step(45, "And user sets \"rowSource\" property of box plot viewer to \"MouseOverGroup\"", () => setProperty(page, "rowSource", el("box plot viewer"), "MouseOverGroup"));
    await session.step(46, "Then box plot viewer should show 0 rows", () => showsRows(page, el("box plot viewer"), 0));
    await session.step(47, "When user hovers over the \"slice Asian\" area of pie chart viewer", () => hoverArea(page, "slice Asian", el("pie chart viewer")));
    await session.step(48, "Then box plot viewer should show 5 rows", () => showsRows(page, el("box plot viewer"), 5));
    await session.step(49, "When user hovers over the \"slice Caucasian\" area of pie chart viewer", () => hoverArea(page, "slice Caucasian", el("pie chart viewer")));
    await session.step(50, "Then box plot viewer should show 463 rows", () => showsRows(page, el("box plot viewer"), 463));
    await session.step(51, "And no errors should have been logged", () => noErrors(page));
  });
  test("pc plot is empty until a group of pie chart is hovered, then shows that group within its filter [viewer=pc plot, source=pie chart, small group=slice Asian, small rows=5, large group=slice Caucasian, large rows=463]", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.line-chart", "@realizes:viewers.histogram", "@realizes:viewers.bar-chart", "@realizes:viewers.pie-chart", "@realizes:viewers.box-plot", "@realizes:viewers.pc-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"],["colorColumnName","RACE"],["filter","${AGE} > 44"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"],["colorColumnName","RACE"],["filter","${AGE} > 44"]]);
    await session.step(21, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","AGE"],["yColumnNames","WEIGHT"],["filter","${AGE} > 44"]]), [["xColumnName","AGE"],["yColumnNames","WEIGHT"],["filter","${AGE} > 44"]]);
    await session.step(25, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["valueColumnName","AGE"],["filter","${AGE} > 44"]]), [["valueColumnName","AGE"],["filter","${AGE} > 44"]]);
    await session.step(28, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["valueColumnName","AGE"],["splitColumnName","RACE"],["filter","${AGE} > 44"]]), [["valueColumnName","AGE"],["splitColumnName","RACE"],["filter","${AGE} > 44"]]);
    await session.step(32, "And user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["categoryColumnName","RACE"],["filter","${AGE} > 44"]]), [["categoryColumnName","RACE"],["filter","${AGE} > 44"]]);
    await session.step(35, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["categoryColumnNames","RACE"],["valueColumnName","AGE"],["filter","${AGE} > 44"]]), [["categoryColumnNames","RACE"],["valueColumnName","AGE"],["filter","${AGE} > 44"]]);
    await session.step(39, "And user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["columnNames","AGE, WEIGHT"],["filter","${AGE} > 44"]]), [["columnNames","AGE, WEIGHT"],["filter","${AGE} > 44"]]);
    await session.step(44, "When user sets \"filter\" property of pie chart viewer to \"\"", () => setProperty(page, "filter", el("pie chart viewer"), ""));
    await session.step(45, "And user sets \"rowSource\" property of pc plot viewer to \"MouseOverGroup\"", () => setProperty(page, "rowSource", el("pc plot viewer"), "MouseOverGroup"));
    await session.step(46, "Then pc plot viewer should show 0 rows", () => showsRows(page, el("pc plot viewer"), 0));
    await session.step(47, "When user hovers over the \"slice Asian\" area of pie chart viewer", () => hoverArea(page, "slice Asian", el("pie chart viewer")));
    await session.step(48, "Then pc plot viewer should show 5 rows", () => showsRows(page, el("pc plot viewer"), 5));
    await session.step(49, "When user hovers over the \"slice Caucasian\" area of pie chart viewer", () => hoverArea(page, "slice Caucasian", el("pie chart viewer")));
    await session.step(50, "Then pc plot viewer should show 463 rows", () => showsRows(page, el("pc plot viewer"), 463));
    await session.step(51, "And no errors should have been logged", () => noErrors(page));
  });
  test("pie chart is empty until a group of bar chart is hovered, then shows that group within its filter [viewer=pie chart, source=bar chart, small group=bar Asian, small rows=5, large group=bar Caucasian, large rows=463]", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.line-chart", "@realizes:viewers.histogram", "@realizes:viewers.bar-chart", "@realizes:viewers.pie-chart", "@realizes:viewers.box-plot", "@realizes:viewers.pc-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"],["colorColumnName","RACE"],["filter","${AGE} > 44"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"],["colorColumnName","RACE"],["filter","${AGE} > 44"]]);
    await session.step(21, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","AGE"],["yColumnNames","WEIGHT"],["filter","${AGE} > 44"]]), [["xColumnName","AGE"],["yColumnNames","WEIGHT"],["filter","${AGE} > 44"]]);
    await session.step(25, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["valueColumnName","AGE"],["filter","${AGE} > 44"]]), [["valueColumnName","AGE"],["filter","${AGE} > 44"]]);
    await session.step(28, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["valueColumnName","AGE"],["splitColumnName","RACE"],["filter","${AGE} > 44"]]), [["valueColumnName","AGE"],["splitColumnName","RACE"],["filter","${AGE} > 44"]]);
    await session.step(32, "And user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["categoryColumnName","RACE"],["filter","${AGE} > 44"]]), [["categoryColumnName","RACE"],["filter","${AGE} > 44"]]);
    await session.step(35, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["categoryColumnNames","RACE"],["valueColumnName","AGE"],["filter","${AGE} > 44"]]), [["categoryColumnNames","RACE"],["valueColumnName","AGE"],["filter","${AGE} > 44"]]);
    await session.step(39, "And user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["columnNames","AGE, WEIGHT"],["filter","${AGE} > 44"]]), [["columnNames","AGE, WEIGHT"],["filter","${AGE} > 44"]]);
    await session.step(44, "When user sets \"filter\" property of bar chart viewer to \"\"", () => setProperty(page, "filter", el("bar chart viewer"), ""));
    await session.step(45, "And user sets \"rowSource\" property of pie chart viewer to \"MouseOverGroup\"", () => setProperty(page, "rowSource", el("pie chart viewer"), "MouseOverGroup"));
    await session.step(46, "Then pie chart viewer should show 0 rows", () => showsRows(page, el("pie chart viewer"), 0));
    await session.step(47, "When user hovers over the \"bar Asian\" area of bar chart viewer", () => hoverArea(page, "bar Asian", el("bar chart viewer")));
    await session.step(48, "Then pie chart viewer should show 5 rows", () => showsRows(page, el("pie chart viewer"), 5));
    await session.step(49, "When user hovers over the \"bar Caucasian\" area of bar chart viewer", () => hoverArea(page, "bar Caucasian", el("bar chart viewer")));
    await session.step(50, "Then pie chart viewer should show 463 rows", () => showsRows(page, el("pie chart viewer"), 463));
    await session.step(51, "And no errors should have been logged", () => noErrors(page));
  });
});
