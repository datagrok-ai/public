/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/annotation-regions/annotation-regions-tools-menu.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot, viewers.density-plot, viewers.line-chart, viewers.histogram, viewers.bar-chart]
--- */
import {test} from '@playwright/test';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, closeContextMenu, menuDoesNotList, menuLists, noErrors, openContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Tools menu offers the region items", () => {
  const session = feature(test, "features/viewers/annotation-regions/annotation-regions-tools-menu.feature", import.meta.url);
  test("The scatter plot keeps the Lasso Tool at the top level", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.density-plot", "@realizes:viewers.line-chart", "@realizes:viewers.histogram", "@realizes:viewers.bar-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "Given user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"]]), [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"]]);
    await session.step(17, "When user opens the context menu of scatter plot viewer", () => openContextMenu(page, el("scatter plot viewer")));
    await session.step(18, "Then the open menu should list \"Tools > Show Annotation Regions\"", () => menuLists(page, "Tools > Show Annotation Regions"));
    await session.step(19, "And the open menu should list \"Tools > Draw Annotation Region\"", () => menuLists(page, "Tools > Draw Annotation Region"));
    await session.step(20, "And the open menu should list \"Tools > Formula Lines...\"", () => menuLists(page, "Tools > Formula Lines..."));
    await session.step(21, "And the open menu should list \"Lasso Tool\"", () => menuLists(page, "Lasso Tool"));
    await session.step(22, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(23, "Then no errors should have been logged", () => noErrors(page));
  });
  test("The density plot has the Lasso Tool under Tools", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.density-plot", "@realizes:viewers.line-chart", "@realizes:viewers.histogram", "@realizes:viewers.bar-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(26, "Given user adds a density plot viewer with:", () => addViewerWith(page, "density plot", [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"]]), [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"]]);
    await session.step(29, "When user opens the context menu of density plot viewer", () => openContextMenu(page, el("density plot viewer")));
    await session.step(30, "Then the open menu should list \"Tools > Show Annotation Regions\"", () => menuLists(page, "Tools > Show Annotation Regions"));
    await session.step(31, "And the open menu should list \"Tools > Draw Annotation Region\"", () => menuLists(page, "Tools > Draw Annotation Region"));
    await session.step(32, "And the open menu should list \"Tools > Formula Lines...\"", () => menuLists(page, "Tools > Formula Lines..."));
    await session.step(33, "And the open menu should list \"Tools > Lasso Tool\"", () => menuLists(page, "Tools > Lasso Tool"));
    await session.step(34, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(35, "Then no errors should have been logged", () => noErrors(page));
  });
  test("The line chart has the Lasso Tool under Tools", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.density-plot", "@realizes:viewers.line-chart", "@realizes:viewers.histogram", "@realizes:viewers.bar-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(38, "Given user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","AGE"],["yColumnNames","HEIGHT"]]), [["xColumnName","AGE"],["yColumnNames","HEIGHT"]]);
    await session.step(41, "When user opens the context menu of line chart viewer", () => openContextMenu(page, el("line chart viewer")));
    await session.step(42, "Then the open menu should list \"Tools > Show Annotation Regions\"", () => menuLists(page, "Tools > Show Annotation Regions"));
    await session.step(43, "And the open menu should list \"Tools > Draw Annotation Region\"", () => menuLists(page, "Tools > Draw Annotation Region"));
    await session.step(44, "And the open menu should list \"Tools > Formula Lines...\"", () => menuLists(page, "Tools > Formula Lines..."));
    await session.step(45, "And the open menu should list \"Tools > Lasso Tool\"", () => menuLists(page, "Tools > Lasso Tool"));
    await session.step(46, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(47, "Then no errors should have been logged", () => noErrors(page));
  });
  test("The histogram offers no lasso", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.density-plot", "@realizes:viewers.line-chart", "@realizes:viewers.histogram", "@realizes:viewers.bar-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(50, "Given user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["valueColumnName","AGE"]]), [["valueColumnName","AGE"]]);
    await session.step(52, "When user opens the context menu of histogram viewer", () => openContextMenu(page, el("histogram viewer")));
    await session.step(53, "Then the open menu should list \"Tools > Show Annotation Regions\"", () => menuLists(page, "Tools > Show Annotation Regions"));
    await session.step(54, "And the open menu should list \"Tools > Draw Annotation Region\"", () => menuLists(page, "Tools > Draw Annotation Region"));
    await session.step(55, "And the open menu should list \"Tools > Formula Lines...\"", () => menuLists(page, "Tools > Formula Lines..."));
    await session.step(56, "And the open menu should not list \"Tools > Lasso Tool\"", () => menuDoesNotList(page, "Tools > Lasso Tool"));
    await session.step(57, "And the open menu should not list \"Lasso Tool\"", () => menuDoesNotList(page, "Lasso Tool"));
    await session.step(58, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(59, "Then no errors should have been logged", () => noErrors(page));
  });
  test("The bar chart offers no lasso", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.density-plot", "@realizes:viewers.line-chart", "@realizes:viewers.histogram", "@realizes:viewers.bar-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(62, "Given user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["splitColumnName","RACE"],["valueColumnName","AGE"],["valueAggrType","avg"],["orientation","vertical"]]), [["splitColumnName","RACE"],["valueColumnName","AGE"],["valueAggrType","avg"],["orientation","vertical"]]);
    await session.step(67, "When user opens the context menu of bar chart viewer", () => openContextMenu(page, el("bar chart viewer")));
    await session.step(68, "Then the open menu should list \"Tools > Show Annotation Regions\"", () => menuLists(page, "Tools > Show Annotation Regions"));
    await session.step(69, "And the open menu should list \"Tools > Draw Annotation Region\"", () => menuLists(page, "Tools > Draw Annotation Region"));
    await session.step(70, "And the open menu should list \"Tools > Formula Lines...\"", () => menuLists(page, "Tools > Formula Lines..."));
    await session.step(71, "And the open menu should not list \"Tools > Lasso Tool\"", () => menuDoesNotList(page, "Tools > Lasso Tool"));
    await session.step(72, "And the open menu should not list \"Lasso Tool\"", () => menuDoesNotList(page, "Lasso Tool"));
    await session.step(73, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(74, "Then no errors should have been logged", () => noErrors(page));
  });
  test("The horizontal bar chart offers no lasso either", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.density-plot", "@realizes:viewers.line-chart", "@realizes:viewers.histogram", "@realizes:viewers.bar-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(77, "Given user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["splitColumnName","RACE"],["valueColumnName","AGE"],["valueAggrType","avg"],["orientation","horizontal"]]), [["splitColumnName","RACE"],["valueColumnName","AGE"],["valueAggrType","avg"],["orientation","horizontal"]]);
    await session.step(82, "When user opens the context menu of bar chart viewer", () => openContextMenu(page, el("bar chart viewer")));
    await session.step(83, "Then the open menu should list \"Tools > Show Annotation Regions\"", () => menuLists(page, "Tools > Show Annotation Regions"));
    await session.step(84, "And the open menu should list \"Tools > Draw Annotation Region\"", () => menuLists(page, "Tools > Draw Annotation Region"));
    await session.step(85, "And the open menu should list \"Tools > Formula Lines...\"", () => menuLists(page, "Tools > Formula Lines..."));
    await session.step(86, "And the open menu should not list \"Tools > Lasso Tool\"", () => menuDoesNotList(page, "Tools > Lasso Tool"));
    await session.step(87, "And the open menu should not list \"Lasso Tool\"", () => menuDoesNotList(page, "Lasso Tool"));
    await session.step(88, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(89, "Then no errors should have been logged", () => noErrors(page));
  });
  test("The box plot offers no lasso", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.density-plot", "@realizes:viewers.line-chart", "@realizes:viewers.histogram", "@realizes:viewers.bar-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(92, "Given user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["category1ColumnName","RACE"],["valueColumnName","AGE"]]), [["category1ColumnName","RACE"],["valueColumnName","AGE"]]);
    await session.step(95, "When user opens the context menu of box plot viewer", () => openContextMenu(page, el("box plot viewer")));
    await session.step(96, "Then the open menu should list \"Tools > Show Annotation Regions\"", () => menuLists(page, "Tools > Show Annotation Regions"));
    await session.step(97, "And the open menu should list \"Tools > Draw Annotation Region\"", () => menuLists(page, "Tools > Draw Annotation Region"));
    await session.step(98, "And the open menu should list \"Tools > Formula Lines...\"", () => menuLists(page, "Tools > Formula Lines..."));
    await session.step(99, "And the open menu should not list \"Tools > Lasso Tool\"", () => menuDoesNotList(page, "Tools > Lasso Tool"));
    await session.step(100, "And the open menu should not list \"Lasso Tool\"", () => menuDoesNotList(page, "Lasso Tool"));
    await session.step(101, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(102, "Then no errors should have been logged", () => noErrors(page));
  });
});
