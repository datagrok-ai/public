/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/category-order/category-order-across-viewers.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.box-plot, viewers.scatter-plot, viewers.trellis-plot, viewers.filters.categorical, viewers.grid]
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
import {clickOn, isExpanded, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnIsCurrentObject, filterPasses, openFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelOpen, contextPanelShows, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, dragAreaToArea, hasNoArea, noErrors, pickFromAreaContextMenu, propertyShouldBe, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {areaLiesBeside} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Category order across viewers", () => {
  const session = feature(test, "features/viewers/category-order/category-order-across-viewers.feature", import.meta.url);
  test("APPLY TO ALL VIEWERS makes the order the column's, and the other viewers follow it", {tag: ["@viewers", "@realizes:viewers.box-plot", "@realizes:viewers.scatter-plot", "@realizes:viewers.trellis-plot", "@realizes:viewers.filters.categorical", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["Value","AGE"],["Category 1","RACE"]]), [["Value","AGE"],["Category 1","RACE"]]);
    await session.step(19, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","RACE"],["Y","AGE"]]), [["X","RACE"],["Y","AGE"]]);
    await session.step(22, "Then the \"category Caucasian\" area of box plot viewer should lie to the left of the \"category Other\" area", () => areaLiesBeside(page, "category Caucasian", el("box plot viewer"), "left", "category Other"));
    await session.step(23, "And the \"marker of row 1\" area of scatter plot viewer should lie to the left of the \"marker of row 2\" area", () => areaLiesBeside(page, "marker of row 1", el("scatter plot viewer"), "left", "marker of row 2"));
    await session.step(26, "Given user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","RACE"],["Viewer Type","Histogram"]]), [["X Column Names","RACE"],["Viewer Type","Histogram"]]);
    await session.step(29, "And user opens the filter panel", () => openFilterPanel(page));
    await session.step(30, "Then the \"categories of RACE\" reading of filter panel should be \"Asian, Black, Caucasian, Other\"", () => readingReads(page, "categories of RACE", el("filter panel"), "Asian, Black, Caucasian, Other"));
    await session.step(31, "And the \"x label Caucasian\" area of trellis plot viewer should lie to the left of the \"x label Other\" area", () => areaLiesBeside(page, "x label Caucasian", el("trellis plot viewer"), "left", "x label Other"));
    await session.step(32, "When user picks \"Order...\" from the context menu of the \"left edge of x axis\" area of scatter plot viewer", () => pickFromAreaContextMenu(page, "Order...", "left edge of x axis", el("scatter plot viewer")));
    await session.step(33, "And user selects \"Row count\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Row count", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
    await session.step(34, "And user selects \"Descending\" in \"Direction\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Descending", el("\"Direction\" input in \"Order categories: RACE\" dialog")));
    await session.step(35, "And user clicks on \"CLOSE\" button in \"Order categories: RACE\" dialog", () => clickOn(page, el("\"CLOSE\" button in \"Order categories: RACE\" dialog")));
    await session.step(36, "Then the \"marker of row 1\" area of scatter plot viewer should lie to the left of the \"marker of row 2\" area", () => areaLiesBeside(page, "marker of row 1", el("scatter plot viewer"), "left", "marker of row 2"));
    await session.step(37, "And the \"marker of row 2\" area of scatter plot viewer should lie to the left of the \"marker of row 10\" area", () => areaLiesBeside(page, "marker of row 2", el("scatter plot viewer"), "left", "marker of row 10"));
    await session.step(38, "When user picks \"Order...\" from the context menu of the \"x axis\" area of box plot viewer", () => pickFromAreaContextMenu(page, "Order...", "x axis", el("box plot viewer")));
    await session.step(39, "And user selects \"Row count\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Row count", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
    await session.step(40, "And user clicks on \"APPLY TO ALL VIEWERS\" button in \"Order categories: RACE\" dialog", () => clickOn(page, el("\"APPLY TO ALL VIEWERS\" button in \"Order categories: RACE\" dialog")));
    await session.step(41, "Then \"APPLY TO ALL VIEWERS\" button in \"Order categories: RACE\" dialog should be disabled", () => shouldBe(page, el("\"APPLY TO ALL VIEWERS\" button in \"Order categories: RACE\" dialog"), "disabled"));
    await session.step(42, "And \"Category Order\" property of box plot viewer should be \"\"", () => propertyShouldBe(page, "Category Order", el("box plot viewer"), ""));
    await session.step(43, "And \"X Category Order\" property of scatter plot viewer should be \"\"", () => propertyShouldBe(page, "X Category Order", el("scatter plot viewer"), ""));
    await session.step(44, "And the \"category Other\" area of box plot viewer should lie to the left of the \"category Caucasian\" area", () => areaLiesBeside(page, "category Other", el("box plot viewer"), "left", "category Caucasian"));
    await session.step(45, "And the \"marker of row 10\" area of scatter plot viewer should lie to the left of the \"marker of row 2\" area", () => areaLiesBeside(page, "marker of row 10", el("scatter plot viewer"), "left", "marker of row 2"));
    await session.step(46, "And the \"marker of row 2\" area of scatter plot viewer should lie to the left of the \"marker of row 1\" area", () => areaLiesBeside(page, "marker of row 2", el("scatter plot viewer"), "left", "marker of row 1"));
    await session.step(47, "And the \"x label Other\" area of trellis plot viewer should lie to the left of the \"x label Caucasian\" area", () => areaLiesBeside(page, "x label Other", el("trellis plot viewer"), "left", "x label Caucasian"));
    await session.step(48, "And the \"x label Black\" area of trellis plot viewer should lie to the left of the \"x label Other\" area", () => areaLiesBeside(page, "x label Black", el("trellis plot viewer"), "left", "x label Other"));
    await session.step(49, "And the \"categories of RACE\" reading of filter panel should be \"Asian, Black, Other, Caucasian\"", () => readingReads(page, "categories of RACE", el("filter panel"), "Asian, Black, Other, Caucasian"));
    await session.step(50, "When user clicks on \"CLOSE\" button in \"Order categories: RACE\" dialog", () => clickOn(page, el("\"CLOSE\" button in \"Order categories: RACE\" dialog")));
    await session.step(51, "Then no errors should have been logged", () => noErrors(page));
  });
  test("A filter keeps the categories it was given when the column is reordered", {tag: ["@viewers", "@realizes:viewers.box-plot", "@realizes:viewers.scatter-plot", "@realizes:viewers.trellis-plot", "@realizes:viewers.filters.categorical", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["Value","AGE"],["Category 1","RACE"]]), [["Value","AGE"],["Category 1","RACE"]]);
    await session.step(19, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","RACE"],["Y","AGE"]]), [["X","RACE"],["Y","AGE"]]);
    await session.step(22, "Then the \"category Caucasian\" area of box plot viewer should lie to the left of the \"category Other\" area", () => areaLiesBeside(page, "category Caucasian", el("box plot viewer"), "left", "category Other"));
    await session.step(23, "And the \"marker of row 1\" area of scatter plot viewer should lie to the left of the \"marker of row 2\" area", () => areaLiesBeside(page, "marker of row 1", el("scatter plot viewer"), "left", "marker of row 2"));
    await session.step(54, "Given user opens the filter panel", () => openFilterPanel(page));
    await session.step(55, "Then the \"selected categories of RACE\" reading of filter panel should be \"Asian, Black, Caucasian, Other\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Asian, Black, Caucasian, Other"));
    await session.step(56, "When user clicks on the \"category Asian of RACE\" area of filter panel", () => clickArea(page, "category Asian of RACE", el("filter panel")));
    await session.step(57, "Then 15 rows should pass the filter", () => filterPasses(page, 15));
    await session.step(58, "When user clicks on the \"checkbox Black of RACE\" area of filter panel", () => clickArea(page, "checkbox Black of RACE", el("filter panel")));
    await session.step(59, "And user clicks on the \"checkbox Other of RACE\" area of filter panel", () => clickArea(page, "checkbox Other of RACE", el("filter panel")));
    await session.step(60, "Then 104 rows should pass the filter", () => filterPasses(page, 104));
    await session.step(61, "And the \"selected categories of RACE\" reading of filter panel should be \"Asian, Black, Other\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Asian, Black, Other"));
    await session.step(62, "When user picks \"Order...\" from the context menu of the \"x axis\" area of box plot viewer", () => pickFromAreaContextMenu(page, "Order...", "x axis", el("box plot viewer")));
    await session.step(63, "And user selects \"Row count\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Row count", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
    await session.step(64, "And user clicks on \"APPLY TO ALL VIEWERS\" button in \"Order categories: RACE\" dialog", () => clickOn(page, el("\"APPLY TO ALL VIEWERS\" button in \"Order categories: RACE\" dialog")));
    await session.step(65, "And user clicks on \"CLOSE\" button in \"Order categories: RACE\" dialog", () => clickOn(page, el("\"CLOSE\" button in \"Order categories: RACE\" dialog")));
    await session.step(66, "Then the \"categories of RACE\" reading of filter panel should be \"Asian, Black, Other, Caucasian\"", () => readingReads(page, "categories of RACE", el("filter panel"), "Asian, Black, Other, Caucasian"));
    await session.step(67, "And the \"selected categories of RACE\" reading of filter panel should be \"Asian, Black, Other\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Asian, Black, Other"));
    await session.step(68, "And 104 rows should pass the filter", () => filterPasses(page, 104));
    await session.step(69, "And box plot viewer should not have a \"category Caucasian\" area", () => hasNoArea(page, el("box plot viewer"), "category Caucasian"));
    await session.step(70, "And the \"category Black\" area of box plot viewer should lie to the left of the \"category Other\" area", () => areaLiesBeside(page, "category Black", el("box plot viewer"), "left", "category Other"));
    await session.step(71, "And no errors should have been logged", () => noErrors(page));
  });
  test("Sort > Custom... in the grid orders the column; a viewer with its own order keeps it", {tag: ["@viewers", "@realizes:viewers.box-plot", "@realizes:viewers.scatter-plot", "@realizes:viewers.trellis-plot", "@realizes:viewers.filters.categorical", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["Value","AGE"],["Category 1","RACE"]]), [["Value","AGE"],["Category 1","RACE"]]);
    await session.step(19, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","RACE"],["Y","AGE"]]), [["X","RACE"],["Y","AGE"]]);
    await session.step(22, "Then the \"category Caucasian\" area of box plot viewer should lie to the left of the \"category Other\" area", () => areaLiesBeside(page, "category Caucasian", el("box plot viewer"), "left", "category Other"));
    await session.step(23, "And the \"marker of row 1\" area of scatter plot viewer should lie to the left of the \"marker of row 2\" area", () => areaLiesBeside(page, "marker of row 1", el("scatter plot viewer"), "left", "marker of row 2"));
    await session.step(74, "When user picks \"Order...\" from the context menu of the \"left edge of x axis\" area of scatter plot viewer", () => pickFromAreaContextMenu(page, "Order...", "left edge of x axis", el("scatter plot viewer")));
    await session.step(75, "And user selects \"Alphabetical\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Alphabetical", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
    await session.step(76, "And user clicks on \"CLOSE\" button in \"Order categories: RACE\" dialog", () => clickOn(page, el("\"CLOSE\" button in \"Order categories: RACE\" dialog")));
    await session.step(77, "And user picks \"Sort > Custom...\" from the context menu of the \"header RACE\" area of grid", () => pickFromAreaContextMenu(page, "Sort > Custom...", "header RACE", el("grid")));
    await session.step(78, "Then \"Order categories\" dialog should be visible", () => shouldBe(page, el("\"Order categories\" dialog"), "visible"));
    await session.step(79, "When user drags the \"cell 4 of RACE\" area of grid in \"Order categories\" dialog to the \"cell 1 of RACE\" area", () => dragAreaToArea(page, "cell 4 of RACE", el("grid in \"Order categories\" dialog"), "cell 1 of RACE"));
    await session.step(80, "And user clicks on OK button in \"Order categories\" dialog", () => clickOn(page, el("OK button in \"Order categories\" dialog")));
    await session.step(81, "Then \"Order categories\" dialog should be absent", () => shouldBe(page, el("\"Order categories\" dialog"), "absent"));
    await session.step(82, "And the \"category Other\" area of box plot viewer should lie to the left of the \"category Asian\" area", () => areaLiesBeside(page, "category Other", el("box plot viewer"), "left", "category Asian"));
    await session.step(83, "And the \"category Asian\" area of box plot viewer should lie to the left of the \"category Black\" area", () => areaLiesBeside(page, "category Asian", el("box plot viewer"), "left", "category Black"));
    await session.step(84, "And the \"marker of row 10\" area of scatter plot viewer should lie to the left of the \"marker of row 1\" area", () => areaLiesBeside(page, "marker of row 10", el("scatter plot viewer"), "left", "marker of row 1"));
    await session.step(85, "And the \"marker of row 1\" area of scatter plot viewer should lie to the left of the \"marker of row 2\" area", () => areaLiesBeside(page, "marker of row 1", el("scatter plot viewer"), "left", "marker of row 2"));
    await session.step(86, "And no errors should have been logged", () => noErrors(page));
  });
  test("CANCEL after APPLY in Sort > Custom... gives the column its default order back", {tag: ["@viewers", "@realizes:viewers.box-plot", "@realizes:viewers.scatter-plot", "@realizes:viewers.trellis-plot", "@realizes:viewers.filters.categorical", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["Value","AGE"],["Category 1","RACE"]]), [["Value","AGE"],["Category 1","RACE"]]);
    await session.step(19, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","RACE"],["Y","AGE"]]), [["X","RACE"],["Y","AGE"]]);
    await session.step(22, "Then the \"category Caucasian\" area of box plot viewer should lie to the left of the \"category Other\" area", () => areaLiesBeside(page, "category Caucasian", el("box plot viewer"), "left", "category Other"));
    await session.step(23, "And the \"marker of row 1\" area of scatter plot viewer should lie to the left of the \"marker of row 2\" area", () => areaLiesBeside(page, "marker of row 1", el("scatter plot viewer"), "left", "marker of row 2"));
    await session.step(89, "Given the \"RACE\" column is the current object", () => columnIsCurrentObject(page, "RACE"));
    await session.step(90, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(91, "And \"Details\" pane in context panel is expanded", () => isExpanded(page, el("\"Details\" pane in context panel")));
    await session.step(92, "Then \"default\" action link in Details pane in context panel should be visible", () => shouldBe(page, el("\"default\" action link in Details pane in context panel"), "visible"));
    await session.step(93, "When user picks \"Sort > Custom...\" from the context menu of the \"header RACE\" area of grid", () => pickFromAreaContextMenu(page, "Sort > Custom...", "header RACE", el("grid")));
    await session.step(94, "And user drags the \"cell 4 of RACE\" area of grid in \"Order categories\" dialog to the \"cell 1 of RACE\" area", () => dragAreaToArea(page, "cell 4 of RACE", el("grid in \"Order categories\" dialog"), "cell 1 of RACE"));
    await session.step(95, "And user clicks on APPLY button in \"Order categories\" dialog", () => clickOn(page, el("APPLY button in \"Order categories\" dialog")));
    await session.step(96, "Then the \"category Other\" area of box plot viewer should lie to the left of the \"category Asian\" area", () => areaLiesBeside(page, "category Other", el("box plot viewer"), "left", "category Asian"));
    await session.step(97, "When user clicks on CANCEL button in \"Order categories\" dialog", () => clickOn(page, el("CANCEL button in \"Order categories\" dialog")));
    await session.step(98, "Then the \"category Asian\" area of box plot viewer should lie to the left of the \"category Other\" area", () => areaLiesBeside(page, "category Asian", el("box plot viewer"), "left", "category Other"));
    await session.step(99, "Given the \"AGE\" column is the current object", () => columnIsCurrentObject(page, "AGE"));
    await session.step(100, "And the \"RACE\" column is the current object", () => columnIsCurrentObject(page, "RACE"));
    await session.step(101, "Then the context panel should show \"RACE\"", () => contextPanelShows(page, "RACE"));
    await session.step(102, "And \"default\" action link in Details pane in context panel should be visible", () => shouldBe(page, el("\"default\" action link in Details pane in context panel"), "visible"));
    await session.step(103, "And no errors should have been logged", () => noErrors(page));
  });
});
