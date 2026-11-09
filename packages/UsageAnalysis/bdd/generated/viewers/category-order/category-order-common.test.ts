/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/category-order/category-order-common.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.box-plot, viewers.scatter-plot, viewers.bar-chart, viewers.pie-chart, viewers.scatter-plot3d]
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
import {clickOn, selectIn, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, noErrors, pickFromAreaContextMenu, readingBetween, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {areaLies, areaLiesBeside} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("One category order shared by every ordering viewer", () => {
  const session = feature(test, "features/viewers/category-order/category-order-common.feature", import.meta.url);
  test("An order of one viewer moves that viewer only; APPLY TO ALL VIEWERS moves them all, and Reset brings them all back", {tag: ["@viewers", "@realizes:viewers.box-plot", "@realizes:viewers.scatter-plot", "@realizes:viewers.bar-chart", "@realizes:viewers.pie-chart", "@realizes:viewers.scatter-plot3d"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["Value","AGE"],["Category 1","RACE"]]), [["Value","AGE"],["Category 1","RACE"]]);
    await session.step(20, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","RACE"],["Y","AGE"]]), [["X","RACE"],["Y","AGE"]]);
    await session.step(23, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"],["Bar Sort Type","custom"]]), [["Split","RACE"],["Bar Sort Type","custom"]]);
    await session.step(26, "And user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"],["Pie Sort Type","custom"]]), [["Category","RACE"],["Pie Sort Type","custom"]]);
    await session.step(29, "And user adds a 3d scatter plot viewer with:", () => addViewerWith(page, "3d scatter plot", [["X","RACE"]]), [["X","RACE"]]);
    await session.step(31, "Then the \"category Caucasian\" area of box plot viewer should lie to the left of the \"category Other\" area", () => areaLiesBeside(page, "category Caucasian", el("box plot viewer"), "left", "category Other"));
    await session.step(32, "And the \"marker of row 1\" area of scatter plot viewer should lie to the left of the \"marker of row 2\" area", () => areaLiesBeside(page, "marker of row 1", el("scatter plot viewer"), "left", "marker of row 2"));
    await session.step(33, "And the \"bar Black\" area of bar chart viewer should lie below the \"bar Caucasian\" area", () => areaLies(page, "bar Black", el("bar chart viewer"), "below", "bar Caucasian"));
    await session.step(34, "And the \"start angle of Asian\" reading of pie chart viewer should be between -0.01 and 0.01", () => readingBetween(page, "start angle of Asian", el("pie chart viewer"), -0.01, 0.01));
    await session.step(35, "And the \"x labels\" reading of 3d scatter plot viewer should be \"Asian, Black, Caucasian, Other\"", () => readingReads(page, "x labels", el("3d scatter plot viewer"), "Asian, Black, Caucasian, Other"));
    await session.step(38, "When user picks \"Order...\" from the context menu of the \"x axis\" area of box plot viewer", () => pickFromAreaContextMenu(page, "Order...", "x axis", el("box plot viewer")));
    await session.step(39, "And user selects \"Row count\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Row count", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
    await session.step(40, "And user selects \"Descending\" in \"Direction\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Descending", el("\"Direction\" input in \"Order categories: RACE\" dialog")));
    await session.step(41, "Then the \"category Other\" area of box plot viewer should lie to the left of the \"category Black\" area", () => areaLiesBeside(page, "category Other", el("box plot viewer"), "left", "category Black"));
    await session.step(42, "And the \"marker of row 1\" area of scatter plot viewer should lie to the left of the \"marker of row 2\" area", () => areaLiesBeside(page, "marker of row 1", el("scatter plot viewer"), "left", "marker of row 2"));
    await session.step(43, "And the \"bar Black\" area of bar chart viewer should lie below the \"bar Caucasian\" area", () => areaLies(page, "bar Black", el("bar chart viewer"), "below", "bar Caucasian"));
    await session.step(44, "And the \"start angle of Asian\" reading of pie chart viewer should be between -0.01 and 0.01", () => readingBetween(page, "start angle of Asian", el("pie chart viewer"), -0.01, 0.01));
    await session.step(45, "And the \"x labels\" reading of 3d scatter plot viewer should be \"Asian, Black, Caucasian, Other\"", () => readingReads(page, "x labels", el("3d scatter plot viewer"), "Asian, Black, Caucasian, Other"));
    await session.step(46, "When user clicks on \"APPLY TO ALL VIEWERS\" button in \"Order categories: RACE\" dialog", () => clickOn(page, el("\"APPLY TO ALL VIEWERS\" button in \"Order categories: RACE\" dialog")));
    await session.step(47, "Then \"Order categories: RACE\" dialog should contain text \"RACE has its own order\"", () => shouldContainText(page, el("\"Order categories: RACE\" dialog"), "RACE has its own order"));
    await session.step(48, "And the \"category Caucasian\" area of box plot viewer should lie to the left of the \"category Other\" area", () => areaLiesBeside(page, "category Caucasian", el("box plot viewer"), "left", "category Other"));
    await session.step(49, "And the \"category Other\" area of box plot viewer should lie to the left of the \"category Black\" area", () => areaLiesBeside(page, "category Other", el("box plot viewer"), "left", "category Black"));
    await session.step(50, "And the \"category Black\" area of box plot viewer should lie to the left of the \"category Asian\" area", () => areaLiesBeside(page, "category Black", el("box plot viewer"), "left", "category Asian"));
    await session.step(51, "And the \"marker of row 1\" area of scatter plot viewer should lie to the left of the \"marker of row 2\" area", () => areaLiesBeside(page, "marker of row 1", el("scatter plot viewer"), "left", "marker of row 2"));
    await session.step(52, "And the \"marker of row 2\" area of scatter plot viewer should lie to the left of the \"marker of row 10\" area", () => areaLiesBeside(page, "marker of row 2", el("scatter plot viewer"), "left", "marker of row 10"));
    await session.step(53, "And the \"bar Caucasian\" area of bar chart viewer should lie below the \"bar Other\" area", () => areaLies(page, "bar Caucasian", el("bar chart viewer"), "below", "bar Other"));
    await session.step(54, "And the \"bar Other\" area of bar chart viewer should lie below the \"bar Black\" area", () => areaLies(page, "bar Other", el("bar chart viewer"), "below", "bar Black"));
    await session.step(55, "And the \"bar Black\" area of bar chart viewer should lie below the \"bar Asian\" area", () => areaLies(page, "bar Black", el("bar chart viewer"), "below", "bar Asian"));
    await session.step(56, "And the \"start angle of Caucasian\" reading of pie chart viewer should be between -0.01 and 0.01", () => readingBetween(page, "start angle of Caucasian", el("pie chart viewer"), -0.01, 0.01));
    await session.step(57, "And the \"start angle of Asian\" reading of pie chart viewer should be between 354.59 and 354.61", () => readingBetween(page, "start angle of Asian", el("pie chart viewer"), 354.59, 354.61));
    await session.step(58, "And the \"x labels\" reading of 3d scatter plot viewer should be \"Caucasian, Other, Black, Asian\"", () => readingReads(page, "x labels", el("3d scatter plot viewer"), "Caucasian, Other, Black, Asian"));
    await session.step(59, "When user clicks on \"Reset\" action link in \"Order categories: RACE\" dialog", () => clickOn(page, el("\"Reset\" action link in \"Order categories: RACE\" dialog")));
    await session.step(60, "Then the \"category Caucasian\" area of box plot viewer should lie to the left of the \"category Other\" area", () => areaLiesBeside(page, "category Caucasian", el("box plot viewer"), "left", "category Other"));
    await session.step(61, "And the \"category Black\" area of box plot viewer should lie to the left of the \"category Caucasian\" area", () => areaLiesBeside(page, "category Black", el("box plot viewer"), "left", "category Caucasian"));
    await session.step(62, "And the \"bar Black\" area of bar chart viewer should lie below the \"bar Caucasian\" area", () => areaLies(page, "bar Black", el("bar chart viewer"), "below", "bar Caucasian"));
    await session.step(63, "And the \"start angle of Asian\" reading of pie chart viewer should be between -0.01 and 0.01", () => readingBetween(page, "start angle of Asian", el("pie chart viewer"), -0.01, 0.01));
    await session.step(64, "And the \"x labels\" reading of 3d scatter plot viewer should be \"Asian, Black, Caucasian, Other\"", () => readingReads(page, "x labels", el("3d scatter plot viewer"), "Asian, Black, Caucasian, Other"));
    await session.step(65, "When user clicks on \"CLOSE\" button in \"Order categories: RACE\" dialog", () => clickOn(page, el("\"CLOSE\" button in \"Order categories: RACE\" dialog")));
    await session.step(66, "Then no errors should have been logged", () => noErrors(page));
  });
});
