/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/box-plot/box-plot-category-order.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.box-plot]
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
import {pickColumnInDialog} from '../../../bindings/filter-panel.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, isExpanded, selectIn, shouldBe, shouldHaveValue, shouldNotOfferChoice} from '@datagrok-libraries/bdd/bindings/common/steps';
import {renameColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelOpen, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaColor, areaNotColor, closeContextMenu, dragAreaToArea, menuDoesNotList, noErrors, pickFromAreaContextMenu, rightClickArea, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {areaLiesBeside} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Box plot category order", () => {
  const session = feature(test, "features/viewers/box-plot/box-plot-category-order.feature", import.meta.url);
  test("Box plot category order", {tag: ["@journey", "@viewers", "@realizes:viewers.box-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(23, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["Value","AGE"],["Category 1","RACE"]]), [["Value","AGE"],["Category 1","RACE"]]);
    await session.step(26, "Then the \"category Asian\" area of box plot viewer should lie to the left of the \"category Black\" area", () => areaLiesBeside(page, "category Asian", el("box plot viewer"), "left", "category Black"));
    await session.step(27, "And the \"category Caucasian\" area of box plot viewer should lie to the left of the \"category Other\" area", () => areaLiesBeside(page, "category Caucasian", el("box plot viewer"), "left", "category Other"));
    await run.scenario("Order... is on the category axis only, and opens on the column's own order", async () => {
      await session.step(30, "When user right-clicks on the \"y axis\" area of box plot viewer", () => rightClickArea(page, "y axis", el("box plot viewer")));
      await session.step(31, "Then the open menu should not list \"Order...\"", () => menuDoesNotList(page, "Order..."));
      await session.step(32, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(33, "And user picks \"Order...\" from the context menu of the \"x axis\" area of box plot viewer", () => pickFromAreaContextMenu(page, "Order...", "x axis", el("box plot viewer")));
      await session.step(34, "Then \"Order categories: RACE\" dialog should be visible", () => shouldBe(page, el("\"Order categories: RACE\" dialog"), "visible"));
      await session.step(35, "And \"Order by\" input in \"Order categories: RACE\" dialog should have value \"Default\"", () => shouldHaveValue(page, el("\"Order by\" input in \"Order categories: RACE\" dialog"), "Default"));
      await session.step(36, "And \"Direction\" input in \"Order categories: RACE\" dialog should be hidden", () => shouldBe(page, el("\"Direction\" input in \"Order categories: RACE\" dialog"), "hidden"));
      await session.step(37, "And \"APPLY TO ALL VIEWERS\" button in \"Order categories: RACE\" dialog should be disabled", () => shouldBe(page, el("\"APPLY TO ALL VIEWERS\" button in \"Order categories: RACE\" dialog"), "disabled"));
      await session.step(38, "When user clicks on \"CLOSE\" button in \"Order categories: RACE\" dialog", () => clickOn(page, el("\"CLOSE\" button in \"Order categories: RACE\" dialog")));
      await session.step(39, "Then \"Order categories: RACE\" dialog should be absent", () => shouldBe(page, el("\"Order categories: RACE\" dialog"), "absent"));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Row count in both directions, then Alphabetical, previewed at once; CLOSE keeps the order", async () => {
      await session.step(43, "When user picks \"Order...\" from the context menu of the \"x axis\" area of box plot viewer", () => pickFromAreaContextMenu(page, "Order...", "x axis", el("box plot viewer")));
      await session.step(44, "And user selects \"Row count\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Row count", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
      await session.step(45, "Then the \"category Other\" area of box plot viewer should lie to the left of the \"category Caucasian\" area", () => areaLiesBeside(page, "category Other", el("box plot viewer"), "left", "category Caucasian"));
      await session.step(46, "And the \"category Black\" area of box plot viewer should lie to the left of the \"category Other\" area", () => areaLiesBeside(page, "category Black", el("box plot viewer"), "left", "category Other"));
      await session.step(47, "And \"Direction\" input in \"Order categories: RACE\" dialog should be visible", () => shouldBe(page, el("\"Direction\" input in \"Order categories: RACE\" dialog"), "visible"));
      await session.step(48, "And \"APPLY TO ALL VIEWERS\" button in \"Order categories: RACE\" dialog should be enabled", () => shouldBe(page, el("\"APPLY TO ALL VIEWERS\" button in \"Order categories: RACE\" dialog"), "enabled"));
      await session.step(49, "When user selects \"Descending\" in \"Direction\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Descending", el("\"Direction\" input in \"Order categories: RACE\" dialog")));
      await session.step(50, "Then the \"category Caucasian\" area of box plot viewer should lie to the left of the \"category Other\" area", () => areaLiesBeside(page, "category Caucasian", el("box plot viewer"), "left", "category Other"));
      await session.step(51, "And the \"category Black\" area of box plot viewer should lie to the left of the \"category Asian\" area", () => areaLiesBeside(page, "category Black", el("box plot viewer"), "left", "category Asian"));
      await session.step(52, "When user selects \"Alphabetical\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Alphabetical", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
      await session.step(53, "Then the \"category Other\" area of box plot viewer should lie to the left of the \"category Caucasian\" area", () => areaLiesBeside(page, "category Other", el("box plot viewer"), "left", "category Caucasian"));
      await session.step(54, "And the \"category Black\" area of box plot viewer should lie to the left of the \"category Asian\" area", () => areaLiesBeside(page, "category Black", el("box plot viewer"), "left", "category Asian"));
      await session.step(55, "When user clicks on \"CLOSE\" button in \"Order categories: RACE\" dialog", () => clickOn(page, el("\"CLOSE\" button in \"Order categories: RACE\" dialog")));
      await session.step(56, "Then the \"category Other\" area of box plot viewer should lie to the left of the \"category Caucasian\" area", () => areaLiesBeside(page, "category Other", el("box plot viewer"), "left", "category Caucasian"));
      await session.step(57, "When user picks \"Order...\" from the context menu of the \"x axis\" area of box plot viewer", () => pickFromAreaContextMenu(page, "Order...", "x axis", el("box plot viewer")));
      await session.step(58, "Then \"Order by\" input in \"Order categories: RACE\" dialog should have value \"Alphabetical\"", () => shouldHaveValue(page, el("\"Order by\" input in \"Order categories: RACE\" dialog"), "Alphabetical"));
      await session.step(59, "And \"Direction\" input in \"Order categories: RACE\" dialog should have value \"Descending\"", () => shouldHaveValue(page, el("\"Direction\" input in \"Order categories: RACE\" dialog"), "Descending"));
      await session.step(60, "When user selects \"Ascending\" in \"Direction\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Ascending", el("\"Direction\" input in \"Order categories: RACE\" dialog")));
      await session.step(61, "And user selects \"Default\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Default", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
      await session.step(62, "Then the \"category Asian\" area of box plot viewer should lie to the left of the \"category Black\" area", () => areaLiesBeside(page, "category Asian", el("box plot viewer"), "left", "category Black"));
      await session.step(63, "And the \"category Caucasian\" area of box plot viewer should lie to the left of the \"category Other\" area", () => areaLiesBeside(page, "category Caucasian", el("box plot viewer"), "left", "category Other"));
      await session.step(64, "And \"APPLY TO ALL VIEWERS\" button in \"Order categories: RACE\" dialog should be disabled", () => shouldBe(page, el("\"APPLY TO ALL VIEWERS\" button in \"Order categories: RACE\" dialog"), "disabled"));
      await session.step(65, "When user clicks on \"CLOSE\" button in \"Order categories: RACE\" dialog", () => clickOn(page, el("\"CLOSE\" button in \"Order categories: RACE\" dialog")));
      await session.step(66, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Aggregated value orders the boxes by another column's aggregate", async () => {
      await session.step(69, "When user picks \"Order...\" from the context menu of the \"x axis\" area of box plot viewer", () => pickFromAreaContextMenu(page, "Order...", "x axis", el("box plot viewer")));
      await session.step(70, "And user selects \"Aggregated value\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Aggregated value", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
      await session.step(71, "Then \"Column\" input in \"Order categories: RACE\" dialog should be visible", () => shouldBe(page, el("\"Column\" input in \"Order categories: RACE\" dialog"), "visible"));
      await session.step(72, "And \"Aggregation\" input in \"Order categories: RACE\" dialog should be visible", () => shouldBe(page, el("\"Aggregation\" input in \"Order categories: RACE\" dialog"), "visible"));
      await session.step(73, "When user picks \"WEIGHT\" in the column selector of the \"Order categories: RACE\" dialog", () => pickColumnInDialog(page, "WEIGHT", "Order categories: RACE"));
      await session.step(74, "And user selects \"Avg\" in \"Aggregation\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Avg", el("\"Aggregation\" input in \"Order categories: RACE\" dialog")));
      await session.step(75, "And user selects \"Ascending\" in \"Direction\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Ascending", el("\"Direction\" input in \"Order categories: RACE\" dialog")));
      await session.step(76, "Then the \"category Asian\" area of box plot viewer should lie to the left of the \"category Other\" area", () => areaLiesBeside(page, "category Asian", el("box plot viewer"), "left", "category Other"));
      await session.step(77, "And the \"category Other\" area of box plot viewer should lie to the left of the \"category Caucasian\" area", () => areaLiesBeside(page, "category Other", el("box plot viewer"), "left", "category Caucasian"));
      await session.step(78, "And the \"category Caucasian\" area of box plot viewer should lie to the left of the \"category Black\" area", () => areaLiesBeside(page, "category Caucasian", el("box plot viewer"), "left", "category Black"));
      await session.step(79, "When user selects \"Descending\" in \"Direction\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Descending", el("\"Direction\" input in \"Order categories: RACE\" dialog")));
      await session.step(80, "Then the \"category Black\" area of box plot viewer should lie to the left of the \"category Caucasian\" area", () => areaLiesBeside(page, "category Black", el("box plot viewer"), "left", "category Caucasian"));
      await session.step(81, "And the \"category Other\" area of box plot viewer should lie to the left of the \"category Asian\" area", () => areaLiesBeside(page, "category Other", el("box plot viewer"), "left", "category Asian"));
      await session.step(82, "When user selects \"Default\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Default", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
      await session.step(83, "Then the \"category Asian\" area of box plot viewer should lie to the left of the \"category Black\" area", () => areaLiesBeside(page, "category Asian", el("box plot viewer"), "left", "category Black"));
      await session.step(84, "When user clicks on \"CLOSE\" button in \"Order categories: RACE\" dialog", () => clickOn(page, el("\"CLOSE\" button in \"Order categories: RACE\" dialog")));
      await session.step(85, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Manual order: a category dragged to the top of the list becomes the first box", async () => {
      await session.step(88, "When user picks \"Order...\" from the context menu of the \"x axis\" area of box plot viewer", () => pickFromAreaContextMenu(page, "Order...", "x axis", el("box plot viewer")));
      await session.step(89, "And user selects \"Manual\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Manual", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
      await session.step(90, "Then \"Direction\" input in \"Order categories: RACE\" dialog should be hidden", () => shouldBe(page, el("\"Direction\" input in \"Order categories: RACE\" dialog"), "hidden"));
      await session.step(91, "And grid in \"Order categories: RACE\" dialog should be visible", () => shouldBe(page, el("grid in \"Order categories: RACE\" dialog"), "visible"));
      await session.step(92, "When user drags the \"cell 4 of RACE\" area of grid in \"Order categories: RACE\" dialog to the \"cell 1 of RACE\" area", () => dragAreaToArea(page, "cell 4 of RACE", el("grid in \"Order categories: RACE\" dialog"), "cell 1 of RACE"));
      await session.step(93, "Then the \"category Other\" area of box plot viewer should lie to the left of the \"category Asian\" area", () => areaLiesBeside(page, "category Other", el("box plot viewer"), "left", "category Asian"));
      await session.step(94, "And the \"category Asian\" area of box plot viewer should lie to the left of the \"category Black\" area", () => areaLiesBeside(page, "category Asian", el("box plot viewer"), "left", "category Black"));
      await session.step(95, "And \"APPLY TO ALL VIEWERS\" button in \"Order categories: RACE\" dialog should be enabled", () => shouldBe(page, el("\"APPLY TO ALL VIEWERS\" button in \"Order categories: RACE\" dialog"), "enabled"));
      await session.step(96, "When user selects \"Default\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Default", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
      await session.step(97, "Then the \"category Asian\" area of box plot viewer should lie to the left of the \"category Other\" area", () => areaLiesBeside(page, "category Asian", el("box plot viewer"), "left", "category Other"));
      await session.step(98, "When user clicks on \"CLOSE\" button in \"Order categories: RACE\" dialog", () => clickOn(page, el("\"CLOSE\" button in \"Order categories: RACE\" dialog")));
      await session.step(99, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Category Order property opens the same dialog from the context panel", async () => {
      await session.step(102, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(103, "When user clicks on settings icon of box plot viewer", () => clickOn(page, el("settings icon of box plot viewer")));
      await session.step(104, "And \"X Axis\" category in context panel is expanded", () => isExpanded(page, el("\"X Axis\" category in context panel")));
      await session.step(105, "Then \"Category Order\" property in context panel should be visible", () => shouldBe(page, el("\"Category Order\" property in context panel"), "visible"));
      await session.step(106, "And \"Category Order\" property in context panel should have value \"\"", () => shouldHaveValue(page, el("\"Category Order\" property in context panel"), ""));
      await session.step(107, "When user clicks on button in \"Category Order\" property in context panel", () => clickOn(page, el("button in \"Category Order\" property in context panel")));
      await session.step(108, "Then \"Order categories: RACE\" dialog should be visible", () => shouldBe(page, el("\"Order categories: RACE\" dialog"), "visible"));
      await session.step(109, "When user selects \"Row count\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Row count", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
      await session.step(110, "Then the \"category Other\" area of box plot viewer should lie to the left of the \"category Caucasian\" area", () => areaLiesBeside(page, "category Other", el("box plot viewer"), "left", "category Caucasian"));
      await session.step(111, "And \"Category Order\" property in context panel should have value '{\"column\":\"RACE\",\"aggr\":\"values\",\"property\":\"\",\"order\":\"asc\",\"manual\":[]}'", () => shouldHaveValue(page, el("\"Category Order\" property in context panel"), "{\"column\":\"RACE\",\"aggr\":\"values\",\"property\":\"\",\"order\":\"asc\",\"manual\":[]}"));
      await session.step(112, "When user selects \"Default\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Default", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
      await session.step(113, "And user clicks on \"CLOSE\" button in \"Order categories: RACE\" dialog", () => clickOn(page, el("\"CLOSE\" button in \"Order categories: RACE\" dialog")));
      await session.step(114, "Then the \"category Caucasian\" area of box plot viewer should lie to the left of the \"category Other\" area", () => areaLiesBeside(page, "category Caucasian", el("box plot viewer"), "left", "category Other"));
      await session.step(115, "And \"Category Order\" property in context panel should have value \"\"", () => shouldHaveValue(page, el("\"Category Order\" property in context panel"), ""));
      await session.step(116, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An order by another column survives that column's rename", async () => {
      await session.step(119, "When user picks \"Order...\" from the context menu of the \"x axis\" area of box plot viewer", () => pickFromAreaContextMenu(page, "Order...", "x axis", el("box plot viewer")));
      await session.step(120, "And user selects \"Aggregated value\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Aggregated value", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
      await session.step(121, "And user picks \"WEIGHT\" in the column selector of the \"Order categories: RACE\" dialog", () => pickColumnInDialog(page, "WEIGHT", "Order categories: RACE"));
      await session.step(122, "And user selects \"Avg\" in \"Aggregation\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Avg", el("\"Aggregation\" input in \"Order categories: RACE\" dialog")));
      await session.step(123, "And user selects \"Ascending\" in \"Direction\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Ascending", el("\"Direction\" input in \"Order categories: RACE\" dialog")));
      await session.step(124, "And user clicks on \"CLOSE\" button in \"Order categories: RACE\" dialog", () => clickOn(page, el("\"CLOSE\" button in \"Order categories: RACE\" dialog")));
      await session.step(125, "Then the \"category Caucasian\" area of box plot viewer should lie to the left of the \"category Black\" area", () => areaLiesBeside(page, "category Caucasian", el("box plot viewer"), "left", "category Black"));
      await session.step(126, "When user renames \"WEIGHT\" column to \"BODY WEIGHT\"", () => renameColumn(page, "WEIGHT", "BODY WEIGHT"));
      await session.step(127, "Then the \"category Caucasian\" area of box plot viewer should lie to the left of the \"category Black\" area", () => areaLiesBeside(page, "category Caucasian", el("box plot viewer"), "left", "category Black"));
      await session.step(128, "And the \"category Asian\" area of box plot viewer should lie to the left of the \"category Other\" area", () => areaLiesBeside(page, "category Asian", el("box plot viewer"), "left", "category Other"));
      await session.step(129, "When user renames \"BODY WEIGHT\" column to \"WEIGHT\"", () => renameColumn(page, "BODY WEIGHT", "WEIGHT"));
      await session.step(130, "And user picks \"Order...\" from the context menu of the \"x axis\" area of box plot viewer", () => pickFromAreaContextMenu(page, "Order...", "x axis", el("box plot viewer")));
      await session.step(131, "And user selects \"Default\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Default", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
      await session.step(132, "And user clicks on \"CLOSE\" button in \"Order categories: RACE\" dialog", () => clickOn(page, el("\"CLOSE\" button in \"Order categories: RACE\" dialog")));
      await session.step(133, "Then the \"category Black\" area of box plot viewer should lie to the left of the \"category Caucasian\" area", () => areaLiesBeside(page, "category Black", el("box plot viewer"), "left", "category Caucasian"));
      await session.step(134, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Bin colors stay with their boxes when the boxes are reordered", async () => {
      await session.step(137, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["binColorColumnName","WEIGHT"],["binColorAggrType","avg"]]), [["binColorColumnName","WEIGHT"],["binColorAggrType","avg"]]);
      await session.step(140, "Then the \"Other values\" area of box plot viewer should contain the color \"#6868ff\"", () => areaColor(page, "Other values", el("box plot viewer"), "#6868ff"));
      await session.step(141, "And the \"Caucasian values\" area of box plot viewer should not contain the color \"#6868ff\"", () => areaNotColor(page, "Caucasian values", el("box plot viewer"), "#6868ff"));
      await session.step(142, "When user picks \"Order...\" from the context menu of the \"x axis\" area of box plot viewer", () => pickFromAreaContextMenu(page, "Order...", "x axis", el("box plot viewer")));
      await session.step(143, "And user selects \"Row count\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Row count", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
      await session.step(144, "Then the \"category Other\" area of box plot viewer should lie to the left of the \"category Caucasian\" area", () => areaLiesBeside(page, "category Other", el("box plot viewer"), "left", "category Caucasian"));
      await session.step(145, "And the \"Other values\" area of box plot viewer should contain the color \"#6868ff\"", () => areaColor(page, "Other values", el("box plot viewer"), "#6868ff"));
      await session.step(146, "And the \"Caucasian values\" area of box plot viewer should not contain the color \"#6868ff\"", () => areaNotColor(page, "Caucasian values", el("box plot viewer"), "#6868ff"));
      await session.step(147, "When user selects \"Default\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Default", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
      await session.step(148, "And user clicks on \"CLOSE\" button in \"Order categories: RACE\" dialog", () => clickOn(page, el("\"CLOSE\" button in \"Order categories: RACE\" dialog")));
      await session.step(149, "And user sets \"binColorColumnName\" property of box plot viewer to \"\"", () => setProperty(page, "binColorColumnName", el("box plot viewer"), ""));
      await session.step(150, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With two category levels the boxes are ordered inside each outer group", async () => {
      await session.step(153, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Category 1","SEX"],["Category 2","RACE"]]), [["Category 1","SEX"],["Category 2","RACE"]]);
      await session.step(156, "Then the \"category F | Black\" area of box plot viewer should lie to the left of the \"category F | Caucasian\" area", () => areaLiesBeside(page, "category F | Black", el("box plot viewer"), "left", "category F | Caucasian"));
      await session.step(157, "When user picks \"Order...\" from the context menu of the \"x axis\" area of box plot viewer", () => pickFromAreaContextMenu(page, "Order...", "x axis", el("box plot viewer")));
      await session.step(158, "Then \"Order categories\" dialog should be visible", () => shouldBe(page, el("\"Order categories\" dialog"), "visible"));
      await session.step(159, "And \"Order by\" input in \"Order categories\" dialog should not offer the choice \"Manual\"", () => shouldNotOfferChoice(page, el("\"Order by\" input in \"Order categories\" dialog"), "Manual"));
      await session.step(160, "When user selects \"Aggregated value\" in \"Order by\" input in \"Order categories\" dialog", () => selectIn(page, "Aggregated value", el("\"Order by\" input in \"Order categories\" dialog")));
      await session.step(161, "Then the \"category F | Caucasian\" area of box plot viewer should lie to the left of the \"category F | Other\" area", () => areaLiesBeside(page, "category F | Caucasian", el("box plot viewer"), "left", "category F | Other"));
      await session.step(162, "And the \"category F | Other\" area of box plot viewer should lie to the left of the \"category F | Black\" area", () => areaLiesBeside(page, "category F | Other", el("box plot viewer"), "left", "category F | Black"));
      await session.step(163, "And the \"category F | Black\" area of box plot viewer should lie to the left of the \"category M | Asian\" area", () => areaLiesBeside(page, "category F | Black", el("box plot viewer"), "left", "category M | Asian"));
      await session.step(164, "And the \"category M | Black\" area of box plot viewer should lie to the left of the \"category M | Caucasian\" area", () => areaLiesBeside(page, "category M | Black", el("box plot viewer"), "left", "category M | Caucasian"));
      await session.step(165, "When user selects \"Default\" in \"Order by\" input in \"Order categories\" dialog", () => selectIn(page, "Default", el("\"Order by\" input in \"Order categories\" dialog")));
      await session.step(166, "And user clicks on \"CLOSE\" button in \"Order categories\" dialog", () => clickOn(page, el("\"CLOSE\" button in \"Order categories\" dialog")));
      await session.step(167, "Then the \"category F | Black\" area of box plot viewer should lie to the left of the \"category F | Caucasian\" area", () => areaLiesBeside(page, "category F | Black", el("box plot viewer"), "left", "category F | Caucasian"));
      await session.step(168, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Category 1","RACE"],["Category 2",""]]), [["Category 1","RACE"],["Category 2",""]]);
      await session.step(171, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
