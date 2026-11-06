/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tile-viewer/tile-viewer-lanes.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tile-viewer]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/file-formats.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {addTileViewer, dragCardIntoLane} from '../../../bindings/tile-viewer.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {filterOut, filterPasses, filterPassesAll, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, hasArea, hasNoArea, menuLists, noErrors, pointerAway, propertyShouldBe, readingAtLeast, readingIs, readingReads, setProperties, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {openViewerMenu, pickFromViewerMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Tile viewer lanes", () => {
  const session = feature(test, "features/viewers/tile-viewer/tile-viewer-lanes.feature", import.meta.url);
  test("Tile viewer lanes", {tag: ["@journey", "@viewers", "@realizes:viewers.tile-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 11, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(22, "And user adds a tile viewer", () => addTileViewer(page));
    await session.step(23, "Then tile viewer should be visible", () => shouldBe(page, el("tile viewer"), "visible"));
    await session.step(24, "And tile viewer should show 1000 rows", () => showsRows(page, el("tile viewer"), 1000));
    await run.scenario("With no lanes column there is one unnamed lane and no region for the viewer's menu", async () => {
      await session.step(27, "Then the \"lanes\" reading of tile viewer should be 1", () => readingIs(page, "lanes", el("tile viewer"), 1));
      await session.step(28, "And the \"lane names\" reading of tile viewer should be \"All rows\"", () => readingReads(page, "lane names", el("tile viewer"), "All rows"));
      await session.step(29, "And the \"single lane\" reading of tile viewer should be \"true\"", () => readingReads(page, "single lane", el("tile viewer"), "true"));
      await session.step(30, "And the \"lanes list\" reading of tile viewer should be \"\"", () => readingReads(page, "lanes list", el("tile viewer"), ""));
      await session.step(31, "And tile viewer should have a \"lane All rows\" area", () => hasArea(page, el("tile viewer"), "lane All rows"));
      await session.step(32, "And tile viewer should have a \"lane content All rows\" area", () => hasArea(page, el("tile viewer"), "lane content All rows"));
      await session.step(33, "And tile viewer should not have a \"lane header All rows\" area", () => hasNoArea(page, el("tile viewer"), "lane header All rows"));
      await session.step(34, "And tile viewer should not have a \"viewer menu\" area", () => hasNoArea(page, el("tile viewer"), "viewer menu"));
      await session.step(35, "And the \"tiles\" reading of tile viewer should be at least 3", () => readingAtLeast(page, "tiles", el("tile viewer"), 3));
      await session.step(36, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Lanes = RACE gives a lane per category, in the column's category order", async () => {
      await session.step(39, "When user sets \"Lanes Column Name\" property of tile viewer to \"RACE\"", () => setProperty(page, "Lanes Column Name", el("tile viewer"), "RACE"));
      await session.step(40, "Then the \"lanes\" reading of tile viewer should be 4", () => readingIs(page, "lanes", el("tile viewer"), 4));
      await session.step(41, "And the \"lane names\" reading of tile viewer should be \"Asian, Black, Caucasian, Other\"", () => readingReads(page, "lane names", el("tile viewer"), "Asian, Black, Caucasian, Other"));
      await session.step(42, "And the \"single lane\" reading of tile viewer should be \"false\"", () => readingReads(page, "single lane", el("tile viewer"), "false"));
      await session.step(43, "And the \"lanes list\" reading of tile viewer should be \"\"", () => readingReads(page, "lanes list", el("tile viewer"), ""));
      await session.step(44, "And tile viewer should have a \"lane header Asian\" area", () => hasArea(page, el("tile viewer"), "lane header Asian"));
      await session.step(45, "And tile viewer should have a \"lane header Caucasian\" area", () => hasArea(page, el("tile viewer"), "lane header Caucasian"));
      await session.step(46, "And tile viewer should have a \"viewer menu\" area", () => hasArea(page, el("tile viewer"), "viewer menu"));
      await session.step(47, "And the \"lane of row 1\" reading of tile viewer should be \"Caucasian\"", () => readingReads(page, "lane of row 1", el("tile viewer"), "Caucasian"));
      await session.step(48, "And the \"lane of row 2\" reading of tile viewer should be \"Other\"", () => readingReads(page, "lane of row 2", el("tile viewer"), "Other"));
      await session.step(49, "And the \"lane of row 10\" reading of tile viewer should be \"Asian\"", () => readingReads(page, "lane of row 10", el("tile viewer"), "Asian"));
      await session.step(50, "And the \"lane of row 101\" reading of tile viewer should be \"Black\"", () => readingReads(page, "lane of row 101", el("tile viewer"), "Black"));
      await session.step(51, "And the \"tiles in lane Asian\" reading of tile viewer should be at least 1", () => readingAtLeast(page, "tiles in lane Asian", el("tile viewer"), 1));
      await session.step(52, "And the \"tiles in lane Black\" reading of tile viewer should be at least 1", () => readingAtLeast(page, "tiles in lane Black", el("tile viewer"), 1));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Lanes = SEX gives the two lanes of that column", async () => {
      await session.step(56, "When user sets \"Lanes Column Name\" property of tile viewer to \"SEX\"", () => setProperty(page, "Lanes Column Name", el("tile viewer"), "SEX"));
      await session.step(57, "Then the \"lanes\" reading of tile viewer should be 2", () => readingIs(page, "lanes", el("tile viewer"), 2));
      await session.step(58, "And the \"lane names\" reading of tile viewer should be \"F, M\"", () => readingReads(page, "lane names", el("tile viewer"), "F, M"));
      await session.step(59, "And tile viewer should have a \"lane header F\" area", () => hasArea(page, el("tile viewer"), "lane header F"));
      await session.step(60, "And tile viewer should have a \"lane header M\" area", () => hasArea(page, el("tile viewer"), "lane header M"));
      await session.step(61, "And tile viewer should not have a \"lane Caucasian\" area", () => hasNoArea(page, el("tile viewer"), "lane Caucasian"));
      await session.step(62, "And the \"lane of row 1\" reading of tile viewer should be \"F\"", () => readingReads(page, "lane of row 1", el("tile viewer"), "F"));
      await session.step(63, "And the \"lane of row 4\" reading of tile viewer should be \"M\"", () => readingReads(page, "lane of row 4", el("tile viewer"), "M"));
      await session.step(64, "When user clicks on settings icon of tile viewer", () => clickOn(page, el("settings icon of tile viewer")));
      await session.step(65, "And user selects \"RACE\" in \"Lanes\" property in context panel", () => selectIn(page, "RACE", el("\"Lanes\" property in context panel")));
      await session.step(66, "Then \"Lanes Column Name\" property of tile viewer should be \"RACE\"", () => propertyShouldBe(page, "Lanes Column Name", el("tile viewer"), "RACE"));
      await session.step(67, "And the \"lane names\" reading of tile viewer should be \"Asian, Black, Caucasian, Other\"", () => readingReads(page, "lane names", el("tile viewer"), "Asian, Black, Caucasian, Other"));
      await session.step(68, "When user selects \"SEX\" in \"Lanes\" property in context panel", () => selectIn(page, "SEX", el("\"Lanes\" property in context panel")));
      await session.step(69, "Then the \"lane names\" reading of tile viewer should be \"F, M\"", () => readingReads(page, "lane names", el("tile viewer"), "F, M"));
      await session.step(70, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An explicit lane list overrides the column's category order", async () => {
      await session.step(73, "When user sets properties of tile viewer:", () => setProperties(page, el("tile viewer"), [["Lanes Column Name","RACE"],["Lanes","Black, Asian"]]), [["Lanes Column Name","RACE"],["Lanes","Black, Asian"]]);
      await session.step(76, "Then the \"lanes list\" reading of tile viewer should be \"Black, Asian\"", () => readingReads(page, "lanes list", el("tile viewer"), "Black, Asian"));
      await session.step(77, "And the \"lanes\" reading of tile viewer should be 2", () => readingIs(page, "lanes", el("tile viewer"), 2));
      await session.step(78, "And the \"lane names\" reading of tile viewer should be \"Black, Asian\"", () => readingReads(page, "lane names", el("tile viewer"), "Black, Asian"));
      await session.step(79, "And tile viewer should have a \"lane Black\" area", () => hasArea(page, el("tile viewer"), "lane Black"));
      await session.step(80, "And tile viewer should have a \"lane Asian\" area", () => hasArea(page, el("tile viewer"), "lane Asian"));
      await session.step(81, "And tile viewer should not have a \"lane Caucasian\" area", () => hasNoArea(page, el("tile viewer"), "lane Caucasian"));
      await session.step(82, "And tile viewer should not have a \"lane Other\" area", () => hasNoArea(page, el("tile viewer"), "lane Other"));
      await session.step(83, "And the \"lane of row 10\" reading of tile viewer should be \"Asian\"", () => readingReads(page, "lane of row 10", el("tile viewer"), "Asian"));
      await session.step(84, "And the \"lane of row 101\" reading of tile viewer should be \"Black\"", () => readingReads(page, "lane of row 101", el("tile viewer"), "Black"));
      await session.step(85, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A filter that empties a lane leaves the lane standing", async () => {
      await session.step(88, "Then the \"tiles in lane Asian\" reading of tile viewer should be at least 1", () => readingAtLeast(page, "tiles in lane Asian", el("tile viewer"), 1));
      await session.step(89, "When user filters out rows where \"RACE\" is \"Asian\"", () => filterOut(page, "RACE", "Asian"));
      await session.step(90, "Then 985 rows should pass the filter", () => filterPasses(page, 985));
      await session.step(91, "And tile viewer should show 985 rows", () => showsRows(page, el("tile viewer"), 985));
      await session.step(92, "And the \"lanes\" reading of tile viewer should be 2", () => readingIs(page, "lanes", el("tile viewer"), 2));
      await session.step(93, "And the \"lane names\" reading of tile viewer should be \"Black, Asian\"", () => readingReads(page, "lane names", el("tile viewer"), "Black, Asian"));
      await session.step(94, "And tile viewer should have a \"lane Asian\" area", () => hasArea(page, el("tile viewer"), "lane Asian"));
      await session.step(95, "And the \"tiles in lane Asian\" reading of tile viewer should be 0", () => readingIs(page, "tiles in lane Asian", el("tile viewer"), 0));
      await session.step(96, "And the \"tiles in lane Black\" reading of tile viewer should be at least 1", () => readingAtLeast(page, "tiles in lane Black", el("tile viewer"), 1));
      await session.step(97, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Empty Lanes off drops the emptied lane and on brings it back in its place", async () => {
      await session.step(100, "Then the \"lanes\" reading of tile viewer should be 2", () => readingIs(page, "lanes", el("tile viewer"), 2));
      await session.step(101, "And \"Show Empty Lanes\" property of tile viewer should be \"true\"", () => propertyShouldBe(page, "Show Empty Lanes", el("tile viewer"), "true"));
      await session.step(102, "When user sets \"Show Empty Lanes\" property of tile viewer to \"false\"", () => setProperty(page, "Show Empty Lanes", el("tile viewer"), "false"));
      await session.step(103, "Then the \"lanes\" reading of tile viewer should be 1", () => readingIs(page, "lanes", el("tile viewer"), 1));
      await session.step(104, "And the \"lane names\" reading of tile viewer should be \"Black\"", () => readingReads(page, "lane names", el("tile viewer"), "Black"));
      await session.step(105, "And tile viewer should not have a \"lane Asian\" area", () => hasNoArea(page, el("tile viewer"), "lane Asian"));
      await session.step(106, "And tile viewer should have a \"lane Black\" area", () => hasArea(page, el("tile viewer"), "lane Black"));
      await session.step(107, "When user sets \"Show Empty Lanes\" property of tile viewer to \"true\"", () => setProperty(page, "Show Empty Lanes", el("tile viewer"), "true"));
      await session.step(108, "Then the \"lanes\" reading of tile viewer should be 2", () => readingIs(page, "lanes", el("tile viewer"), 2));
      await session.step(109, "And the \"lane names\" reading of tile viewer should be \"Black, Asian\"", () => readingReads(page, "lane names", el("tile viewer"), "Black, Asian"));
      await session.step(110, "And tile viewer should have a \"lane Asian\" area", () => hasArea(page, el("tile viewer"), "lane Asian"));
      await session.step(111, "And the \"tiles in lane Asian\" reading of tile viewer should be 0", () => readingIs(page, "tiles in lane Asian", el("tile viewer"), 0));
      await session.step(112, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Clearing the lanes column leaves one lane holding the rows the filter left", async () => {
      await session.step(115, "When user sets \"Lanes Column Name\" property of tile viewer to \"\"", () => setProperty(page, "Lanes Column Name", el("tile viewer"), ""));
      await session.step(116, "Then the \"lanes\" reading of tile viewer should be 1", () => readingIs(page, "lanes", el("tile viewer"), 1));
      await session.step(117, "And the \"lane names\" reading of tile viewer should be \"All rows\"", () => readingReads(page, "lane names", el("tile viewer"), "All rows"));
      await session.step(118, "And the \"single lane\" reading of tile viewer should be \"true\"", () => readingReads(page, "single lane", el("tile viewer"), "true"));
      await session.step(119, "And tile viewer should show 985 rows", () => showsRows(page, el("tile viewer"), 985));
      await session.step(120, "And tile viewer should have a \"tile of row 1\" area", () => hasArea(page, el("tile viewer"), "tile of row 1"));
      await session.step(121, "And tile viewer should not have a \"tile of row 10\" area", () => hasNoArea(page, el("tile viewer"), "tile of row 10"));
      await session.step(122, "When user resets the filter", () => resetFilter(page));
      await session.step(123, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(124, "And tile viewer should show 1000 rows", () => showsRows(page, el("tile viewer"), 1000));
      await session.step(125, "When user sets properties of tile viewer:", () => setProperties(page, el("tile viewer"), [["Lanes",""],["Lanes Column Name","RACE"]]), [["Lanes",""],["Lanes Column Name","RACE"]]);
      await session.step(128, "Then the \"lanes\" reading of tile viewer should be 4", () => readingIs(page, "lanes", el("tile viewer"), 4));
      await session.step(129, "And the \"lanes list\" reading of tile viewer should be \"\"", () => readingReads(page, "lanes list", el("tile viewer"), ""));
      await session.step(130, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Dragging a card into another lane writes that lane's category into its row", async () => {
      await session.step(133, "Then the \"lane of row 10\" reading of tile viewer should be \"Asian\"", () => readingReads(page, "lane of row 10", el("tile viewer"), "Asian"));
      await session.step(134, "And the value of \"RACE\" column in row 10 should be \"Asian\"", () => valueInRow(page, "RACE", 10, "Asian"));
      await session.step(135, "When user moves the pointer away from tile viewer", () => pointerAway(page, el("tile viewer")));
      await session.step(136, "And user drags the card of row 10 of tile viewer into lane \"Black\"", () => dragCardIntoLane(page, 10, el("tile viewer"), "Black"));
      await session.step(137, "Then the value of \"RACE\" column in row 10 should be \"Black\"", () => valueInRow(page, "RACE", 10, "Black"));
      await session.step(138, "And the \"lane of row 10\" reading of tile viewer should be \"Black\"", () => readingReads(page, "lane of row 10", el("tile viewer"), "Black"));
      await session.step(139, "And the \"current row\" reading of tile viewer should be 10", () => readingIs(page, "current row", el("tile viewer"), 10));
      await session.step(140, "And the \"RACE of row 10\" reading of tile viewer should be \"Black\"", () => readingReads(page, "RACE of row 10", el("tile viewer"), "Black"));
      await session.step(141, "When user moves the pointer away from tile viewer", () => pointerAway(page, el("tile viewer")));
      await session.step(142, "And user drags the card of row 10 of tile viewer into lane \"Asian\"", () => dragCardIntoLane(page, 10, el("tile viewer"), "Asian"));
      await session.step(143, "Then the value of \"RACE\" column in row 10 should be \"Asian\"", () => valueInRow(page, "RACE", 10, "Asian"));
      await session.step(144, "And the \"lane of row 10\" reading of tile viewer should be \"Asian\"", () => readingReads(page, "lane of row 10", el("tile viewer"), "Asian"));
      await session.step(145, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A drop in the card's own lane changes nothing", async () => {
      await session.step(148, "Then the \"lane of row 101\" reading of tile viewer should be \"Black\"", () => readingReads(page, "lane of row 101", el("tile viewer"), "Black"));
      await session.step(149, "And the value of \"RACE\" column in row 101 should be \"Black\"", () => valueInRow(page, "RACE", 101, "Black"));
      await session.step(150, "When user moves the pointer away from tile viewer", () => pointerAway(page, el("tile viewer")));
      await session.step(151, "And user drags the card of row 101 of tile viewer into lane \"Black\"", () => dragCardIntoLane(page, 101, el("tile viewer"), "Black"));
      await session.step(152, "Then the value of \"RACE\" column in row 101 should be \"Black\"", () => valueInRow(page, "RACE", 101, "Black"));
      await session.step(153, "And the \"lane of row 101\" reading of tile viewer should be \"Black\"", () => readingReads(page, "lane of row 101", el("tile viewer"), "Black"));
      await session.step(154, "And the \"RACE of row 101\" reading of tile viewer should be \"Black\"", () => readingReads(page, "RACE of row 101", el("tile viewer"), "Black"));
      await session.step(155, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Allow Drag Between Lanes off blocks the write that the same drag makes when it is on", async () => {
      await session.step(158, "Then the \"drag between lanes\" reading of tile viewer should be \"true\"", () => readingReads(page, "drag between lanes", el("tile viewer"), "true"));
      await session.step(159, "And the \"lane of row 10\" reading of tile viewer should be \"Asian\"", () => readingReads(page, "lane of row 10", el("tile viewer"), "Asian"));
      await session.step(160, "When user sets \"Allow Drag Between Lanes\" property of tile viewer to \"false\"", () => setProperty(page, "Allow Drag Between Lanes", el("tile viewer"), "false"));
      await session.step(161, "Then the \"drag between lanes\" reading of tile viewer should be \"false\"", () => readingReads(page, "drag between lanes", el("tile viewer"), "false"));
      await session.step(162, "When user moves the pointer away from tile viewer", () => pointerAway(page, el("tile viewer")));
      await session.step(163, "And user drags the card of row 10 of tile viewer into lane \"Black\"", () => dragCardIntoLane(page, 10, el("tile viewer"), "Black"));
      await session.step(164, "Then the value of \"RACE\" column in row 10 should be \"Asian\"", () => valueInRow(page, "RACE", 10, "Asian"));
      await session.step(165, "And the \"lane of row 10\" reading of tile viewer should be \"Asian\"", () => readingReads(page, "lane of row 10", el("tile viewer"), "Asian"));
      await session.step(166, "When user sets \"Allow Drag Between Lanes\" property of tile viewer to \"true\"", () => setProperty(page, "Allow Drag Between Lanes", el("tile viewer"), "true"));
      await session.step(167, "Then the \"drag between lanes\" reading of tile viewer should be \"true\"", () => readingReads(page, "drag between lanes", el("tile viewer"), "true"));
      await session.step(168, "When user moves the pointer away from tile viewer", () => pointerAway(page, el("tile viewer")));
      await session.step(169, "And user drags the card of row 10 of tile viewer into lane \"Black\"", () => dragCardIntoLane(page, 10, el("tile viewer"), "Black"));
      await session.step(170, "Then the value of \"RACE\" column in row 10 should be \"Black\"", () => valueInRow(page, "RACE", 10, "Black"));
      await session.step(171, "And the \"lane of row 10\" reading of tile viewer should be \"Black\"", () => readingReads(page, "lane of row 10", el("tile viewer"), "Black"));
      await session.step(172, "When user moves the pointer away from tile viewer", () => pointerAway(page, el("tile viewer")));
      await session.step(173, "And user drags the card of row 10 of tile viewer into lane \"Asian\"", () => dragCardIntoLane(page, 10, el("tile viewer"), "Asian"));
      await session.step(174, "Then the value of \"RACE\" column in row 10 should be \"Asian\"", () => valueInRow(page, "RACE", 10, "Asian"));
      await session.step(175, "And the \"lane of row 10\" reading of tile viewer should be \"Asian\"", () => readingReads(page, "lane of row 10", el("tile viewer"), "Asian"));
      await session.step(176, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The viewer's menu is reached where no card is", async () => {
      await session.step(179, "Then tile viewer should have a \"viewer menu\" area", () => hasArea(page, el("tile viewer"), "viewer menu"));
      await session.step(180, "When user opens the viewer menu of tile viewer", () => openViewerMenu(page, el("tile viewer")));
      await session.step(181, "Then the open menu should list \"Edit Form...\"", () => menuLists(page, "Edit Form..."));
      await session.step(182, "And the open menu should list \"Lanes\"", () => menuLists(page, "Lanes"));
      await session.step(183, "And the open menu should list \"Show Empty Lanes\"", () => menuLists(page, "Show Empty Lanes"));
      await session.step(184, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(185, "And user picks \"Show Empty Lanes\" from the viewer menu of tile viewer", () => pickFromViewerMenu(page, "Show Empty Lanes", el("tile viewer")));
      await session.step(186, "Then \"Show Empty Lanes\" property of tile viewer should be \"false\"", () => propertyShouldBe(page, "Show Empty Lanes", el("tile viewer"), "false"));
      await session.step(187, "When user picks \"Show Empty Lanes\" from the viewer menu of tile viewer", () => pickFromViewerMenu(page, "Show Empty Lanes", el("tile viewer")));
      await session.step(188, "Then \"Show Empty Lanes\" property of tile viewer should be \"true\"", () => propertyShouldBe(page, "Show Empty Lanes", el("tile viewer"), "true"));
      await session.step(189, "And the \"lanes\" reading of tile viewer should be 4", () => readingIs(page, "lanes", el("tile viewer"), 4));
      await session.step(190, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
