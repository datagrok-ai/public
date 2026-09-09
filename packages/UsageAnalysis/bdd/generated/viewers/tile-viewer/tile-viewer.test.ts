/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tile-viewer/tile-viewer.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tile-viewer]
--- */
import {test} from '@playwright/test';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {addTileViewer, descriptionAbove, descriptionBelow, everyTileBetween, everyTileShows, readingContains, readingNotContains} from '../../../bindings/tile-viewer.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnCount} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {clearSelection, filterPasses, filterPassesAll, filterTo, noneSelected, resetFilter, selectWhereIs, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {areaShorter, areaTaller, boundTable, closeContextMenu, hasArea, hasNoArea, menuDoesNotList, menuLists, noErrors, propertyShouldBe, readingAtLeast, readingDoesNotRead, readingIs, readingReads, readingsEqual, rightClickArea, setProperties, setProperty, showsFewerRows, showsMoreRows, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Tile viewer property surface", () => {
  const session = feature(test, "features/viewers/tile-viewer/tile-viewer.feature", import.meta.url);
  test("Tile viewer property surface", {tag: ["@journey", "@viewers", "@realizes:viewers.tile-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(31, "And user adds a tile viewer", () => addTileViewer(page));
    await session.step(32, "Then tile viewer should be visible", () => shouldBe(page, el("tile viewer"), "visible"));
    await session.step(33, "And the \"lanes\" reading of tile viewer should be 1", () => readingIs(page, "lanes", el("tile viewer"), 1));
    await session.step(34, "And the \"single lane\" reading of tile viewer should be \"true\"", () => readingReads(page, "single lane", el("tile viewer"), "true"));
    await session.step(35, "And tile viewer should show 1000 rows", () => showsRows(page, el("tile viewer"), 1000));
    await session.step(36, "And the \"fields shown\" reading of tile viewer should be 10", () => readingIs(page, "fields shown", el("tile viewer"), 10));
    await run.scenario("A card per row, showing the grid's display string for every field it holds", async () => {
      await session.step(39, "Then the table should have 11 columns", () => columnCount(page, 11));
      await session.step(40, "And the \"table\" reading of tile viewer should be \"demog-1000\"", () => readingReads(page, "table", el("tile viewer"), "demog-1000"));
      await session.step(41, "And the \"auto generate\" reading of tile viewer should be \"true\"", () => readingReads(page, "auto generate", el("tile viewer"), "true"));
      await session.step(42, "And the \"form designed\" reading of tile viewer should be \"false\"", () => readingReads(page, "form designed", el("tile viewer"), "false"));
      await session.step(43, "And the \"fields\" reading of tile viewer should contain \"AGE\"", () => readingContains(page, "fields", el("tile viewer"), "AGE"));
      await session.step(44, "And the \"fields\" reading of tile viewer should contain \"USUBJID\"", () => readingContains(page, "fields", el("tile viewer"), "USUBJID"));
      await session.step(45, "And the \"fields\" reading of tile viewer should contain \"WEIGHT\"", () => readingContains(page, "fields", el("tile viewer"), "WEIGHT"));
      await session.step(46, "And the \"fields\" reading of tile viewer should not contain \"SEVERITY\"", () => readingNotContains(page, "fields", el("tile viewer"), "SEVERITY"));
      await session.step(47, "And the \"tiles\" reading of tile viewer should be at least 3", () => readingAtLeast(page, "tiles", el("tile viewer"), 3));
      await session.step(48, "And the \"tiles\" and \"tiles in lane All rows\" readings of tile viewer should be the same", () => readingsEqual(page, "tiles", "tiles in lane All rows", el("tile viewer")));
      await session.step(49, "And tile viewer should have a \"tile of row 1\" area", () => hasArea(page, el("tile viewer"), "tile of row 1"));
      await session.step(50, "And tile viewer should have a \"field AGE of row 1\" area", () => hasArea(page, el("tile viewer"), "field AGE of row 1"));
      await session.step(51, "And tile viewer should have a \"label AGE of row 1\" area", () => hasArea(page, el("tile viewer"), "label AGE of row 1"));
      await session.step(52, "And tile viewer should not have a \"field SEVERITY of row 1\" area", () => hasNoArea(page, el("tile viewer"), "field SEVERITY of row 1"));
      await session.step(53, "And the \"current row\" reading of tile viewer should be 1", () => readingIs(page, "current row", el("tile viewer"), 1));
      await session.step(54, "And the \"lane of row 1\" reading of tile viewer should be \"All rows\"", () => readingReads(page, "lane of row 1", el("tile viewer"), "All rows"));
      await session.step(55, "And the \"USUBJID of row 1\" reading of tile viewer should be \"X0273T21000300003\"", () => readingReads(page, "USUBJID of row 1", el("tile viewer"), "X0273T21000300003"));
      await session.step(56, "And the \"AGE of row 1\" reading of tile viewer should be \"26\"", () => readingReads(page, "AGE of row 1", el("tile viewer"), "26"));
      await session.step(57, "And the \"SEX of row 1\" reading of tile viewer should be \"F\"", () => readingReads(page, "SEX of row 1", el("tile viewer"), "F"));
      await session.step(58, "And the \"RACE of row 1\" reading of tile viewer should be \"Caucasian\"", () => readingReads(page, "RACE of row 1", el("tile viewer"), "Caucasian"));
      await session.step(59, "And the \"DIS_POP of row 1\" reading of tile viewer should be \"Indigestion\"", () => readingReads(page, "DIS_POP of row 1", el("tile viewer"), "Indigestion"));
      await session.step(60, "And the \"AGE of row 2\" reading of tile viewer should be \"30\"", () => readingReads(page, "AGE of row 2", el("tile viewer"), "30"));
      await session.step(61, "And the \"RACE of row 2\" reading of tile viewer should be \"Other\"", () => readingReads(page, "RACE of row 2", el("tile viewer"), "Other"));
      await session.step(62, "And the \"WEIGHT of row 1\" reading of tile viewer should be \"74.10\"", () => readingReads(page, "WEIGHT of row 1", el("tile viewer"), "74.10"));
      await session.step(63, "And the \"WEIGHT of row 1\" reading of tile viewer should not be \"74.1\"", () => readingDoesNotRead(page, "WEIGHT of row 1", el("tile viewer"), "74.1"));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A right-click on a card opens the column's menu, and the viewer has no menu region here", async () => {
      await session.step(67, "Then tile viewer should not have a \"viewer menu\" area", () => hasNoArea(page, el("tile viewer"), "viewer menu"));
      await session.step(68, "When user right-clicks on the \"field AGE of row 1\" area of tile viewer", () => rightClickArea(page, "field AGE of row 1", el("tile viewer")));
      await session.step(69, "Then the open menu should list \"Remove\"", () => menuLists(page, "Remove"));
      await session.step(70, "And the open menu should list \"Rename...\"", () => menuLists(page, "Rename..."));
      await session.step(71, "And the open menu should not list \"Edit Form...\"", () => menuDoesNotList(page, "Edit Form..."));
      await session.step(72, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(73, "And user right-clicks on the \"label AGE of row 1\" area of tile viewer", () => rightClickArea(page, "label AGE of row 1", el("tile viewer")));
      await session.step(74, "Then the open menu should list \"Rename...\"", () => menuLists(page, "Rename..."));
      await session.step(75, "And the open menu should not list \"Show Empty Lanes\"", () => menuDoesNotList(page, "Show Empty Lanes"));
      await session.step(76, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(77, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Tiles Font sizes the lane headers", async () => {
      await session.step(80, "Given user sets \"Lanes Column Name\" property of tile viewer to \"SEX\"", () => setProperty(page, "Lanes Column Name", el("tile viewer"), "SEX"));
      await session.step(81, "Then the \"lane names\" reading of tile viewer should be \"F, M\"", () => readingReads(page, "lane names", el("tile viewer"), "F, M"));
      await session.step(82, "And tile viewer should have a \"lane header F\" area", () => hasArea(page, el("tile viewer"), "lane header F"));
      await session.step(83, "And the \"tiles font\" reading of tile viewer should be 'normal normal 13px \"Roboto\"'", () => readingReads(page, "tiles font", el("tile viewer"), "normal normal 13px \"Roboto\""));
      await session.step(84, "When user sets \"Tiles Font\" property of tile viewer to 'normal normal 18px \"Roboto\"'", () => setProperty(page, "Tiles Font", el("tile viewer"), "normal normal 18px \"Roboto\""));
      await session.step(85, "Then the \"tiles font\" reading of tile viewer should be 'normal normal 18px \"Roboto\"'", () => readingReads(page, "tiles font", el("tile viewer"), "normal normal 18px \"Roboto\""));
      await session.step(86, "And the \"lane header F\" area of tile viewer should be taller than before", () => areaTaller(page, "lane header F", el("tile viewer")));
      await session.step(87, "When user sets \"Tiles Font\" property of tile viewer to 'normal normal 13px \"Roboto\"'", () => setProperty(page, "Tiles Font", el("tile viewer"), "normal normal 13px \"Roboto\""));
      await session.step(88, "Then the \"tiles font\" reading of tile viewer should be 'normal normal 13px \"Roboto\"'", () => readingReads(page, "tiles font", el("tile viewer"), "normal normal 13px \"Roboto\""));
      await session.step(89, "And the \"lane header F\" area of tile viewer should be shorter than before", () => areaShorter(page, "lane header F", el("tile viewer")));
      await session.step(90, "When user sets \"Lanes Column Name\" property of tile viewer to \"\"", () => setProperty(page, "Lanes Column Name", el("tile viewer"), ""));
      await session.step(91, "Then the \"single lane\" reading of tile viewer should be \"true\"", () => readingReads(page, "single lane", el("tile viewer"), "true"));
      await session.step(92, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Title and description", async () => {
      await session.step(95, "When user sets properties of tile viewer:", () => setProperties(page, el("tile viewer"), [["Show Title","true"],["Title","Cards"]]));
      await session.step(98, "Then title of tile viewer should be visible", () => shouldBe(page, el("title of tile viewer"), "visible"));
      await session.step(99, "And title of tile viewer should have text \"Cards\"", () => shouldHaveText(page, el("title of tile viewer"), "Cards"));
      await session.step(100, "When user sets properties of tile viewer:", () => setProperties(page, el("tile viewer"), [["Description","A card per patient"],["Description Visibility Mode","Always"]]));
      await session.step(103, "Then description of tile viewer should have text \"A card per patient\"", () => shouldHaveText(page, el("description of tile viewer"), "A card per patient"));
      await session.step(104, "And the description of tile viewer should be above its content", () => descriptionAbove(page, el("tile viewer")));
      await session.step(105, "When user sets \"Description Position\" property of tile viewer to \"Bottom\"", () => setProperty(page, "Description Position", el("tile viewer"), "Bottom"));
      await session.step(106, "Then description of tile viewer should be visible", () => shouldBe(page, el("description of tile viewer"), "visible"));
      await session.step(107, "And the description of tile viewer should be below its content", () => descriptionBelow(page, el("tile viewer")));
      await session.step(108, "When user sets \"Description Position\" property of tile viewer to \"Top\"", () => setProperty(page, "Description Position", el("tile viewer"), "Top"));
      await session.step(109, "Then the description of tile viewer should be above its content", () => descriptionAbove(page, el("tile viewer")));
      await session.step(110, "When user sets \"Description Visibility Mode\" property of tile viewer to \"Never\"", () => setProperty(page, "Description Visibility Mode", el("tile viewer"), "Never"));
      await session.step(111, "Then description of tile viewer should be absent", () => shouldBe(page, el("description of tile viewer"), "absent"));
      await session.step(112, "When user sets \"Title\" property of tile viewer to \"\"", () => setProperty(page, "Title", el("tile viewer"), ""));
      await session.step(113, "Then title of tile viewer should be hidden", () => shouldBe(page, el("title of tile viewer"), "hidden"));
      await session.step(114, "When user sets properties of tile viewer:", () => setProperties(page, el("tile viewer"), [["Show Title","false"],["Description",""],["Description Visibility Mode","Auto"]]));
      await session.step(118, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Row Source picks the rows the cards are built from", async () => {
      await session.step(121, "Then \"Row Source\" property of tile viewer should be \"Filtered\"", () => propertyShouldBe(page, "Row Source", el("tile viewer"), "Filtered"));
      await session.step(122, "And tile viewer should show 1000 rows", () => showsRows(page, el("tile viewer"), 1000));
      await session.step(123, "When user selects rows where \"RACE\" is \"Asian\"", () => selectWhereIs(page, "RACE", "Asian"));
      await session.step(124, "Then 15 rows should be selected", () => selectedRowCount(page, 15));
      await session.step(125, "When user sets \"Row Source\" property of tile viewer to \"Selected\"", () => setProperty(page, "Row Source", el("tile viewer"), "Selected"));
      await session.step(126, "Then tile viewer should show 15 rows", () => showsRows(page, el("tile viewer"), 15));
      await session.step(127, "And tile viewer should show fewer rows than before", () => showsFewerRows(page, el("tile viewer")));
      await session.step(128, "And every tile of tile viewer should show \"Asian\" in \"RACE\"", () => everyTileShows(page, el("tile viewer"), "Asian", "RACE"));
      await session.step(129, "And the \"selected rows shown\" reading of tile viewer should be \"false\"", () => readingReads(page, "selected rows shown", el("tile viewer"), "false"));
      await session.step(130, "When user sets \"Row Source\" property of tile viewer to \"All\"", () => setProperty(page, "Row Source", el("tile viewer"), "All"));
      await session.step(131, "Then tile viewer should show 1000 rows", () => showsRows(page, el("tile viewer"), 1000));
      await session.step(132, "And tile viewer should show more rows than before", () => showsMoreRows(page, el("tile viewer")));
      await session.step(133, "And the \"selected rows shown\" reading of tile viewer should be \"true\"", () => readingReads(page, "selected rows shown", el("tile viewer"), "true"));
      await session.step(134, "And the \"RACE of row 1\" reading of tile viewer should be \"Caucasian\"", () => readingReads(page, "RACE of row 1", el("tile viewer"), "Caucasian"));
      await session.step(135, "When user sets \"Row Source\" property of tile viewer to \"Filtered\"", () => setProperty(page, "Row Source", el("tile viewer"), "Filtered"));
      await session.step(136, "And user clears the row selection", () => clearSelection(page));
      await session.step(137, "Then tile viewer should show 1000 rows", () => showsRows(page, el("tile viewer"), 1000));
      await session.step(138, "And no rows should be selected", () => noneSelected(page));
      await session.step(139, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The viewer's own Filter formula narrows the cards and leaves the table's filter alone", async () => {
      await session.step(142, "Then tile viewer should show 1000 rows", () => showsRows(page, el("tile viewer"), 1000));
      await session.step(143, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(144, "When user sets \"Filter\" property of tile viewer to \"${AGE} > 50\"", () => setProperty(page, "Filter", el("tile viewer"), "${AGE} > 50"));
      await session.step(145, "Then tile viewer should show 367 rows", () => showsRows(page, el("tile viewer"), 367));
      await session.step(146, "And tile viewer should show fewer rows than before", () => showsFewerRows(page, el("tile viewer")));
      await session.step(147, "And every tile of tile viewer should show a value between 51 and 89 in \"AGE\"", () => everyTileBetween(page, el("tile viewer"), 51, 89, "AGE"));
      await session.step(148, "And 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(149, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(150, "When user sets \"Filter\" property of tile viewer to \"\"", () => setProperty(page, "Filter", el("tile viewer"), ""));
      await session.step(151, "Then tile viewer should show 1000 rows", () => showsRows(page, el("tile viewer"), 1000));
      await session.step(152, "And tile viewer should show more rows than before", () => showsMoreRows(page, el("tile viewer")));
      await session.step(153, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A filter on the table reaches the cards", async () => {
      await session.step(156, "Then tile viewer should show 1000 rows", () => showsRows(page, el("tile viewer"), 1000));
      await session.step(157, "When user filters rows where \"SEX\" is \"M\"", () => filterTo(page, "SEX", "M"));
      await session.step(158, "Then 447 rows should pass the filter", () => filterPasses(page, 447));
      await session.step(159, "And tile viewer should show 447 rows", () => showsRows(page, el("tile viewer"), 447));
      await session.step(160, "And tile viewer should show fewer rows than before", () => showsFewerRows(page, el("tile viewer")));
      await session.step(161, "And every tile of tile viewer should show \"M\" in \"SEX\"", () => everyTileShows(page, el("tile viewer"), "M", "SEX"));
      await session.step(162, "When user resets the filter", () => resetFilter(page));
      await session.step(163, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(164, "And tile viewer should show 1000 rows", () => showsRows(page, el("tile viewer"), 1000));
      await session.step(165, "And tile viewer should show more rows than before", () => showsMoreRows(page, el("tile viewer")));
      await session.step(166, "And the \"SEX of row 1\" reading of tile viewer should be \"F\"", () => readingReads(page, "SEX of row 1", el("tile viewer"), "F"));
      await session.step(167, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Table property rebinds the cards to another table", async () => {
      await session.step(170, "Given user opens spgi dataset", () => openDataset(page, ds("spgi")));
      await session.step(171, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
      await session.step(172, "When user sets \"Table\" property of tile viewer to \"spgi-100\"", () => setProperty(page, "Table", el("tile viewer"), "spgi-100"));
      await session.step(173, "Then tile viewer should be bound to table \"spgi-100\"", () => boundTable(page, el("tile viewer"), "spgi-100"));
      await session.step(174, "And the \"table\" reading of tile viewer should be \"spgi-100\"", () => readingReads(page, "table", el("tile viewer"), "spgi-100"));
      await session.step(175, "And tile viewer should show 100 rows", () => showsRows(page, el("tile viewer"), 100));
      await session.step(176, "And the \"fields\" reading of tile viewer should contain \"Primary Series Name\"", () => readingContains(page, "fields", el("tile viewer"), "Primary Series Name"));
      await session.step(177, "And the \"fields\" reading of tile viewer should not contain \"SEX\"", () => readingNotContains(page, "fields", el("tile viewer"), "SEX"));
      await session.step(178, "And tile viewer should not have a \"field SEX of row 1\" area", () => hasNoArea(page, el("tile viewer"), "field SEX of row 1"));
      await session.step(179, "And tile viewer should have a \"tile of row 1\" area", () => hasArea(page, el("tile viewer"), "tile of row 1"));
      await session.step(180, "When user sets \"Table\" property of tile viewer to \"demog-1000\"", () => setProperty(page, "Table", el("tile viewer"), "demog-1000"));
      await session.step(181, "Then tile viewer should be bound to table \"demog-1000\"", () => boundTable(page, el("tile viewer"), "demog-1000"));
      await session.step(182, "And tile viewer should show 1000 rows", () => showsRows(page, el("tile viewer"), 1000));
      await session.step(183, "And the \"fields shown\" reading of tile viewer should be 10", () => readingIs(page, "fields shown", el("tile viewer"), 10));
      await session.step(184, "And the \"fields\" reading of tile viewer should contain \"SEX\"", () => readingContains(page, "fields", el("tile viewer"), "SEX"));
      await session.step(185, "And the \"SEX of row 1\" reading of tile viewer should be \"F\"", () => readingReads(page, "SEX of row 1", el("tile viewer"), "F"));
      await session.step(186, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
