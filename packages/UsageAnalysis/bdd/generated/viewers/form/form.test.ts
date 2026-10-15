/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/form/form.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.form]
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
import {clickOn, downloadThrough, downloadedContains, downloadedNotContains, focusOn, pressKeyIn, shouldBe, shouldNotBe, uploadDownloaded} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs, hasColumn, hasNoColumn, makeRowCurrent, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {addCalculated, clearSelection, colorLinear, colorOff, filterTo, noneSelected, onlyOfSelected, removeColumn, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, dragAreaBy, enterIntoArea, hasArea, hasNoArea, hoverArea, loadLayout, noErrors, pickFromContextMenu, pointerAway, propertyShouldBe, readingDiffers, readingIs, readingReads, saveLayout, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {areaHangsBelow} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Form viewer", () => {
  const session = feature(test, "features/viewers/form/form.feature", import.meta.url);
  test("Form viewer", {tag: ["@journey", "@viewers", "@realizes:viewers.form"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 16, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user adds a form viewer", () => addViewer(page, "form"));
    await session.step(20, "Then form viewer should be visible", () => shouldBe(page, el("form viewer"), "visible"));
    await session.step(21, "And the \"fields shown\" reading of form viewer should be 11", () => readingIs(page, "fields shown", el("form viewer"), 11));
    await session.step(22, "And form viewer should have a \"field USUBJID\" area", () => hasArea(page, el("form viewer"), "field USUBJID"));
    await session.step(23, "And form viewer should have a \"field SEVERITY\" area", () => hasArea(page, el("form viewer"), "field SEVERITY"));
    await run.scenario("The fields show the current row", async () => {
      await session.step(26, "Then the \"row\" reading of form viewer should be 1", () => readingIs(page, "row", el("form viewer"), 1));
      await session.step(27, "And the \"USUBJID\" reading of form viewer should be \"X0273T21000300003\"", () => readingReads(page, "USUBJID", el("form viewer"), "X0273T21000300003"));
      await session.step(28, "And the \"AGE\" reading of form viewer should be \"26\"", () => readingReads(page, "AGE", el("form viewer"), "26"));
      await session.step(29, "And the \"SEX\" reading of form viewer should be \"F\"", () => readingReads(page, "SEX", el("form viewer"), "F"));
      await session.step(30, "And the \"RACE\" reading of form viewer should be \"Caucasian\"", () => readingReads(page, "RACE", el("form viewer"), "Caucasian"));
      await session.step(31, "And the \"DIS_POP\" reading of form viewer should be \"Indigestion\"", () => readingReads(page, "DIS_POP", el("form viewer"), "Indigestion"));
      await session.step(32, "And form viewer should have a \"field AGE\" area", () => hasArea(page, el("form viewer"), "field AGE"));
      await session.step(33, "And form viewer should have a \"label AGE\" area", () => hasArea(page, el("form viewer"), "label AGE"));
      await session.step(34, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column's colour coding reaches its field", async () => {
      await session.step(37, "Then the \"background of AGE\" reading of form viewer should be \"\"", () => readingReads(page, "background of AGE", el("form viewer"), ""));
      await session.step(38, "When user colors \"AGE\" column linearly from \"#FF0000\" to \"#0000FF\"", () => colorLinear(page, "AGE", "#FF0000", "#0000FF"));
      await session.step(39, "Then the \"background of AGE\" reading of form viewer should differ from before", () => readingDiffers(page, "background of AGE", el("form viewer")));
      await session.step(40, "When user removes the coloring of \"AGE\" column", () => colorOff(page, "AGE"));
      await session.step(41, "Then the \"background of AGE\" reading of form viewer should be \"\"", () => readingReads(page, "background of AGE", el("form viewer"), ""));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The arrows walk the rows", async () => {
      await session.step(45, "When user clicks on the \"next row\" area of form viewer", () => clickArea(page, "next row", el("form viewer")));
      await session.step(46, "Then row 2 should be current", () => currentRowIs(page, 2));
      await session.step(47, "And the \"USUBJID\" reading of form viewer should be \"X0273T21000300005\"", () => readingReads(page, "USUBJID", el("form viewer"), "X0273T21000300005"));
      await session.step(48, "And the \"AGE\" reading of form viewer should be \"30\"", () => readingReads(page, "AGE", el("form viewer"), "30"));
      await session.step(49, "When user clicks on the \"previous row\" area of form viewer", () => clickArea(page, "previous row", el("form viewer")));
      await session.step(50, "Then row 1 should be current", () => currentRowIs(page, 1));
      await session.step(51, "And the \"USUBJID\" reading of form viewer should be \"X0273T21000300003\"", () => readingReads(page, "USUBJID", el("form viewer"), "X0273T21000300003"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The form follows the grid's current row", async () => {
      await session.step(55, "When user makes row 43 current", () => makeRowCurrent(page, 43));
      await session.step(56, "Then the \"row\" reading of form viewer should be 43", () => readingIs(page, "row", el("form viewer"), 43));
      await session.step(57, "And the \"USUBJID\" reading of form viewer should be \"X0273T21002400007\"", () => readingReads(page, "USUBJID", el("form viewer"), "X0273T21002400007"));
      await session.step(58, "And the \"AGE\" reading of form viewer should be \"54\"", () => readingReads(page, "AGE", el("form viewer"), "54"));
      await session.step(59, "When user makes row 1 current", () => makeRowCurrent(page, 1));
      await session.step(60, "Then the \"row\" reading of form viewer should be 1", () => readingIs(page, "row", el("form viewer"), 1));
    });
    await run.scenario("The row selector toggles the row's selection", async () => {
      await session.step(63, "Given user clears the row selection", () => clearSelection(page));
      await session.step(64, "When user clicks on the \"row selector\" area of form viewer", () => clickArea(page, "row selector", el("form viewer")));
      await session.step(65, "Then only rows where \"USUBJID\" is \"X0273T21000300003\" should be selected", () => onlyOfSelected(page, "USUBJID", "X0273T21000300003"));
      await session.step(66, "And the \"row selected\" reading of form viewer should be \"true\"", () => readingReads(page, "row selected", el("form viewer"), "true"));
      await session.step(67, "And \"square\" icon in form viewer should be selected", () => shouldBe(page, el("\"square\" icon in form viewer"), "selected"));
      await session.step(68, "When user clicks on the \"row selector\" area of form viewer", () => clickArea(page, "row selector", el("form viewer")));
      await session.step(69, "Then no rows should be selected", () => noneSelected(page));
      await session.step(70, "And the \"row selected\" reading of form viewer should be \"false\"", () => readingReads(page, "row selected", el("form viewer"), "false"));
      await session.step(71, "And \"square\" icon in form viewer should not be selected", () => shouldNotBe(page, el("\"square\" icon in form viewer"), "selected"));
    });
    await run.scenario("Track Row: the form follows the mouse-over row, or nothing", async () => {
      await session.step(74, "When user picks \"Track Row > Mouse Over\" from the context menu of form viewer", () => pickFromContextMenu(page, "Track Row > Mouse Over", el("form viewer")));
      await session.step(75, "Then \"Sync Mode\" property of form viewer should be \"Mouse Over\"", () => propertyShouldBe(page, "Sync Mode", el("form viewer"), "Mouse Over"));
      await session.step(76, "When user hovers over the \"cell 2 of AGE\" area of grid", () => hoverArea(page, "cell 2 of AGE", el("grid")));
      await session.step(77, "Then the \"USUBJID\" reading of form viewer should be \"X0273T21000300005\"", () => readingReads(page, "USUBJID", el("form viewer"), "X0273T21000300005"));
      await session.step(78, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(79, "And user picks \"Track Row > None\" from the context menu of form viewer", () => pickFromContextMenu(page, "Track Row > None", el("form viewer")));
      await session.step(80, "Then \"Sync Mode\" property of form viewer should be \"None\"", () => propertyShouldBe(page, "Sync Mode", el("form viewer"), "None"));
      await session.step(81, "When user makes row 5 current", () => makeRowCurrent(page, 5));
      await session.step(82, "Then the \"USUBJID\" reading of form viewer should be \"X0273T21000300005\"", () => readingReads(page, "USUBJID", el("form viewer"), "X0273T21000300005"));
      await session.step(83, "When user picks \"Track Row > Current\" from the context menu of form viewer", () => pickFromContextMenu(page, "Track Row > Current", el("form viewer")));
      await session.step(84, "Then \"Sync Mode\" property of form viewer should be \"Current\"", () => propertyShouldBe(page, "Sync Mode", el("form viewer"), "Current"));
      await session.step(85, "And the \"USUBJID\" reading of form viewer should be \"X0273T21000500006\"", () => readingReads(page, "USUBJID", el("form viewer"), "X0273T21000500006"));
      await session.step(86, "When user makes row 1 current", () => makeRowCurrent(page, 1));
    });
    await run.scenario("The keyboard walks and selects rows", async () => {
      await session.step(89, "Given user clears the row selection", () => clearSelection(page));
      await session.step(90, "When user focuses on form viewer", () => focusOn(page, el("form viewer")));
      await session.step(91, "And user presses ArrowRight in form viewer", () => pressKeyIn(page, "ArrowRight", el("form viewer")));
      await session.step(92, "Then row 2 should be current", () => currentRowIs(page, 2));
      await session.step(93, "When user presses ArrowDown in form viewer", () => pressKeyIn(page, "ArrowDown", el("form viewer")));
      await session.step(94, "Then row 3 should be current", () => currentRowIs(page, 3));
      await session.step(95, "When user presses ArrowLeft in form viewer", () => pressKeyIn(page, "ArrowLeft", el("form viewer")));
      await session.step(96, "Then row 2 should be current", () => currentRowIs(page, 2));
      await session.step(97, "When user presses ArrowUp in form viewer", () => pressKeyIn(page, "ArrowUp", el("form viewer")));
      await session.step(98, "Then row 1 should be current", () => currentRowIs(page, 1));
      await session.step(99, "When user presses Space in form viewer", () => pressKeyIn(page, "Space", el("form viewer")));
      await session.step(100, "Then only rows where \"USUBJID\" is \"X0273T21000300003\" should be selected", () => onlyOfSelected(page, "USUBJID", "X0273T21000300003"));
      await session.step(101, "When user clears the row selection", () => clearSelection(page));
    });
    await run.scenario("A field writes its cell only while the form is editable", async () => {
      await session.step(104, "Then the \"editable\" reading of form viewer should be \"false\"", () => readingReads(page, "editable", el("form viewer"), "false"));
      await session.step(105, "When user clicks on the \"edit\" area of form viewer", () => clickArea(page, "edit", el("form viewer")));
      await session.step(106, "Then the \"editable\" reading of form viewer should be \"true\"", () => readingReads(page, "editable", el("form viewer"), "true"));
      await session.step(107, "When user enters \"31\" into the \"field AGE\" area of form viewer", () => enterIntoArea(page, "31", "field AGE", el("form viewer")));
      await session.step(108, "Then the value of \"AGE\" column in row 1 should be \"31\"", () => valueInRow(page, "AGE", 1, "31"));
      await session.step(109, "And the \"AGE\" reading of form viewer should be \"31\"", () => readingReads(page, "AGE", el("form viewer"), "31"));
      await session.step(110, "When user enters \"26\" into the \"field AGE\" area of form viewer", () => enterIntoArea(page, "26", "field AGE", el("form viewer")));
      await session.step(111, "And user clicks on the \"edit\" area of form viewer", () => clickArea(page, "edit", el("form viewer")));
      await session.step(112, "Then the \"editable\" reading of form viewer should be \"false\"", () => readingReads(page, "editable", el("form viewer"), "false"));
      await session.step(113, "And the value of \"AGE\" column in row 1 should be \"26\"", () => valueInRow(page, "AGE", 1, "26"));
      await session.step(114, "When user enters \"99\" into the \"field AGE\" area of form viewer", () => enterIntoArea(page, "99", "field AGE", el("form viewer")));
      await session.step(115, "Then the value of \"AGE\" column in row 1 should be \"26\"", () => valueInRow(page, "AGE", 1, "26"));
    });
    await run.scenario("The toolbar shows what the properties allow", async () => {
      await session.step(118, "When user sets \"Show Next Row Arrow\" property of form viewer to \"false\"", () => setProperty(page, "Show Next Row Arrow", el("form viewer"), "false"));
      await session.step(119, "Then form viewer should not have a \"next row\" area", () => hasNoArea(page, el("form viewer"), "next row"));
      await session.step(120, "And form viewer should have a \"previous row\" area", () => hasArea(page, el("form viewer"), "previous row"));
      await session.step(121, "When user sets \"Show Row Selector\" property of form viewer to \"false\"", () => setProperty(page, "Show Row Selector", el("form viewer"), "false"));
      await session.step(122, "Then form viewer should not have a \"row selector\" area", () => hasNoArea(page, el("form viewer"), "row selector"));
      await session.step(123, "When user sets properties of form viewer:", () => setProperties(page, el("form viewer"), [["Show Next Row Arrow","true"],["Show Row Selector","true"]]), [["Show Next Row Arrow","true"],["Show Row Selector","true"]]);
      await session.step(126, "Then form viewer should have a \"next row\" area", () => hasArea(page, el("form viewer"), "next row"));
      await session.step(127, "And form viewer should have a \"row selector\" area", () => hasArea(page, el("form viewer"), "row selector"));
      await session.step(128, "When user sets \"Show Navigation\" property of form viewer to \"false\"", () => setProperty(page, "Show Navigation", el("form viewer"), "false"));
      await session.step(129, "Then form viewer should not have a \"next row\" area", () => hasNoArea(page, el("form viewer"), "next row"));
      await session.step(130, "And form viewer should not have a \"select columns\" area", () => hasNoArea(page, el("form viewer"), "select columns"));
      await session.step(131, "When user sets \"Show Navigation\" property of form viewer to \"true\"", () => setProperty(page, "Show Navigation", el("form viewer"), "true"));
      await session.step(132, "Then form viewer should have a \"next row\" area", () => hasArea(page, el("form viewer"), "next row"));
    });
    await run.scenario("The field set follows the columns it is given", async () => {
      await session.step(135, "When user sets \"columnNames\" property of form viewer to \"AGE, SEX\"", () => setProperty(page, "columnNames", el("form viewer"), "AGE, SEX"));
      await session.step(136, "Then the \"fields shown\" reading of form viewer should be 2", () => readingIs(page, "fields shown", el("form viewer"), 2));
      await session.step(137, "And form viewer should have a \"field SEX\" area", () => hasArea(page, el("form viewer"), "field SEX"));
      await session.step(138, "And form viewer should have a \"field AGE\" area", () => hasArea(page, el("form viewer"), "field AGE"));
      await session.step(139, "And form viewer should not have a \"field RACE\" area", () => hasNoArea(page, el("form viewer"), "field RACE"));
      await session.step(140, "And the \"AGE\" reading of form viewer should be \"26\"", () => readingReads(page, "AGE", el("form viewer"), "26"));
    });
    await run.scenario("The field set survives a layout round-trip", async () => {
      await session.step(143, "Given user sets \"columnNames\" property of form viewer to \"AGE, HEIGHT, WEIGHT\"", () => setProperty(page, "columnNames", el("form viewer"), "AGE, HEIGHT, WEIGHT"));
      await session.step(144, "And the \"fields shown\" reading of form viewer should be 3", () => readingIs(page, "fields shown", el("form viewer"), 3));
      await session.step(145, "When user saves the layout of the current table view", () => saveLayout(page));
      await session.step(146, "And user clicks on close icon of form viewer", () => clickOn(page, el("close icon of form viewer")));
      await session.step(147, "Then form viewer should be absent", () => shouldBe(page, el("form viewer"), "absent"));
      await session.step(148, "When user loads the saved layout", () => loadLayout(page));
      await session.step(149, "Then form viewer should be visible", () => shouldBe(page, el("form viewer"), "visible"));
      await session.step(150, "And the \"fields shown\" reading of form viewer should be 3", () => readingIs(page, "fields shown", el("form viewer"), 3));
      await session.step(151, "And form viewer should have a \"field HEIGHT\" area", () => hasArea(page, el("form viewer"), "field HEIGHT"));
      await session.step(152, "And form viewer should have a \"field WEIGHT\" area", () => hasArea(page, el("form viewer"), "field WEIGHT"));
      await session.step(153, "And the \"AGE\" reading of form viewer should be \"26\"", () => readingReads(page, "AGE", el("form viewer"), "26"));
      await session.step(154, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A form saved to a file comes back from the file", async () => {
      await session.step(157, "Given user sets \"columnNames\" property of form viewer to \"AGE, SEX, RACE\"", () => setProperty(page, "columnNames", el("form viewer"), "AGE, SEX, RACE"));
      await session.step(158, "Then the \"fields shown\" reading of form viewer should be 3", () => readingIs(page, "fields shown", el("form viewer"), 3));
      await session.step(159, "When user downloads a file through arrow to bottom icon in form viewer", () => downloadThrough(page, el("arrow to bottom icon in form viewer")));
      await session.step(160, "Then the downloaded file should contain \"RACE\"", () => downloadedContains(page, "RACE"));
      await session.step(161, "And the downloaded file should not contain \"HEIGHT\"", () => downloadedNotContains(page, "HEIGHT"));
      await session.step(162, "When user sets \"columnNames\" property of form viewer to \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY\"", () => setProperty(page, "columnNames", el("form viewer"), "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"));
      await session.step(163, "Then the \"fields shown\" reading of form viewer should be 11", () => readingIs(page, "fields shown", el("form viewer"), 11));
      await session.step(164, "When user uploads the downloaded file through folder open icon in form viewer", () => uploadDownloaded(page, el("folder open icon in form viewer")));
      await session.step(165, "Then the \"fields shown\" reading of form viewer should be 3", () => readingIs(page, "fields shown", el("form viewer"), 3));
      await session.step(166, "And form viewer should have a \"field RACE\" area", () => hasArea(page, el("form viewer"), "field RACE"));
      await session.step(167, "And form viewer should have a \"field SEX\" area", () => hasArea(page, el("form viewer"), "field SEX"));
      await session.step(168, "And form viewer should not have a \"field HEIGHT\" area", () => hasNoArea(page, el("form viewer"), "field HEIGHT"));
      await session.step(169, "And the \"AGE\" reading of form viewer should be \"26\"", () => readingReads(page, "AGE", el("form viewer"), "26"));
      await session.step(170, "When user sets \"columnNames\" property of form viewer to \"AGE, HEIGHT, WEIGHT\"", () => setProperty(page, "columnNames", el("form viewer"), "AGE, HEIGHT, WEIGHT"));
      await session.step(171, "Then the \"fields shown\" reading of form viewer should be 3", () => readingIs(page, "fields shown", el("form viewer"), 3));
      await session.step(172, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A field dragged in design mode keeps its place, and the saved file brings the place back", async () => {
      await session.step(175, "Given user sets \"columnNames\" property of form viewer to \"AGE, SEX, RACE\"", () => setProperty(page, "columnNames", el("form viewer"), "AGE, SEX, RACE"));
      await session.step(176, "Then the \"fields shown\" reading of form viewer should be 3", () => readingIs(page, "fields shown", el("form viewer"), 3));
      await session.step(177, "And the \"design mode\" reading of form viewer should be \"false\"", () => readingReads(page, "design mode", el("form viewer"), "false"));
      await session.step(178, "And the \"field RACE\" area of form viewer should hang below the \"field AGE\" area", () => areaHangsBelow(page, "field RACE", el("form viewer"), "field AGE"));
      await session.step(179, "When user clicks on the \"design\" area of form viewer", () => clickArea(page, "design", el("form viewer")));
      await session.step(180, "Then the \"design mode\" reading of form viewer should be \"true\"", () => readingReads(page, "design mode", el("form viewer"), "true"));
      await session.step(181, "When user drags the \"field AGE\" area of form viewer by 160 pixels to the down", () => dragAreaBy(page, "field AGE", el("form viewer"), 160, "down"));
      await session.step(182, "Then the \"field AGE\" area of form viewer should hang below the \"field RACE\" area", () => areaHangsBelow(page, "field AGE", el("form viewer"), "field RACE"));
      await session.step(183, "And the \"field AGE\" area of form viewer should hang below the \"field SEX\" area", () => areaHangsBelow(page, "field AGE", el("form viewer"), "field SEX"));
      await session.step(184, "When user clicks on the \"design\" area of form viewer", () => clickArea(page, "design", el("form viewer")));
      await session.step(185, "Then the \"design mode\" reading of form viewer should be \"false\"", () => readingReads(page, "design mode", el("form viewer"), "false"));
      await session.step(186, "And the \"field AGE\" area of form viewer should hang below the \"field RACE\" area", () => areaHangsBelow(page, "field AGE", el("form viewer"), "field RACE"));
      await session.step(187, "When user downloads a file through arrow to bottom icon in form viewer", () => downloadThrough(page, el("arrow to bottom icon in form viewer")));
      await session.step(188, "And user sets \"columnNames\" property of form viewer to \"AGE, SEX, RACE, HEIGHT\"", () => setProperty(page, "columnNames", el("form viewer"), "AGE, SEX, RACE, HEIGHT"));
      await session.step(189, "Then the \"fields shown\" reading of form viewer should be 4", () => readingIs(page, "fields shown", el("form viewer"), 4));
      await session.step(190, "And the \"field RACE\" area of form viewer should hang below the \"field AGE\" area", () => areaHangsBelow(page, "field RACE", el("form viewer"), "field AGE"));
      await session.step(191, "When user uploads the downloaded file through folder open icon in form viewer", () => uploadDownloaded(page, el("folder open icon in form viewer")));
      await session.step(192, "Then the \"fields shown\" reading of form viewer should be 3", () => readingIs(page, "fields shown", el("form viewer"), 3));
      await session.step(193, "And the \"field AGE\" area of form viewer should hang below the \"field RACE\" area", () => areaHangsBelow(page, "field AGE", el("form viewer"), "field RACE"));
      await session.step(194, "When user sets \"columnNames\" property of form viewer to \"AGE, HEIGHT, WEIGHT\"", () => setProperty(page, "columnNames", el("form viewer"), "AGE, HEIGHT, WEIGHT"));
      await session.step(195, "Then the \"fields shown\" reading of form viewer should be 3", () => readingIs(page, "fields shown", el("form viewer"), 3));
      await session.step(196, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A removed column empties its field but keeps it on the form", async () => {
      await session.step(199, "When user removes \"WEIGHT\" column", () => removeColumn(page, "WEIGHT"));
      await session.step(200, "Then the table should not have a column \"WEIGHT\"", () => hasNoColumn(page, "WEIGHT"));
      await session.step(201, "And the \"fields shown\" reading of form viewer should be 3", () => readingIs(page, "fields shown", el("form viewer"), 3));
      await session.step(202, "And the \"WEIGHT\" reading of form viewer should be \"\"", () => readingReads(page, "WEIGHT", el("form viewer"), ""));
      await session.step(203, "And the \"AGE\" reading of form viewer should be \"26\"", () => readingReads(page, "AGE", el("form viewer"), "26"));
      await session.step(204, "And no errors should have been logged", () => noErrors(page));
      await session.step(205, "When user adds a calculated column \"WEIGHT\" with formula \"0\"", () => addCalculated(page, "WEIGHT", "0"));
      await session.step(206, "Then the table should have a column \"WEIGHT\"", () => hasColumn(page, "WEIGHT"));
    });
    await run.scenario("The arrows walk the filtered rows", async () => {
      await session.step(209, "Given user sets \"columnNames\" property of form viewer to \"USUBJID, AGE, SEX\"", () => setProperty(page, "columnNames", el("form viewer"), "USUBJID, AGE, SEX"));
      await session.step(210, "When user filters rows where \"SEX\" is \"M\"", () => filterTo(page, "SEX", "M"));
      await session.step(211, "And user makes row 4 current", () => makeRowCurrent(page, 4));
      await session.step(212, "Then the \"USUBJID\" reading of form viewer should be \"X0273T21000400002\"", () => readingReads(page, "USUBJID", el("form viewer"), "X0273T21000400002"));
      await session.step(213, "When user clicks on the \"next row\" area of form viewer", () => clickArea(page, "next row", el("form viewer")));
      await session.step(214, "Then the \"SEX\" reading of form viewer should be \"M\"", () => readingReads(page, "SEX", el("form viewer"), "M"));
      await session.step(215, "And the \"USUBJID\" reading of form viewer should be \"X0273T21000500008\"", () => readingReads(page, "USUBJID", el("form viewer"), "X0273T21000500008"));
      await session.step(216, "When user resets the filter", () => resetFilter(page));
      await session.step(217, "And user makes row 1 current", () => makeRowCurrent(page, 1));
    });
    await run.scenario("The viewer closes from its title bar", async () => {
      await session.step(220, "When user clicks on close icon of form viewer", () => clickOn(page, el("close icon of form viewer")));
      await session.step(221, "Then form viewer should be absent", () => shouldBe(page, el("form viewer"), "absent"));
      await session.step(222, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
