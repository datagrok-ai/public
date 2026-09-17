/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/form/form.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.form]
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
import {clickOn, focusOn, pressKeyIn, shouldBe, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs, hasColumn, hasNoColumn, makeRowCurrent, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {addCalculated, clearSelection, colorLinear, colorOff, filterTo, noneSelected, onlyOfSelected, removeColumn, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, enterIntoArea, hasArea, hasNoArea, hoverArea, loadLayout, noErrors, pickFromContextMenu, pointerAway, propertyShouldBe, readingDiffers, readingIs, readingReads, saveLayout, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Form viewer", () => {
  const session = feature(test, "features/viewers/form/form.feature", import.meta.url);
  test("Form viewer", {tag: ["@journey", "@viewers", "@realizes:viewers.form"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 14, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a form viewer", () => addViewer(page, "form"));
    await session.step(16, "Then form viewer should be visible", () => shouldBe(page, el("form viewer"), "visible"));
    await session.step(17, "And the \"fields shown\" reading of form viewer should be 11", () => readingIs(page, "fields shown", el("form viewer"), 11));
    await session.step(18, "And form viewer should have a \"field USUBJID\" area", () => hasArea(page, el("form viewer"), "field USUBJID"));
    await session.step(19, "And form viewer should have a \"field SEVERITY\" area", () => hasArea(page, el("form viewer"), "field SEVERITY"));
    await run.scenario("The fields show the current row", async () => {
      await session.step(22, "Then the \"row\" reading of form viewer should be 1", () => readingIs(page, "row", el("form viewer"), 1));
      await session.step(23, "And the \"USUBJID\" reading of form viewer should be \"X0273T21000300003\"", () => readingReads(page, "USUBJID", el("form viewer"), "X0273T21000300003"));
      await session.step(24, "And the \"AGE\" reading of form viewer should be \"26\"", () => readingReads(page, "AGE", el("form viewer"), "26"));
      await session.step(25, "And the \"SEX\" reading of form viewer should be \"F\"", () => readingReads(page, "SEX", el("form viewer"), "F"));
      await session.step(26, "And the \"RACE\" reading of form viewer should be \"Caucasian\"", () => readingReads(page, "RACE", el("form viewer"), "Caucasian"));
      await session.step(27, "And the \"DIS_POP\" reading of form viewer should be \"Indigestion\"", () => readingReads(page, "DIS_POP", el("form viewer"), "Indigestion"));
      await session.step(28, "And form viewer should have a \"field AGE\" area", () => hasArea(page, el("form viewer"), "field AGE"));
      await session.step(29, "And form viewer should have a \"label AGE\" area", () => hasArea(page, el("form viewer"), "label AGE"));
      await session.step(30, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column's colour coding reaches its field", async () => {
      await session.step(33, "Then the \"background of AGE\" reading of form viewer should be \"\"", () => readingReads(page, "background of AGE", el("form viewer"), ""));
      await session.step(34, "When user colors \"AGE\" column linearly from \"#FF0000\" to \"#0000FF\"", () => colorLinear(page, "AGE", "#FF0000", "#0000FF"));
      await session.step(35, "Then the \"background of AGE\" reading of form viewer should differ from before", () => readingDiffers(page, "background of AGE", el("form viewer")));
      await session.step(36, "When user removes the coloring of \"AGE\" column", () => colorOff(page, "AGE"));
      await session.step(37, "Then the \"background of AGE\" reading of form viewer should be \"\"", () => readingReads(page, "background of AGE", el("form viewer"), ""));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The arrows walk the rows", async () => {
      await session.step(41, "When user clicks on the \"next row\" area of form viewer", () => clickArea(page, "next row", el("form viewer")));
      await session.step(42, "Then row 2 should be current", () => currentRowIs(page, 2));
      await session.step(43, "And the \"USUBJID\" reading of form viewer should be \"X0273T21000300005\"", () => readingReads(page, "USUBJID", el("form viewer"), "X0273T21000300005"));
      await session.step(44, "And the \"AGE\" reading of form viewer should be \"30\"", () => readingReads(page, "AGE", el("form viewer"), "30"));
      await session.step(45, "When user clicks on the \"previous row\" area of form viewer", () => clickArea(page, "previous row", el("form viewer")));
      await session.step(46, "Then row 1 should be current", () => currentRowIs(page, 1));
      await session.step(47, "And the \"USUBJID\" reading of form viewer should be \"X0273T21000300003\"", () => readingReads(page, "USUBJID", el("form viewer"), "X0273T21000300003"));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The form follows the grid's current row", async () => {
      await session.step(51, "When user makes row 43 current", () => makeRowCurrent(page, 43));
      await session.step(52, "Then the \"row\" reading of form viewer should be 43", () => readingIs(page, "row", el("form viewer"), 43));
      await session.step(53, "And the \"USUBJID\" reading of form viewer should be \"X0273T21002400007\"", () => readingReads(page, "USUBJID", el("form viewer"), "X0273T21002400007"));
      await session.step(54, "And the \"AGE\" reading of form viewer should be \"54\"", () => readingReads(page, "AGE", el("form viewer"), "54"));
      await session.step(55, "When user makes row 1 current", () => makeRowCurrent(page, 1));
      await session.step(56, "Then the \"row\" reading of form viewer should be 1", () => readingIs(page, "row", el("form viewer"), 1));
    });
    await run.scenario("The row selector toggles the row's selection", async () => {
      await session.step(59, "Given user clears the row selection", () => clearSelection(page));
      await session.step(60, "When user clicks on the \"row selector\" area of form viewer", () => clickArea(page, "row selector", el("form viewer")));
      await session.step(61, "Then only rows where \"USUBJID\" is \"X0273T21000300003\" should be selected", () => onlyOfSelected(page, "USUBJID", "X0273T21000300003"));
      await session.step(62, "And the \"row selected\" reading of form viewer should be \"true\"", () => readingReads(page, "row selected", el("form viewer"), "true"));
      await session.step(63, "And \"square\" icon in form viewer should be selected", () => shouldBe(page, el("\"square\" icon in form viewer"), "selected"));
      await session.step(64, "When user clicks on the \"row selector\" area of form viewer", () => clickArea(page, "row selector", el("form viewer")));
      await session.step(65, "Then no rows should be selected", () => noneSelected(page));
      await session.step(66, "And the \"row selected\" reading of form viewer should be \"false\"", () => readingReads(page, "row selected", el("form viewer"), "false"));
      await session.step(67, "And \"square\" icon in form viewer should not be selected", () => shouldNotBe(page, el("\"square\" icon in form viewer"), "selected"));
    });
    await run.scenario("Track Row: the form follows the mouse-over row, or nothing", async () => {
      await session.step(70, "When user picks \"Track Row > Mouse Over\" from the context menu of form viewer", () => pickFromContextMenu(page, "Track Row > Mouse Over", el("form viewer")));
      await session.step(71, "Then \"Sync Mode\" property of form viewer should be \"Mouse Over\"", () => propertyShouldBe(page, "Sync Mode", el("form viewer"), "Mouse Over"));
      await session.step(72, "When user hovers over the \"cell 2 of AGE\" area of grid", () => hoverArea(page, "cell 2 of AGE", el("grid")));
      await session.step(73, "Then the \"USUBJID\" reading of form viewer should be \"X0273T21000300005\"", () => readingReads(page, "USUBJID", el("form viewer"), "X0273T21000300005"));
      await session.step(74, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(75, "And user picks \"Track Row > None\" from the context menu of form viewer", () => pickFromContextMenu(page, "Track Row > None", el("form viewer")));
      await session.step(76, "Then \"Sync Mode\" property of form viewer should be \"None\"", () => propertyShouldBe(page, "Sync Mode", el("form viewer"), "None"));
      await session.step(77, "When user makes row 5 current", () => makeRowCurrent(page, 5));
      await session.step(78, "Then the \"USUBJID\" reading of form viewer should be \"X0273T21000300005\"", () => readingReads(page, "USUBJID", el("form viewer"), "X0273T21000300005"));
      await session.step(79, "When user picks \"Track Row > Current\" from the context menu of form viewer", () => pickFromContextMenu(page, "Track Row > Current", el("form viewer")));
      await session.step(80, "Then \"Sync Mode\" property of form viewer should be \"Current\"", () => propertyShouldBe(page, "Sync Mode", el("form viewer"), "Current"));
      await session.step(81, "And the \"USUBJID\" reading of form viewer should be \"X0273T21000500006\"", () => readingReads(page, "USUBJID", el("form viewer"), "X0273T21000500006"));
      await session.step(82, "When user makes row 1 current", () => makeRowCurrent(page, 1));
    });
    await run.scenario("The keyboard walks and selects rows", async () => {
      await session.step(85, "Given user clears the row selection", () => clearSelection(page));
      await session.step(86, "When user focuses on form viewer", () => focusOn(page, el("form viewer")));
      await session.step(87, "And user presses ArrowRight in form viewer", () => pressKeyIn(page, "ArrowRight", el("form viewer")));
      await session.step(88, "Then row 2 should be current", () => currentRowIs(page, 2));
      await session.step(89, "When user presses ArrowDown in form viewer", () => pressKeyIn(page, "ArrowDown", el("form viewer")));
      await session.step(90, "Then row 3 should be current", () => currentRowIs(page, 3));
      await session.step(91, "When user presses ArrowLeft in form viewer", () => pressKeyIn(page, "ArrowLeft", el("form viewer")));
      await session.step(92, "Then row 2 should be current", () => currentRowIs(page, 2));
      await session.step(93, "When user presses ArrowUp in form viewer", () => pressKeyIn(page, "ArrowUp", el("form viewer")));
      await session.step(94, "Then row 1 should be current", () => currentRowIs(page, 1));
      await session.step(95, "When user presses Space in form viewer", () => pressKeyIn(page, "Space", el("form viewer")));
      await session.step(96, "Then only rows where \"USUBJID\" is \"X0273T21000300003\" should be selected", () => onlyOfSelected(page, "USUBJID", "X0273T21000300003"));
      await session.step(97, "When user clears the row selection", () => clearSelection(page));
    });
    await run.scenario("A field writes its cell only while the form is editable", async () => {
      await session.step(100, "Then the \"editable\" reading of form viewer should be \"false\"", () => readingReads(page, "editable", el("form viewer"), "false"));
      await session.step(101, "When user clicks on the \"edit\" area of form viewer", () => clickArea(page, "edit", el("form viewer")));
      await session.step(102, "Then the \"editable\" reading of form viewer should be \"true\"", () => readingReads(page, "editable", el("form viewer"), "true"));
      await session.step(103, "When user enters \"31\" into the \"field AGE\" area of form viewer", () => enterIntoArea(page, "31", "field AGE", el("form viewer")));
      await session.step(104, "Then the value of \"AGE\" column in row 1 should be \"31\"", () => valueInRow(page, "AGE", 1, "31"));
      await session.step(105, "And the \"AGE\" reading of form viewer should be \"31\"", () => readingReads(page, "AGE", el("form viewer"), "31"));
      await session.step(106, "When user enters \"26\" into the \"field AGE\" area of form viewer", () => enterIntoArea(page, "26", "field AGE", el("form viewer")));
      await session.step(107, "And user clicks on the \"edit\" area of form viewer", () => clickArea(page, "edit", el("form viewer")));
      await session.step(108, "Then the \"editable\" reading of form viewer should be \"false\"", () => readingReads(page, "editable", el("form viewer"), "false"));
      await session.step(109, "And the value of \"AGE\" column in row 1 should be \"26\"", () => valueInRow(page, "AGE", 1, "26"));
      await session.step(110, "When user enters \"99\" into the \"field AGE\" area of form viewer", () => enterIntoArea(page, "99", "field AGE", el("form viewer")));
      await session.step(111, "Then the value of \"AGE\" column in row 1 should be \"26\"", () => valueInRow(page, "AGE", 1, "26"));
    });
    await run.scenario("The toolbar shows what the properties allow", async () => {
      await session.step(114, "When user sets \"Show Next Row Arrow\" property of form viewer to \"false\"", () => setProperty(page, "Show Next Row Arrow", el("form viewer"), "false"));
      await session.step(115, "Then form viewer should not have a \"next row\" area", () => hasNoArea(page, el("form viewer"), "next row"));
      await session.step(116, "And form viewer should have a \"previous row\" area", () => hasArea(page, el("form viewer"), "previous row"));
      await session.step(117, "When user sets \"Show Row Selector\" property of form viewer to \"false\"", () => setProperty(page, "Show Row Selector", el("form viewer"), "false"));
      await session.step(118, "Then form viewer should not have a \"row selector\" area", () => hasNoArea(page, el("form viewer"), "row selector"));
      await session.step(119, "When user sets properties of form viewer:", () => setProperties(page, el("form viewer"), [["Show Next Row Arrow","true"],["Show Row Selector","true"]]));
      await session.step(122, "Then form viewer should have a \"next row\" area", () => hasArea(page, el("form viewer"), "next row"));
      await session.step(123, "And form viewer should have a \"row selector\" area", () => hasArea(page, el("form viewer"), "row selector"));
      await session.step(124, "When user sets \"Show Navigation\" property of form viewer to \"false\"", () => setProperty(page, "Show Navigation", el("form viewer"), "false"));
      await session.step(125, "Then form viewer should not have a \"next row\" area", () => hasNoArea(page, el("form viewer"), "next row"));
      await session.step(126, "And form viewer should not have a \"select columns\" area", () => hasNoArea(page, el("form viewer"), "select columns"));
      await session.step(127, "When user sets \"Show Navigation\" property of form viewer to \"true\"", () => setProperty(page, "Show Navigation", el("form viewer"), "true"));
      await session.step(128, "Then form viewer should have a \"next row\" area", () => hasArea(page, el("form viewer"), "next row"));
    });
    await run.scenario("The field set follows the columns it is given", async () => {
      await session.step(131, "When user sets \"columnNames\" property of form viewer to \"AGE, SEX\"", () => setProperty(page, "columnNames", el("form viewer"), "AGE, SEX"));
      await session.step(132, "Then the \"fields shown\" reading of form viewer should be 2", () => readingIs(page, "fields shown", el("form viewer"), 2));
      await session.step(133, "And form viewer should have a \"field SEX\" area", () => hasArea(page, el("form viewer"), "field SEX"));
      await session.step(134, "And form viewer should have a \"field AGE\" area", () => hasArea(page, el("form viewer"), "field AGE"));
      await session.step(135, "And form viewer should not have a \"field RACE\" area", () => hasNoArea(page, el("form viewer"), "field RACE"));
      await session.step(136, "And the \"AGE\" reading of form viewer should be \"26\"", () => readingReads(page, "AGE", el("form viewer"), "26"));
    });
    await run.scenario("The field set survives a layout round-trip", async () => {
      await session.step(139, "Given user sets \"columnNames\" property of form viewer to \"AGE, HEIGHT, WEIGHT\"", () => setProperty(page, "columnNames", el("form viewer"), "AGE, HEIGHT, WEIGHT"));
      await session.step(140, "And the \"fields shown\" reading of form viewer should be 3", () => readingIs(page, "fields shown", el("form viewer"), 3));
      await session.step(141, "When user saves the layout of the current table view", () => saveLayout(page));
      await session.step(142, "And user clicks on close icon of form viewer", () => clickOn(page, el("close icon of form viewer")));
      await session.step(143, "Then form viewer should be absent", () => shouldBe(page, el("form viewer"), "absent"));
      await session.step(144, "When user loads the saved layout", () => loadLayout(page));
      await session.step(145, "Then form viewer should be visible", () => shouldBe(page, el("form viewer"), "visible"));
      await session.step(146, "And the \"fields shown\" reading of form viewer should be 3", () => readingIs(page, "fields shown", el("form viewer"), 3));
      await session.step(147, "And form viewer should have a \"field HEIGHT\" area", () => hasArea(page, el("form viewer"), "field HEIGHT"));
      await session.step(148, "And form viewer should have a \"field WEIGHT\" area", () => hasArea(page, el("form viewer"), "field WEIGHT"));
      await session.step(149, "And the \"AGE\" reading of form viewer should be \"26\"", () => readingReads(page, "AGE", el("form viewer"), "26"));
      await session.step(150, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A removed column empties its field but keeps it on the form", async () => {
      await session.step(153, "When user removes \"WEIGHT\" column", () => removeColumn(page, "WEIGHT"));
      await session.step(154, "Then the table should not have a column \"WEIGHT\"", () => hasNoColumn(page, "WEIGHT"));
      await session.step(155, "And the \"fields shown\" reading of form viewer should be 3", () => readingIs(page, "fields shown", el("form viewer"), 3));
      await session.step(156, "And the \"WEIGHT\" reading of form viewer should be \"\"", () => readingReads(page, "WEIGHT", el("form viewer"), ""));
      await session.step(157, "And the \"AGE\" reading of form viewer should be \"26\"", () => readingReads(page, "AGE", el("form viewer"), "26"));
      await session.step(158, "And no errors should have been logged", () => noErrors(page));
      await session.step(159, "When user adds a calculated column \"WEIGHT\" with formula \"0\"", () => addCalculated(page, "WEIGHT", "0"));
      await session.step(160, "Then the table should have a column \"WEIGHT\"", () => hasColumn(page, "WEIGHT"));
    });
    await run.scenario("The arrows walk the filtered rows", async () => {
      await session.step(163, "Given user sets \"columnNames\" property of form viewer to \"USUBJID, AGE, SEX\"", () => setProperty(page, "columnNames", el("form viewer"), "USUBJID, AGE, SEX"));
      await session.step(164, "When user filters rows where \"SEX\" is \"M\"", () => filterTo(page, "SEX", "M"));
      await session.step(165, "And user makes row 4 current", () => makeRowCurrent(page, 4));
      await session.step(166, "Then the \"USUBJID\" reading of form viewer should be \"X0273T21000400002\"", () => readingReads(page, "USUBJID", el("form viewer"), "X0273T21000400002"));
      await session.step(167, "When user clicks on the \"next row\" area of form viewer", () => clickArea(page, "next row", el("form viewer")));
      await session.step(168, "Then the \"SEX\" reading of form viewer should be \"M\"", () => readingReads(page, "SEX", el("form viewer"), "M"));
      await session.step(169, "And the \"USUBJID\" reading of form viewer should be \"X0273T21000500008\"", () => readingReads(page, "USUBJID", el("form viewer"), "X0273T21000500008"));
      await session.step(170, "When user resets the filter", () => resetFilter(page));
      await session.step(171, "And user makes row 1 current", () => makeRowCurrent(page, 1));
    });
    await run.scenario("The viewer closes from its title bar", async () => {
      await session.step(174, "When user clicks on close icon of form viewer", () => clickOn(page, el("close icon of form viewer")));
      await session.step(175, "Then form viewer should be absent", () => shouldBe(page, el("form viewer"), "absent"));
      await session.step(176, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
