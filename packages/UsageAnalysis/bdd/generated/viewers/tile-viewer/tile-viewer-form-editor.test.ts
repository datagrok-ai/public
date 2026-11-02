/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tile-viewer/tile-viewer-form-editor.feature
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
import {addTileViewerWith, deleteLabelField, deleteValueField, designerLabelFields, designerValueFields} from '../../../bindings/tile-viewer.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, pressKey, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnCount, hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {addCalculated, removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {eventFired, hasArea, hasNoArea, listenFor, noErrors, pickFromAreaContextMenu, readingIs, readingReads, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {fieldsAsRemembered, fieldsRefilled, pickFromViewerMenu, readingContains, readingNotContains, rememberFields} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Tile viewer form designer", () => {
  const session = feature(test, "features/viewers/tile-viewer/tile-viewer-form-editor.feature", import.meta.url);
  test("Tile viewer form designer", {tag: ["@journey", "@viewers", "@realizes:viewers.tile-viewer", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(30, "Given user is logged in", () => loggedIn(page));
    await session.step(31, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(32, "And user adds a tile viewer with:", () => addTileViewerWith(page, [["Lanes Column Name","SEX"]]), [["Lanes Column Name","SEX"]]);
    await session.step(34, "Then tile viewer should be visible", () => shouldBe(page, el("tile viewer"), "visible"));
    await session.step(35, "And the \"lane names\" reading of tile viewer should be \"F, M\"", () => readingReads(page, "lane names", el("tile viewer"), "F, M"));
    await session.step(36, "And the \"fields shown\" reading of tile viewer should be 10", () => readingIs(page, "fields shown", el("tile viewer"), 10));
    await run.scenario("The card the designer opens on is the auto-generated one", async () => {
      await session.step(39, "Then the \"auto generate\" reading of tile viewer should be \"true\"", () => readingReads(page, "auto generate", el("tile viewer"), "true"));
      await session.step(40, "And the \"form designed\" reading of tile viewer should be \"false\"", () => readingReads(page, "form designed", el("tile viewer"), "false"));
      await session.step(41, "And the \"table\" reading of tile viewer should be \"demog-1000\"", () => readingReads(page, "table", el("tile viewer"), "demog-1000"));
      await session.step(42, "And the table should have 11 columns", () => columnCount(page, 11));
      await session.step(43, "And the \"fields\" reading of tile viewer should contain \"AGE\"", () => readingContains(page, "fields", el("tile viewer"), "AGE"));
      await session.step(44, "And the \"fields\" reading of tile viewer should contain \"DEMOG\"", () => readingContains(page, "fields", el("tile viewer"), "DEMOG"));
      await session.step(45, "And the \"fields\" reading of tile viewer should not contain \"SEVERITY\"", () => readingNotContains(page, "fields", el("tile viewer"), "SEVERITY"));
      await session.step(46, "And form designer should be absent", () => shouldBe(page, el("form designer"), "absent"));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Edit Form opens the designer on the viewer's own table", async () => {
      await session.step(50, "Given user listens for \"d4-tile-viewer-form-edit-request\" event on tile viewer", () => listenFor(page, "d4-tile-viewer-form-edit-request", el("tile viewer")));
      await session.step(51, "When user picks \"Edit Form...\" from the viewer menu of tile viewer", () => pickFromViewerMenu(page, "Edit Form...", el("tile viewer")));
      await session.step(52, "Then form designer should be visible", () => shouldBe(page, el("form designer"), "visible"));
      await session.step(53, "And the form designer should show value fields \"AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT\"", () => designerValueFields(page, "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"));
      await session.step(54, "And the form designer should show label fields \"AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT\"", () => designerLabelFields(page, "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"));
      await session.step(55, "When user clicks on \"CLOSE AND APPLY\" button", () => clickOn(page, el("\"CLOSE AND APPLY\" button")));
      await session.step(56, "Then form designer should be absent", () => shouldBe(page, el("form designer"), "absent"));
      await session.step(57, "And \"d4-tile-viewer-form-edit-request\" event should have fired on tile viewer", () => eventFired(page, "d4-tile-viewer-form-edit-request", el("tile viewer")));
      await session.step(58, "And the \"fields shown\" reading of tile viewer should be 10", () => readingIs(page, "fields shown", el("tile viewer"), 10));
      await session.step(59, "And the \"AGE of row 1\" reading of tile viewer should be \"26\"", () => readingReads(page, "AGE of row 1", el("tile viewer"), "26"));
      await session.step(60, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The column chooser counts the fields the card shows, and CANCEL keeps them", async () => {
      await session.step(63, "When user picks \"Edit Form...\" from the viewer menu of tile viewer", () => pickFromViewerMenu(page, "Edit Form...", el("tile viewer")));
      await session.step(64, "Then form designer should be visible", () => shouldBe(page, el("form designer"), "visible"));
      await session.step(65, "When user clicks on \"EDIT\" button", () => clickOn(page, el("\"EDIT\" button")));
      await session.step(66, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
      await session.step(67, "And \"All\" link in \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"All\" link in \"Select columns...\" dialog"), "visible"));
      await session.step(68, "And \"None\" link in \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"None\" link in \"Select columns...\" dialog"), "visible"));
      await session.step(69, "And \"10 checked\" text in \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"10 checked\" text in \"Select columns...\" dialog"), "visible"));
      await session.step(70, "When user clicks on \"CANCEL\" button in \"Select columns...\" dialog", () => clickOn(page, el("\"CANCEL\" button in \"Select columns...\" dialog")));
      await session.step(71, "Then \"Select columns...\" dialog should be absent", () => shouldBe(page, el("\"Select columns...\" dialog"), "absent"));
      await session.step(72, "And the form designer should show value fields \"AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT\"", () => designerValueFields(page, "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"));
      await session.step(73, "When user clicks on \"CLOSE AND APPLY\" button", () => clickOn(page, el("\"CLOSE AND APPLY\" button")));
      await session.step(74, "Then form designer should be absent", () => shouldBe(page, el("form designer"), "absent"));
      await session.step(75, "And the \"fields shown\" reading of tile viewer should be 10", () => readingIs(page, "fields shown", el("tile viewer"), 10));
      await session.step(76, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A deleted caption comes back with RESET, and goes with CLOSE AND APPLY", async () => {
      await session.step(79, "When user picks \"Edit Form...\" from the viewer menu of tile viewer", () => pickFromViewerMenu(page, "Edit Form...", el("tile viewer")));
      await session.step(80, "Then form designer should be visible", () => shouldBe(page, el("form designer"), "visible"));
      await session.step(81, "And the form designer should show label fields \"AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT\"", () => designerLabelFields(page, "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"));
      await session.step(82, "When user deletes the \"SEX\" label field in the form designer", () => deleteLabelField(page, "SEX"));
      await session.step(83, "Then the form designer should show label fields \"AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, STARTED, USUBJID, WEIGHT\"", () => designerLabelFields(page, "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, STARTED, USUBJID, WEIGHT"));
      await session.step(84, "And the form designer should show value fields \"AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT\"", () => designerValueFields(page, "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"));
      await session.step(85, "When user clicks on \"RESET\" button", () => clickOn(page, el("\"RESET\" button")));
      await session.step(86, "Then the form designer should show label fields \"AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT\"", () => designerLabelFields(page, "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"));
      await session.step(87, "And the form designer should show value fields \"AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT\"", () => designerValueFields(page, "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"));
      await session.step(88, "When user deletes the \"SEX\" label field in the form designer", () => deleteLabelField(page, "SEX"));
      await session.step(89, "And user clicks on \"CLOSE AND APPLY\" button", () => clickOn(page, el("\"CLOSE AND APPLY\" button")));
      await session.step(90, "Then form designer should be absent", () => shouldBe(page, el("form designer"), "absent"));
      await session.step(91, "And the \"fields shown\" reading of tile viewer should be 10", () => readingIs(page, "fields shown", el("tile viewer"), 10));
      await session.step(92, "And tile viewer should have a \"field SEX of row 1\" area", () => hasArea(page, el("tile viewer"), "field SEX of row 1"));
      await session.step(93, "And tile viewer should not have a \"label SEX of row 1\" area", () => hasNoArea(page, el("tile viewer"), "label SEX of row 1"));
      await session.step(94, "And the \"SEX of row 1\" reading of tile viewer should be \"F\"", () => readingReads(page, "SEX of row 1", el("tile viewer"), "F"));
      await session.step(95, "And the \"auto generate\" reading of tile viewer should be \"false\"", () => readingReads(page, "auto generate", el("tile viewer"), "false"));
      await session.step(96, "And the \"form designed\" reading of tile viewer should be \"true\"", () => readingReads(page, "form designed", el("tile viewer"), "true"));
      await session.step(97, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("RESET reverts to the state the designer opened in, not to the factory card", async () => {
      await session.step(100, "When user picks \"Edit Form...\" from the viewer menu of tile viewer", () => pickFromViewerMenu(page, "Edit Form...", el("tile viewer")));
      await session.step(101, "Then form designer should be visible", () => shouldBe(page, el("form designer"), "visible"));
      await session.step(102, "And the form designer should show label fields \"AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, STARTED, USUBJID, WEIGHT\"", () => designerLabelFields(page, "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, STARTED, USUBJID, WEIGHT"));
      await session.step(103, "When user deletes the \"AGE\" value field in the form designer", () => deleteValueField(page, "AGE"));
      await session.step(104, "Then the form designer should show value fields \"CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT\"", () => designerValueFields(page, "CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"));
      await session.step(105, "When user clicks on \"RESET\" button", () => clickOn(page, el("\"RESET\" button")));
      await session.step(106, "Then the form designer should show value fields \"AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT\"", () => designerValueFields(page, "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"));
      await session.step(107, "And the form designer should show label fields \"AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, STARTED, USUBJID, WEIGHT\"", () => designerLabelFields(page, "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, STARTED, USUBJID, WEIGHT"));
      await session.step(108, "When user deletes the \"AGE\" value field in the form designer", () => deleteValueField(page, "AGE"));
      await session.step(109, "And user clicks on \"CLOSE AND APPLY\" button", () => clickOn(page, el("\"CLOSE AND APPLY\" button")));
      await session.step(110, "Then form designer should be absent", () => shouldBe(page, el("form designer"), "absent"));
      await session.step(111, "And the \"fields shown\" reading of tile viewer should be 9", () => readingIs(page, "fields shown", el("tile viewer"), 9));
      await session.step(112, "And the \"fields\" reading of tile viewer should not contain \"AGE\"", () => readingNotContains(page, "fields", el("tile viewer"), "AGE"));
      await session.step(113, "And tile viewer should not have a \"field AGE of row 1\" area", () => hasNoArea(page, el("tile viewer"), "field AGE of row 1"));
      await session.step(115, "And tile viewer should not have a \"field AGE of row 4\" area", () => hasNoArea(page, el("tile viewer"), "field AGE of row 4"));
      await session.step(116, "And tile viewer should have a \"field DEMOG of row 4\" area", () => hasArea(page, el("tile viewer"), "field DEMOG of row 4"));
      await session.step(117, "And the \"form designed\" reading of tile viewer should be \"true\"", () => readingReads(page, "form designed", el("tile viewer"), "true"));
      await session.step(118, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("On a designed card a departing column empties its field and no column takes the slot", async () => {
      await session.step(121, "Then the \"form designed\" reading of tile viewer should be \"true\"", () => readingReads(page, "form designed", el("tile viewer"), "true"));
      await session.step(122, "And the \"fields shown\" reading of tile viewer should be 9", () => readingIs(page, "fields shown", el("tile viewer"), 9));
      await session.step(123, "And the \"fields\" reading of tile viewer should contain \"DEMOG\"", () => readingContains(page, "fields", el("tile viewer"), "DEMOG"));
      await session.step(124, "And the \"fields\" reading of tile viewer should not contain \"SEVERITY\"", () => readingNotContains(page, "fields", el("tile viewer"), "SEVERITY"));
      await session.step(125, "And the \"DEMOG of row 1\" reading of tile viewer should be \"26 C F\"", () => readingReads(page, "DEMOG of row 1", el("tile viewer"), "26 C F"));
      await session.step(126, "And the table should have 11 columns", () => columnCount(page, 11));
      await session.step(127, "When user remembers the fields of tile viewer", () => rememberFields(page, el("tile viewer")));
      await session.step(128, "And user picks \"Remove\" from the context menu of the \"field DEMOG of row 1\" area of tile viewer", () => pickFromAreaContextMenu(page, "Remove", "field DEMOG of row 1", el("tile viewer")));
      await session.step(129, "Then the table should not have a column \"DEMOG\"", () => hasNoColumn(page, "DEMOG"));
      await session.step(130, "And the table should have 10 columns", () => columnCount(page, 10));
      await session.step(131, "And the \"fields shown\" reading of tile viewer should be 9", () => readingIs(page, "fields shown", el("tile viewer"), 9));
      await session.step(132, "And the fields of tile viewer should be as remembered", () => fieldsAsRemembered(page, el("tile viewer")));
      await session.step(133, "And the \"fields\" reading of tile viewer should not contain \"SEVERITY\"", () => readingNotContains(page, "fields", el("tile viewer"), "SEVERITY"));
      await session.step(134, "And tile viewer should have a \"field DEMOG of row 1\" area", () => hasArea(page, el("tile viewer"), "field DEMOG of row 1"));
      await session.step(135, "And the \"DEMOG of row 1\" reading of tile viewer should be \"\"", () => readingReads(page, "DEMOG of row 1", el("tile viewer"), ""));
      await session.step(136, "And the \"auto generate\" reading of tile viewer should be \"false\"", () => readingReads(page, "auto generate", el("tile viewer"), "false"));
      await session.step(137, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Auto Generate rebuilds the card, and then a departing column's slot is refilled", async () => {
      await session.step(140, "Given user adds a calculated column \"DEMOG\" with formula \"${AGE}\"", () => addCalculated(page, "DEMOG", "${AGE}"));
      await session.step(141, "Then the table should have 11 columns", () => columnCount(page, 11));
      await session.step(142, "When user sets \"Auto Generate\" property of tile viewer to \"true\"", () => setProperty(page, "Auto Generate", el("tile viewer"), "true"));
      await session.step(143, "And user adds a calculated column \"TMP\" with formula \"1\"", () => addCalculated(page, "TMP", "1"));
      await session.step(144, "And user removes \"TMP\" column", () => removeColumn(page, "TMP"));
      await session.step(145, "Then the table should have 11 columns", () => columnCount(page, 11));
      await session.step(146, "And the \"auto generate\" reading of tile viewer should be \"true\"", () => readingReads(page, "auto generate", el("tile viewer"), "true"));
      await session.step(147, "And the \"form designed\" reading of tile viewer should be \"false\"", () => readingReads(page, "form designed", el("tile viewer"), "false"));
      await session.step(148, "And the \"fields shown\" reading of tile viewer should be 10", () => readingIs(page, "fields shown", el("tile viewer"), 10));
      await session.step(149, "And the \"fields\" reading of tile viewer should contain \"AGE\"", () => readingContains(page, "fields", el("tile viewer"), "AGE"));
      await session.step(150, "When user remembers the fields of tile viewer", () => rememberFields(page, el("tile viewer")));
      await session.step(151, "And user picks \"Remove\" from the context menu of the \"field AGE of row 1\" area of tile viewer", () => pickFromAreaContextMenu(page, "Remove", "field AGE of row 1", el("tile viewer")));
      await session.step(152, "Then the table should not have a column \"AGE\"", () => hasNoColumn(page, "AGE"));
      await session.step(153, "And the table should have 10 columns", () => columnCount(page, 10));
      await session.step(154, "And the \"fields shown\" reading of tile viewer should be 10", () => readingIs(page, "fields shown", el("tile viewer"), 10));
      await session.step(155, "And the fields of tile viewer should have refilled the freed slot", () => fieldsRefilled(page, el("tile viewer")));
      await session.step(156, "And the \"fields\" reading of tile viewer should not contain \"AGE\"", () => readingNotContains(page, "fields", el("tile viewer"), "AGE"));
      await session.step(157, "And tile viewer should not have a \"field AGE of row 1\" area", () => hasNoArea(page, el("tile viewer"), "field AGE of row 1"));
      await session.step(158, "And the \"auto generate\" reading of tile viewer should be \"true\"", () => readingReads(page, "auto generate", el("tile viewer"), "true"));
      await session.step(159, "And no errors should have been logged", () => noErrors(page));
      await session.step(160, "When user adds a calculated column \"AGE\" with formula \"${HEIGHT}\"", () => addCalculated(page, "AGE", "${HEIGHT}"));
      await session.step(161, "Then the table should have 11 columns", () => columnCount(page, 11));
    });
    await run.scenario("Escape in the designer leaves it open with the edits it holds", async () => {
      await session.step(165, "When user picks \"Edit Form...\" from the viewer menu of tile viewer", () => pickFromViewerMenu(page, "Edit Form...", el("tile viewer")));
      await session.step(166, "Then form designer should be visible", () => shouldBe(page, el("form designer"), "visible"));
      await session.step(167, "When user deletes the \"WEIGHT\" value field in the form designer", () => deleteValueField(page, "WEIGHT"));
      await session.step(168, "And user presses Escape", () => pressKey(page, "Escape"));
      await session.step(169, "Then form designer should be visible", () => shouldBe(page, el("form designer"), "visible"));
      await session.step(170, "When user clicks on \"CLOSE AND APPLY\" button", () => clickOn(page, el("\"CLOSE AND APPLY\" button")));
      await session.step(171, "Then form designer should be absent", () => shouldBe(page, el("form designer"), "absent"));
      await session.step(172, "And the \"fields\" reading of tile viewer should not contain \"WEIGHT\"", () => readingNotContains(page, "fields", el("tile viewer"), "WEIGHT"));
      await session.step(173, "And the \"form designed\" reading of tile viewer should be \"true\"", () => readingReads(page, "form designed", el("tile viewer"), "true"));
      await session.step(174, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
