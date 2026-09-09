/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tile-viewer/tile-viewer-form-editor.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tile-viewer]
--- */
import {test} from '@playwright/test';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {addTileViewerWith, deleteLabelField, deleteValueField, designerLabelFields, designerValueFields, fieldsAsRemembered, fieldsRefilled, pickFromViewerMenu, readingContains, readingNotContains, rememberFields} from '../../../bindings/tile-viewer.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnCount, hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {addCalculated, removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {eventFired, hasArea, hasNoArea, listenFor, noErrors, pickFromAreaContextMenu, readingIs, readingReads, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Tile viewer form designer", () => {
  const session = feature(test, "features/viewers/tile-viewer/tile-viewer-form-editor.feature", import.meta.url);
  test("Tile viewer form designer", {tag: ["@journey", "@viewers", "@realizes:viewers.tile-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(26, "Given user is logged in", () => loggedIn(page));
    await session.step(27, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(28, "And user adds a tile viewer with:", () => addTileViewerWith(page, [["Lanes Column Name","SEX"]]));
    await session.step(30, "Then tile viewer should be visible", () => shouldBe(page, el("tile viewer"), "visible"));
    await session.step(31, "And the \"lane names\" reading of tile viewer should be \"F, M\"", () => readingReads(page, "lane names", el("tile viewer"), "F, M"));
    await session.step(32, "And the \"fields shown\" reading of tile viewer should be 10", () => readingIs(page, "fields shown", el("tile viewer"), 10));
    await run.scenario("The card the designer opens on is the auto-generated one", async () => {
      await session.step(35, "Then the \"auto generate\" reading of tile viewer should be \"true\"", () => readingReads(page, "auto generate", el("tile viewer"), "true"));
      await session.step(36, "And the \"form designed\" reading of tile viewer should be \"false\"", () => readingReads(page, "form designed", el("tile viewer"), "false"));
      await session.step(37, "And the \"table\" reading of tile viewer should be \"demog-1000\"", () => readingReads(page, "table", el("tile viewer"), "demog-1000"));
      await session.step(38, "And the table should have 11 columns", () => columnCount(page, 11));
      await session.step(39, "And the \"fields\" reading of tile viewer should contain \"AGE\"", () => readingContains(page, "fields", el("tile viewer"), "AGE"));
      await session.step(40, "And the \"fields\" reading of tile viewer should contain \"DEMOG\"", () => readingContains(page, "fields", el("tile viewer"), "DEMOG"));
      await session.step(41, "And the \"fields\" reading of tile viewer should not contain \"SEVERITY\"", () => readingNotContains(page, "fields", el("tile viewer"), "SEVERITY"));
      await session.step(42, "And form designer should be absent", () => shouldBe(page, el("form designer"), "absent"));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Edit Form opens the designer on the viewer's own table", async () => {
      await session.step(46, "Given user listens for \"d4-tile-viewer-form-edit-request\" event on tile viewer", () => listenFor(page, "d4-tile-viewer-form-edit-request", el("tile viewer")));
      await session.step(47, "When user picks \"Edit Form...\" from the viewer menu of tile viewer", () => pickFromViewerMenu(page, "Edit Form...", el("tile viewer")));
      await session.step(48, "Then form designer should be visible", () => shouldBe(page, el("form designer"), "visible"));
      await session.step(49, "And the form designer should show value fields \"AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT\"", () => designerValueFields(page, "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"));
      await session.step(50, "And the form designer should show label fields \"AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT\"", () => designerLabelFields(page, "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"));
      await session.step(51, "When user clicks on \"CLOSE AND APPLY\" button", () => clickOn(page, el("\"CLOSE AND APPLY\" button")));
      await session.step(52, "Then form designer should be absent", () => shouldBe(page, el("form designer"), "absent"));
      await session.step(53, "And \"d4-tile-viewer-form-edit-request\" event should have fired on tile viewer", () => eventFired(page, "d4-tile-viewer-form-edit-request", el("tile viewer")));
      await session.step(54, "And the \"fields shown\" reading of tile viewer should be 10", () => readingIs(page, "fields shown", el("tile viewer"), 10));
      await session.step(55, "And the \"AGE of row 1\" reading of tile viewer should be \"26\"", () => readingReads(page, "AGE of row 1", el("tile viewer"), "26"));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The column chooser counts the fields the card shows, and CANCEL keeps them", async () => {
      await session.step(59, "When user picks \"Edit Form...\" from the viewer menu of tile viewer", () => pickFromViewerMenu(page, "Edit Form...", el("tile viewer")));
      await session.step(60, "Then form designer should be visible", () => shouldBe(page, el("form designer"), "visible"));
      await session.step(61, "When user clicks on \"EDIT\" button", () => clickOn(page, el("\"EDIT\" button")));
      await session.step(62, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
      await session.step(63, "And \"All\" link in \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"All\" link in \"Select columns...\" dialog"), "visible"));
      await session.step(64, "And \"None\" link in \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"None\" link in \"Select columns...\" dialog"), "visible"));
      await session.step(65, "And \"10 checked\" text in \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"10 checked\" text in \"Select columns...\" dialog"), "visible"));
      await session.step(66, "When user clicks on \"CANCEL\" button in \"Select columns...\" dialog", () => clickOn(page, el("\"CANCEL\" button in \"Select columns...\" dialog")));
      await session.step(67, "Then \"Select columns...\" dialog should be absent", () => shouldBe(page, el("\"Select columns...\" dialog"), "absent"));
      await session.step(68, "And the form designer should show value fields \"AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT\"", () => designerValueFields(page, "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"));
      await session.step(69, "When user clicks on \"CLOSE AND APPLY\" button", () => clickOn(page, el("\"CLOSE AND APPLY\" button")));
      await session.step(70, "Then form designer should be absent", () => shouldBe(page, el("form designer"), "absent"));
      await session.step(71, "And the \"fields shown\" reading of tile viewer should be 10", () => readingIs(page, "fields shown", el("tile viewer"), 10));
      await session.step(72, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A deleted caption comes back with RESET, and goes with CLOSE AND APPLY", async () => {
      await session.step(75, "When user picks \"Edit Form...\" from the viewer menu of tile viewer", () => pickFromViewerMenu(page, "Edit Form...", el("tile viewer")));
      await session.step(76, "Then form designer should be visible", () => shouldBe(page, el("form designer"), "visible"));
      await session.step(77, "And the form designer should show label fields \"AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT\"", () => designerLabelFields(page, "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"));
      await session.step(78, "When user deletes the \"SEX\" label field in the form designer", () => deleteLabelField(page, "SEX"));
      await session.step(79, "Then the form designer should show label fields \"AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, STARTED, USUBJID, WEIGHT\"", () => designerLabelFields(page, "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, STARTED, USUBJID, WEIGHT"));
      await session.step(80, "And the form designer should show value fields \"AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT\"", () => designerValueFields(page, "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"));
      await session.step(81, "When user clicks on \"RESET\" button", () => clickOn(page, el("\"RESET\" button")));
      await session.step(82, "Then the form designer should show label fields \"AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT\"", () => designerLabelFields(page, "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"));
      await session.step(83, "And the form designer should show value fields \"AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT\"", () => designerValueFields(page, "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"));
      await session.step(84, "When user deletes the \"SEX\" label field in the form designer", () => deleteLabelField(page, "SEX"));
      await session.step(85, "And user clicks on \"CLOSE AND APPLY\" button", () => clickOn(page, el("\"CLOSE AND APPLY\" button")));
      await session.step(86, "Then form designer should be absent", () => shouldBe(page, el("form designer"), "absent"));
      await session.step(87, "And the \"fields shown\" reading of tile viewer should be 10", () => readingIs(page, "fields shown", el("tile viewer"), 10));
      await session.step(88, "And tile viewer should have a \"field SEX of row 1\" area", () => hasArea(page, el("tile viewer"), "field SEX of row 1"));
      await session.step(89, "And tile viewer should not have a \"label SEX of row 1\" area", () => hasNoArea(page, el("tile viewer"), "label SEX of row 1"));
      await session.step(90, "And the \"SEX of row 1\" reading of tile viewer should be \"F\"", () => readingReads(page, "SEX of row 1", el("tile viewer"), "F"));
      await session.step(91, "And the \"auto generate\" reading of tile viewer should be \"false\"", () => readingReads(page, "auto generate", el("tile viewer"), "false"));
      await session.step(92, "And the \"form designed\" reading of tile viewer should be \"true\"", () => readingReads(page, "form designed", el("tile viewer"), "true"));
      await session.step(93, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("RESET reverts to the state the designer opened in, not to the factory card", async () => {
      await session.step(96, "When user picks \"Edit Form...\" from the viewer menu of tile viewer", () => pickFromViewerMenu(page, "Edit Form...", el("tile viewer")));
      await session.step(97, "Then form designer should be visible", () => shouldBe(page, el("form designer"), "visible"));
      await session.step(98, "And the form designer should show label fields \"AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, STARTED, USUBJID, WEIGHT\"", () => designerLabelFields(page, "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, STARTED, USUBJID, WEIGHT"));
      await session.step(99, "When user deletes the \"AGE\" value field in the form designer", () => deleteValueField(page, "AGE"));
      await session.step(100, "Then the form designer should show value fields \"CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT\"", () => designerValueFields(page, "CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"));
      await session.step(101, "When user clicks on \"RESET\" button", () => clickOn(page, el("\"RESET\" button")));
      await session.step(102, "Then the form designer should show value fields \"AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT\"", () => designerValueFields(page, "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"));
      await session.step(103, "And the form designer should show label fields \"AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, STARTED, USUBJID, WEIGHT\"", () => designerLabelFields(page, "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, STARTED, USUBJID, WEIGHT"));
      await session.step(104, "When user deletes the \"AGE\" value field in the form designer", () => deleteValueField(page, "AGE"));
      await session.step(105, "And user clicks on \"CLOSE AND APPLY\" button", () => clickOn(page, el("\"CLOSE AND APPLY\" button")));
      await session.step(106, "Then form designer should be absent", () => shouldBe(page, el("form designer"), "absent"));
      await session.step(107, "And the \"fields shown\" reading of tile viewer should be 9", () => readingIs(page, "fields shown", el("tile viewer"), 9));
      await session.step(108, "And the \"fields\" reading of tile viewer should not contain \"AGE\"", () => readingNotContains(page, "fields", el("tile viewer"), "AGE"));
      await session.step(109, "And tile viewer should not have a \"field AGE of row 1\" area", () => hasNoArea(page, el("tile viewer"), "field AGE of row 1"));
      await session.step(110, "And the \"form designed\" reading of tile viewer should be \"true\"", () => readingReads(page, "form designed", el("tile viewer"), "true"));
      await session.step(111, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("On a designed card a departing column empties its field and no column takes the slot", async () => {
      await session.step(114, "Then the \"form designed\" reading of tile viewer should be \"true\"", () => readingReads(page, "form designed", el("tile viewer"), "true"));
      await session.step(115, "And the \"fields shown\" reading of tile viewer should be 9", () => readingIs(page, "fields shown", el("tile viewer"), 9));
      await session.step(116, "And the \"fields\" reading of tile viewer should contain \"DEMOG\"", () => readingContains(page, "fields", el("tile viewer"), "DEMOG"));
      await session.step(117, "And the \"fields\" reading of tile viewer should not contain \"SEVERITY\"", () => readingNotContains(page, "fields", el("tile viewer"), "SEVERITY"));
      await session.step(118, "And the \"DEMOG of row 1\" reading of tile viewer should be \"26 C F\"", () => readingReads(page, "DEMOG of row 1", el("tile viewer"), "26 C F"));
      await session.step(119, "And the table should have 11 columns", () => columnCount(page, 11));
      await session.step(120, "When user remembers the fields of tile viewer", () => rememberFields(page, el("tile viewer")));
      await session.step(121, "And user picks \"Remove\" from the context menu of the \"field DEMOG of row 1\" area of tile viewer", () => pickFromAreaContextMenu(page, "Remove", "field DEMOG of row 1", el("tile viewer")));
      await session.step(122, "Then the table should not have a column \"DEMOG\"", () => hasNoColumn(page, "DEMOG"));
      await session.step(123, "And the table should have 10 columns", () => columnCount(page, 10));
      await session.step(124, "And the \"fields shown\" reading of tile viewer should be 9", () => readingIs(page, "fields shown", el("tile viewer"), 9));
      await session.step(125, "And the fields of tile viewer should be as remembered", () => fieldsAsRemembered(page, el("tile viewer")));
      await session.step(126, "And the \"fields\" reading of tile viewer should not contain \"SEVERITY\"", () => readingNotContains(page, "fields", el("tile viewer"), "SEVERITY"));
      await session.step(127, "And tile viewer should have a \"field DEMOG of row 1\" area", () => hasArea(page, el("tile viewer"), "field DEMOG of row 1"));
      await session.step(128, "And the \"DEMOG of row 1\" reading of tile viewer should be \"\"", () => readingReads(page, "DEMOG of row 1", el("tile viewer"), ""));
      await session.step(129, "And the \"auto generate\" reading of tile viewer should be \"false\"", () => readingReads(page, "auto generate", el("tile viewer"), "false"));
      await session.step(130, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Auto Generate rebuilds the card, and then a departing column's slot is refilled", async () => {
      await session.step(133, "Given user adds a calculated column \"DEMOG\" with formula \"${AGE}\"", () => addCalculated(page, "DEMOG", "${AGE}"));
      await session.step(134, "Then the table should have 11 columns", () => columnCount(page, 11));
      await session.step(135, "When user sets \"Auto Generate\" property of tile viewer to \"true\"", () => setProperty(page, "Auto Generate", el("tile viewer"), "true"));
      await session.step(136, "And user adds a calculated column \"TMP\" with formula \"1\"", () => addCalculated(page, "TMP", "1"));
      await session.step(137, "And user removes \"TMP\" column", () => removeColumn(page, "TMP"));
      await session.step(138, "Then the table should have 11 columns", () => columnCount(page, 11));
      await session.step(139, "And the \"auto generate\" reading of tile viewer should be \"true\"", () => readingReads(page, "auto generate", el("tile viewer"), "true"));
      await session.step(140, "And the \"form designed\" reading of tile viewer should be \"false\"", () => readingReads(page, "form designed", el("tile viewer"), "false"));
      await session.step(141, "And the \"fields shown\" reading of tile viewer should be 10", () => readingIs(page, "fields shown", el("tile viewer"), 10));
      await session.step(142, "And the \"fields\" reading of tile viewer should contain \"AGE\"", () => readingContains(page, "fields", el("tile viewer"), "AGE"));
      await session.step(143, "When user remembers the fields of tile viewer", () => rememberFields(page, el("tile viewer")));
      await session.step(144, "And user removes \"AGE\" column", () => removeColumn(page, "AGE"));
      await session.step(145, "Then the table should not have a column \"AGE\"", () => hasNoColumn(page, "AGE"));
      await session.step(146, "And the table should have 10 columns", () => columnCount(page, 10));
      await session.step(147, "And the \"fields shown\" reading of tile viewer should be 10", () => readingIs(page, "fields shown", el("tile viewer"), 10));
      await session.step(148, "And the fields of tile viewer should have refilled the freed slot", () => fieldsRefilled(page, el("tile viewer")));
      await session.step(149, "And the \"fields\" reading of tile viewer should not contain \"AGE\"", () => readingNotContains(page, "fields", el("tile viewer"), "AGE"));
      await session.step(150, "And tile viewer should not have a \"field AGE of row 1\" area", () => hasNoArea(page, el("tile viewer"), "field AGE of row 1"));
      await session.step(151, "And the \"auto generate\" reading of tile viewer should be \"true\"", () => readingReads(page, "auto generate", el("tile viewer"), "true"));
      await session.step(152, "And no errors should have been logged", () => noErrors(page));
      await session.step(153, "When user adds a calculated column \"AGE\" with formula \"${HEIGHT}\"", () => addCalculated(page, "AGE", "${HEIGHT}"));
      await session.step(154, "Then the table should have 11 columns", () => columnCount(page, 11));
    });
    run.finish();
  });
});
