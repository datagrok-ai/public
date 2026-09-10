/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/forms/forms-layout.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.forms]
--- */
import {test} from '@playwright/test';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {pinnedCardRows, recordCardRows} from '../../../bindings/forms.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {selectWhereIs} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, openDataset, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, addViewerWith, hasArea, loadLayout, noBalloons, noErrors, pickFromAreaContextMenu, readingReads, readingsEqual, saveLayoutToServer, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {fieldsAsRemembered, rememberFields} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Forms viewer layout and project round-trips", () => {
  const session = feature(test, "features/viewers/forms/forms-layout.feature", import.meta.url);
  test("Forms viewer layout and project round-trips", {tag: ["@journey", "@viewers", "@realizes:viewers.forms"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a forms viewer with:", () => addViewerWith(page, "forms", [["fieldsColumnNames","SEVERITY, STARTED, CONTROL, DEMOG, WEIGHT, HEIGHT, DIS_POP, RACE, SEX, AGE, USUBJID"],["sortByColumnName","AGE"],["showMouseOverRow","false"]]));
    await session.step(21, "And user selects rows where \"SEVERITY\" is \"Critical\"", () => selectWhereIs(page, "SEVERITY", "Critical"));
    await session.step(22, "Then forms viewer should be visible", () => shouldBe(page, el("forms viewer"), "visible"));
    await session.step(23, "And the \"fields\" reading of forms viewer should be \"SEVERITY, STARTED, CONTROL, DEMOG, WEIGHT, HEIGHT, DIS_POP, RACE, SEX, AGE, USUBJID\"", () => readingReads(page, "fields", el("forms viewer"), "SEVERITY, STARTED, CONTROL, DEMOG, WEIGHT, HEIGHT, DIS_POP, RACE, SEX, AGE, USUBJID"));
    await session.step(24, "And the \"sort column\" reading of forms viewer should be \"AGE\"", () => readingReads(page, "sort column", el("forms viewer"), "AGE"));
    await session.step(25, "And the record cards of forms viewer should show rows \"304, 512, 428, 430, 215\"", () => recordCardRows(page, el("forms viewer"), "304, 512, 428, 430, 215"));
    await run.scenario("The field set is drawn in the order it was given, not in table order", async () => {
      await session.step(28, "Then the \"header labels\" reading of forms viewer should be \"SEVERITY, STARTED, CONTROL, DEMOG, WEIGHT, HEIGHT, DIS_POP, RACE, SEX, AGE, USUBJID\"", () => readingReads(page, "header labels", el("forms viewer"), "SEVERITY, STARTED, CONTROL, DEMOG, WEIGHT, HEIGHT, DIS_POP, RACE, SEX, AGE, USUBJID"));
      await session.step(29, "And the \"fields\" and \"header labels\" readings of forms viewer should be the same", () => readingsEqual(page, "fields", "header labels", el("forms viewer")));
      await session.step(30, "And forms viewer should have a \"sort indicator AGE\" area", () => hasArea(page, el("forms viewer"), "sort indicator AGE"));
      await session.step(31, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A row is pinned by value, and the pinned pane holds it", async () => {
      await session.step(34, "Then the \"pinned pane shown\" reading of forms viewer should be \"false\"", () => readingReads(page, "pinned pane shown", el("forms viewer"), "false"));
      await session.step(35, "When user picks \"Pin Row\" from the context menu of the \"field USUBJID of card 2\" area of forms viewer", () => pickFromAreaContextMenu(page, "Pin Row", "field USUBJID of card 2", el("forms viewer")));
      await session.step(36, "Then the \"pinned pane shown\" reading of forms viewer should be \"true\"", () => readingReads(page, "pinned pane shown", el("forms viewer"), "true"));
      await session.step(37, "And the \"pinned by\" reading of forms viewer should be \"USUBJID\"", () => readingReads(page, "pinned by", el("forms viewer"), "USUBJID"));
      await session.step(38, "And the \"pinned values\" reading of forms viewer should be \"X0273T29012500105\"", () => readingReads(page, "pinned values", el("forms viewer"), "X0273T29012500105"));
      await session.step(39, "And the pinned cards of forms viewer should show rows \"304\"", () => pinnedCardRows(page, el("forms viewer"), "304"));
      await session.step(40, "And the record cards of forms viewer should show rows \"512, 428, 430, 215\"", () => recordCardRows(page, el("forms viewer"), "512, 428, 430, 215"));
      await session.step(41, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Re-applying the saved layout over a corrupted view restores the viewer and drops a foreign one", async () => {
      await session.step(45, "Given user remembers the fields of forms viewer", () => rememberFields(page, el("forms viewer")));
      await session.step(46, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(47, "And user clicks on close icon of forms viewer", () => clickOn(page, el("close icon of forms viewer")));
      await session.step(48, "Then forms viewer should be absent", () => shouldBe(page, el("forms viewer"), "absent"));
      await session.step(49, "When user adds a histogram viewer", () => addViewer(page, "histogram"));
      await session.step(50, "Then histogram viewer should be visible", () => shouldBe(page, el("histogram viewer"), "visible"));
      await session.step(51, "When user loads the saved layout", () => loadLayout(page));
      await session.step(52, "Then forms viewer should be visible", () => shouldBe(page, el("forms viewer"), "visible"));
      await session.step(53, "And histogram viewer should be absent", () => shouldBe(page, el("histogram viewer"), "absent"));
      await session.step(54, "And the open tableview should have 0 histogram viewers", () => viewerCount(page, 0, "histogram"));
      await session.step(55, "And the fields of forms viewer should be as remembered", () => fieldsAsRemembered(page, el("forms viewer")));
      await session.step(56, "And the \"fields\" reading of forms viewer should be \"SEVERITY, STARTED, CONTROL, DEMOG, WEIGHT, HEIGHT, DIS_POP, RACE, SEX, AGE, USUBJID\"", () => readingReads(page, "fields", el("forms viewer"), "SEVERITY, STARTED, CONTROL, DEMOG, WEIGHT, HEIGHT, DIS_POP, RACE, SEX, AGE, USUBJID"));
      await session.step(57, "And the \"header labels\" reading of forms viewer should be \"SEVERITY, STARTED, CONTROL, DEMOG, WEIGHT, HEIGHT, DIS_POP, RACE, SEX, AGE, USUBJID\"", () => readingReads(page, "header labels", el("forms viewer"), "SEVERITY, STARTED, CONTROL, DEMOG, WEIGHT, HEIGHT, DIS_POP, RACE, SEX, AGE, USUBJID"));
      await session.step(58, "And the \"sort column\" reading of forms viewer should be \"AGE\"", () => readingReads(page, "sort column", el("forms viewer"), "AGE"));
      await session.step(59, "And forms viewer should have a \"sort indicator AGE\" area", () => hasArea(page, el("forms viewer"), "sort indicator AGE"));
      await session.step(60, "And the \"pinned by\" reading of forms viewer should be \"USUBJID\"", () => readingReads(page, "pinned by", el("forms viewer"), "USUBJID"));
      await session.step(61, "And the \"pinned values\" reading of forms viewer should be \"X0273T29012500105\"", () => readingReads(page, "pinned values", el("forms viewer"), "X0273T29012500105"));
      await session.step(62, "And the \"pinned pane shown\" reading of forms viewer should be \"true\"", () => readingReads(page, "pinned pane shown", el("forms viewer"), "true"));
      await session.step(63, "And the pinned cards of forms viewer should show rows \"304\"", () => pinnedCardRows(page, el("forms viewer"), "304"));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A project round-trip brings the field set and the pinned row back across a session", async () => {
      await session.step(67, "When user saves the current view as project \"bdd forms layout\"", () => saveAsProject(page, "bdd forms layout"));
      await session.step(68, "And user closes all views", () => closeAllViews(page));
      await session.step(69, "And user opens the \"bdd forms layout\" project", () => openProject(page, "bdd forms layout"));
      await session.step(70, "Then forms viewer should be visible", () => shouldBe(page, el("forms viewer"), "visible"));
      await session.step(71, "And the \"fields\" reading of forms viewer should be \"SEVERITY, STARTED, CONTROL, DEMOG, WEIGHT, HEIGHT, DIS_POP, RACE, SEX, AGE, USUBJID\"", () => readingReads(page, "fields", el("forms viewer"), "SEVERITY, STARTED, CONTROL, DEMOG, WEIGHT, HEIGHT, DIS_POP, RACE, SEX, AGE, USUBJID"));
      await session.step(72, "And the \"header labels\" reading of forms viewer should be \"SEVERITY, STARTED, CONTROL, DEMOG, WEIGHT, HEIGHT, DIS_POP, RACE, SEX, AGE, USUBJID\"", () => readingReads(page, "header labels", el("forms viewer"), "SEVERITY, STARTED, CONTROL, DEMOG, WEIGHT, HEIGHT, DIS_POP, RACE, SEX, AGE, USUBJID"));
      await session.step(73, "And the \"sort column\" reading of forms viewer should be \"AGE\"", () => readingReads(page, "sort column", el("forms viewer"), "AGE"));
      await session.step(74, "And forms viewer should have a \"sort indicator AGE\" area", () => hasArea(page, el("forms viewer"), "sort indicator AGE"));
      await session.step(75, "And the \"pinned by\" reading of forms viewer should be \"USUBJID\"", () => readingReads(page, "pinned by", el("forms viewer"), "USUBJID"));
      await session.step(76, "And the \"pinned values\" reading of forms viewer should be \"X0273T29012500105\"", () => readingReads(page, "pinned values", el("forms viewer"), "X0273T29012500105"));
      await session.step(77, "And the \"pinned pane shown\" reading of forms viewer should be \"true\"", () => readingReads(page, "pinned pane shown", el("forms viewer"), "true"));
      await session.step(78, "And the pinned cards of forms viewer should show rows \"304\"", () => pinnedCardRows(page, el("forms viewer"), "304"));
      await session.step(79, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
