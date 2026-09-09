/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/forms/forms-persistence.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.forms]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {selectWhereOneOf} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, openDataset, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, addViewerWith, loadLayout, noErrors, pickFromAreaContextMenu, propertyShouldBe, readingIs, readingReads, saveLayoutToServer, setProperties} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Forms viewer layout and project persistence", () => {
  const session = feature(test, "features/viewers/forms/forms-persistence.feature", import.meta.url);
  test("Forms viewer layout and project persistence", {tag: ["@journey", "@viewers", "@realizes:viewers.forms"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(12, "And user adds a forms viewer with:", () => addViewerWith(page, "forms", [["Show Mouse Over Row","false"]]));
    await session.step(14, "And user sets properties of forms viewer:", () => setProperties(page, el("forms viewer"), [["Fields","AGE, USUBJID, HEIGHT, SEX"],["Sort By","AGE"]]));
    await session.step(17, "And user selects rows where \"USUBJID\" is one of \"X0273T21000400001, X0273T21000500008, X0273T21001500015\"", () => selectWhereOneOf(page, "USUBJID", "X0273T21000400001, X0273T21000500008, X0273T21001500015"));
    await session.step(18, "And user picks \"Pin Row\" from the context menu of the \"field USUBJID of card 2\" area of forms viewer", () => pickFromAreaContextMenu(page, "Pin Row", "field USUBJID of card 2", el("forms viewer")));
    await session.step(19, "Then the \"fields\" reading of forms viewer should be \"AGE, USUBJID, HEIGHT, SEX\"", () => readingReads(page, "fields", el("forms viewer"), "AGE, USUBJID, HEIGHT, SEX"));
    await session.step(20, "And the \"sort column\" reading of forms viewer should be \"AGE\"", () => readingReads(page, "sort column", el("forms viewer"), "AGE"));
    await session.step(21, "And the \"pinned records\" reading of forms viewer should be 1", () => readingIs(page, "pinned records", el("forms viewer"), 1));
    await session.step(22, "And the \"USUBJID of pinned card 1\" reading of forms viewer should be \"X0273T21000400001\"", () => readingReads(page, "USUBJID of pinned card 1", el("forms viewer"), "X0273T21000400001"));
    await run.scenario("A saved layout restores the viewer over a changed view", async () => {
      await session.step(25, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(26, "And user clicks on close icon of forms viewer", () => clickOn(page, el("close icon of forms viewer")));
      await session.step(27, "Then forms viewer should be absent", () => shouldBe(page, el("forms viewer"), "absent"));
      await session.step(28, "When user adds a histogram viewer", () => addViewer(page, "histogram"));
      await session.step(29, "And user loads the saved layout", () => loadLayout(page));
      await session.step(30, "Then histogram viewer should be absent", () => shouldBe(page, el("histogram viewer"), "absent"));
      await session.step(31, "And forms viewer should be visible", () => shouldBe(page, el("forms viewer"), "visible"));
      await session.step(32, "And the \"fields\" reading of forms viewer should be \"AGE, USUBJID, HEIGHT, SEX\"", () => readingReads(page, "fields", el("forms viewer"), "AGE, USUBJID, HEIGHT, SEX"));
      await session.step(33, "And \"Fields\" property of forms viewer should be \"AGE, USUBJID, HEIGHT, SEX\"", () => propertyShouldBe(page, "Fields", el("forms viewer"), "AGE, USUBJID, HEIGHT, SEX"));
      await session.step(34, "And the \"sort column\" reading of forms viewer should be \"AGE\"", () => readingReads(page, "sort column", el("forms viewer"), "AGE"));
      await session.step(35, "And the \"pinned records\" reading of forms viewer should be 1", () => readingIs(page, "pinned records", el("forms viewer"), 1));
      await session.step(36, "And the \"USUBJID of pinned card 1\" reading of forms viewer should be \"X0273T21000400001\"", () => readingReads(page, "USUBJID of pinned card 1", el("forms viewer"), "X0273T21000400001"));
      await session.step(37, "And \"pinnedRowValues\" property of forms viewer should be \"X0273T21000400001\"", () => propertyShouldBe(page, "pinnedRowValues", el("forms viewer"), "X0273T21000400001"));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A project round-trip keeps the field set, the sort column and the pinned row", async () => {
      await session.step(41, "When user saves the current view as project \"zz-forms-persistence\"", () => saveAsProject(page, "zz-forms-persistence"));
      await session.step(42, "And user closes all views", () => closeAllViews(page));
      await session.step(43, "And user opens the \"zz-forms-persistence\" project", () => openProject(page, "zz-forms-persistence"));
      await session.step(44, "Then forms viewer should be visible", () => shouldBe(page, el("forms viewer"), "visible"));
      await session.step(45, "And the \"fields\" reading of forms viewer should be \"AGE, USUBJID, HEIGHT, SEX\"", () => readingReads(page, "fields", el("forms viewer"), "AGE, USUBJID, HEIGHT, SEX"));
      await session.step(46, "And the \"sort column\" reading of forms viewer should be \"AGE\"", () => readingReads(page, "sort column", el("forms viewer"), "AGE"));
      await session.step(47, "And the \"pinned records\" reading of forms viewer should be 1", () => readingIs(page, "pinned records", el("forms viewer"), 1));
      await session.step(48, "And the \"USUBJID of pinned card 1\" reading of forms viewer should be \"X0273T21000400001\"", () => readingReads(page, "USUBJID of pinned card 1", el("forms viewer"), "X0273T21000400001"));
      await session.step(49, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
