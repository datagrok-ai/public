/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/navigation/direct-link-loading.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/add-new-column.js';
import '../../bindings/enrichment.js';
import '../../bindings/formula-lines.js';
import '../../bindings/home.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {refreshBrowseTree} from '../../bindings/io.js';
import {noLoader, openDirectLink, projectIsOpen, projectNotOpen, viewNotCurrent} from '../../bindings/navigation.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, expand, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeAllViews, noProjectOnServer, openDataset, projectsOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A project opened by its direct link", () => {
  const session = feature(test, "features/navigation/direct-link-loading.feature", import.meta.url);
  test("A project opened by its direct link", {tag: ["@journey", "@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await run.scenario("A project with the demog table is saved from the toolbar", async () => {
      await session.step(23, "Given user opens demog dataset", () => openDataset(page, ds("demog")));
      await session.step(24, "And no project named \"bdd-direct-link-{time}\" is on the server", () => noProjectOnServer(page, session.text("bdd-direct-link-{time}")));
      await session.step(25, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(26, "When user clicks on Save button in toolbar", () => clickOn(page, el("Save button in toolbar")));
      await session.step(27, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(28, "When user enters \"bdd-direct-link-{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("bdd-direct-link-{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(29, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(30, "Then \"Save project\" dialog should be hidden", () => shouldBe(page, el("\"Save project\" dialog"), "hidden"));
      await session.step(31, "And 1 project named \"bdd-direct-link-{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("bdd-direct-link-{time}")));
      await session.step(32, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The direct link opens the project's table view", async () => {
      await session.step(35, "When user closes all views", () => closeAllViews(page));
      await session.step(36, "Then the project \"bdd-direct-link-{time}\" should not be open", () => projectNotOpen(page, session.text("bdd-direct-link-{time}")));
      await session.step(37, "When user loads the direct link of project \"bdd-direct-link-{time}\"", () => openDirectLink(page, session.text("bdd-direct-link-{time}")));
      await session.step(38, "Then the project \"bdd-direct-link-{time}\" should be open", () => projectIsOpen(page, session.text("bdd-direct-link-{time}")));
      await session.step(39, "And the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(40, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(41, "And the table should have a column \"AGE\"", () => hasColumn(page, "AGE"));
      await session.step(42, "And the table should have a column \"RACE\"", () => hasColumn(page, "RACE"));
      await session.step(43, "And grid should show 5850 rows", () => showsRows(page, el("grid"), 5850));
      await session.step(44, "And the \"Home\" view should not be current", () => viewNotCurrent(page, "Home"));
      await session.step(45, "And no loading indicator should be visible", () => noLoader(page));
      await session.step(46, "And no errors should have been logged", () => noErrors(page));
      await session.step(47, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The same project opened from the Browse tree shows the same view", async () => {
      await session.step(50, "When user closes all views", () => closeAllViews(page));
      await session.step(51, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(52, "And the project \"bdd-direct-link-{time}\" should not be open", () => projectNotOpen(page, session.text("bdd-direct-link-{time}")));
      await session.step(53, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(54, "When user refreshes the browse tree", () => refreshBrowseTree(page));
      await session.step(55, "And user expands \"My stuff\" tree node inside browse tree", () => expand(page, el("\"My stuff\" tree node inside browse tree")));
      await session.step(56, "And user double-clicks on \"My stuff > bdd-direct-link-{time}\" tree node inside browse tree", () => doubleClickOn(page, el(session.text("\"My stuff > bdd-direct-link-{time}\" tree node inside browse tree"))));
      await session.step(57, "Then the project \"bdd-direct-link-{time}\" should be open", () => projectIsOpen(page, session.text("bdd-direct-link-{time}")));
      await session.step(58, "And the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(59, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(60, "And the table should have a column \"AGE\"", () => hasColumn(page, "AGE"));
      await session.step(61, "And grid should show 5850 rows", () => showsRows(page, el("grid"), 5850));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
      await session.step(63, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
