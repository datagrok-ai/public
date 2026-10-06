/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sar/sequence-space.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {peptidesInitialized, sarReady, sarSetting} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, enterInto, expand, shouldBe, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasColumn, hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {listenCustom} from '@datagrok-libraries/bdd/bindings/platform/events';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, painted, readingIs, viewerAdded, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Sequence space from the settings dialog", () => {
  const session = feature(test, "features/sar/sequence-space.feature", import.meta.url);
  test("Sequence space from the settings dialog", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(8, "Given user is logged in", () => loggedIn(page));
    await session.step(9, "And the Peptides package is initialized", () => peptidesInitialized(page));
    await session.step(10, "And user opens peptides dataset", () => openDataset(page, ds("peptides")));
    await session.step(11, "When user picks \"Bio > Analyze > SAR...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > SAR..."));
    await session.step(12, "Then \"Analyze Peptides\" dialog should be visible", () => shouldBe(page, el("\"Analyze Peptides\" dialog"), "visible"));
    await session.step(13, "And \"Generate clusters\" checkbox in \"Analyze Peptides\" dialog should be checked", () => shouldBe(page, el("\"Generate clusters\" checkbox in \"Analyze Peptides\" dialog"), "checked"));
    await session.step(14, "When user clicks on \"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog", () => clickOn(page, el("\"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog")));
    await session.step(15, "And user enters \"93\" into \"Similarity Threshold\" input in \"Analyze Peptides\" dialog", () => enterInto(page, "93", el("\"Similarity Threshold\" input in \"Analyze Peptides\" dialog")));
    await session.step(16, "Given user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
    await session.step(17, "When user clicks on OK button in \"Analyze Peptides\" dialog", () => clickOn(page, el("OK button in \"Analyze Peptides\" dialog")));
    await session.step(18, "Then the SAR analysis should be ready", () => sarReady(page));
    await session.step(19, "And MCL viewer should be added to the open tableview", () => viewerAdded(page, "MCL"));
    await session.step(20, "And the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
    await session.step(21, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(22, "And no errors should have been logged", () => noErrors(page));
    await run.scenario("Checking Sequence space adds the embedding scatter plot", async () => {
      await session.step(28, "When user clicks on \"Peptides analysis settings\" icon", () => clickOn(page, el("\"Peptides analysis settings\" icon")));
      await session.step(29, "And user expands Viewers pane in \"Peptides settings\" dialog", () => expand(page, el("Viewers pane in \"Peptides settings\" dialog")));
      await session.step(30, "Then \"Sequence space\" checkbox in \"Peptides settings\" dialog should be unchecked", () => shouldBe(page, el("\"Sequence space\" checkbox in \"Peptides settings\" dialog"), "unchecked"));
      await session.step(31, "When user checks \"Sequence space\" checkbox in \"Peptides settings\" dialog", () => check(page, el("\"Sequence space\" checkbox in \"Peptides settings\" dialog")));
      await session.step(32, "Given user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
      await session.step(33, "When user clicks on OK button in \"Peptides settings\" dialog", () => clickOn(page, el("OK button in \"Peptides settings\" dialog")));
      await session.step(34, "Then \"Peptides settings\" dialog should be hidden", () => shouldBe(page, el("\"Peptides settings\" dialog"), "hidden"));
      await session.step(35, "And the SAR analysis should be ready", () => sarReady(page));
      await session.step(36, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(37, "And the SAR setting \"showSequenceSpace\" should be \"true\"", () => sarSetting(page, "showSequenceSpace", "true"));
      await session.step(38, "And scatter plot viewer should be added to the open tableview", () => viewerAdded(page, "scatter plot"));
      await session.step(39, "And the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
      await session.step(40, "And the table should have a column \"Embed_X_1\"", () => hasColumn(page, "Embed_X_1"));
      await session.step(41, "And the table should have a column \"Embed_Y_1\"", () => hasColumn(page, "Embed_Y_1"));
      await session.step(42, "And the table should have a column \"Cluster (DBSCAN)\"", () => hasColumn(page, "Cluster (DBSCAN)"));
      await session.step(43, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
      await session.step(44, "And the \"rows shown\" reading of scatter plot viewer should be 647", () => readingIs(page, "rows shown", el("scatter plot viewer"), 647));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Unchecking Sequence space removes the scatter plot and its columns", async () => {
      await session.step(48, "When user clicks on \"Peptides analysis settings\" icon", () => clickOn(page, el("\"Peptides analysis settings\" icon")));
      await session.step(49, "And user expands Viewers pane in \"Peptides settings\" dialog", () => expand(page, el("Viewers pane in \"Peptides settings\" dialog")));
      await session.step(50, "Then \"Sequence space\" checkbox in \"Peptides settings\" dialog should be checked", () => shouldBe(page, el("\"Sequence space\" checkbox in \"Peptides settings\" dialog"), "checked"));
      await session.step(51, "When user unchecks \"Sequence space\" checkbox in \"Peptides settings\" dialog", () => uncheck(page, el("\"Sequence space\" checkbox in \"Peptides settings\" dialog")));
      await session.step(52, "Given user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
      await session.step(53, "When user clicks on OK button in \"Peptides settings\" dialog", () => clickOn(page, el("OK button in \"Peptides settings\" dialog")));
      await session.step(54, "Then \"Peptides settings\" dialog should be hidden", () => shouldBe(page, el("\"Peptides settings\" dialog"), "hidden"));
      await session.step(55, "And the SAR analysis should be ready", () => sarReady(page));
      await session.step(56, "And the SAR setting \"showSequenceSpace\" should be \"false\"", () => sarSetting(page, "showSequenceSpace", "false"));
      await session.step(57, "And the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
      await session.step(58, "And the table should not have a column \"Embed_X_1\"", () => hasNoColumn(page, "Embed_X_1"));
      await session.step(59, "And the table should not have a column \"Cluster (DBSCAN)\"", () => hasNoColumn(page, "Cluster (DBSCAN)"));
      await session.step(60, "And the open tableview should have 1 MCL viewer", () => viewerCount(page, 1, "MCL"));
      await session.step(61, "And no errors should have been logged", () => noErrors(page));
      await session.step(62, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Checking Sequence space again adds one scatter plot", async () => {
      await session.step(65, "When user clicks on \"Peptides analysis settings\" icon", () => clickOn(page, el("\"Peptides analysis settings\" icon")));
      await session.step(66, "And user expands Viewers pane in \"Peptides settings\" dialog", () => expand(page, el("Viewers pane in \"Peptides settings\" dialog")));
      await session.step(67, "And user checks \"Sequence space\" checkbox in \"Peptides settings\" dialog", () => check(page, el("\"Sequence space\" checkbox in \"Peptides settings\" dialog")));
      await session.step(68, "Given user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
      await session.step(69, "When user clicks on OK button in \"Peptides settings\" dialog", () => clickOn(page, el("OK button in \"Peptides settings\" dialog")));
      await session.step(70, "Then the SAR analysis should be ready", () => sarReady(page));
      await session.step(71, "And the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
      await session.step(72, "And the table should have a column \"Embed_X_1\"", () => hasColumn(page, "Embed_X_1"));
      await session.step(73, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
      await session.step(74, "And no errors should have been logged", () => noErrors(page));
      await session.step(75, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
