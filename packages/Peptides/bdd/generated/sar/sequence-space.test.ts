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
import {check, clickOn, enterInto, expand, selectIn, shouldBe, shouldHaveValue, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
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
    const run = journey(test, 4, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And the Peptides package is initialized", () => peptidesInitialized(page));
    await session.step(13, "And user opens peptides dataset", () => openDataset(page, ds("peptides")));
    await session.step(14, "When user picks \"Bio > Analyze > SAR...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > SAR..."));
    await session.step(15, "Then \"Analyze Peptides\" dialog should be visible", () => shouldBe(page, el("\"Analyze Peptides\" dialog"), "visible"));
    await session.step(16, "And \"Generate clusters\" checkbox in \"Analyze Peptides\" dialog should be checked", () => shouldBe(page, el("\"Generate clusters\" checkbox in \"Analyze Peptides\" dialog"), "checked"));
    await session.step(17, "When user clicks on \"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog", () => clickOn(page, el("\"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog")));
    await session.step(18, "And user enters \"93\" into \"Similarity Threshold\" input in \"Analyze Peptides\" dialog", () => enterInto(page, "93", el("\"Similarity Threshold\" input in \"Analyze Peptides\" dialog")));
    await session.step(19, "Given user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
    await session.step(20, "When user clicks on OK button in \"Analyze Peptides\" dialog", () => clickOn(page, el("OK button in \"Analyze Peptides\" dialog")));
    await session.step(21, "Then the SAR analysis should be ready", () => sarReady(page));
    await session.step(22, "And MCL viewer should be added to the open tableview", () => viewerAdded(page, "MCL"));
    await session.step(23, "And the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
    await session.step(24, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(25, "And no errors should have been logged", () => noErrors(page));
    await run.scenario("Checking Sequence space adds the embedding scatter plot", async () => {
      await session.step(31, "When user clicks on \"Peptides analysis settings\" icon", () => clickOn(page, el("\"Peptides analysis settings\" icon")));
      await session.step(32, "And user expands Viewers pane in \"Peptides settings\" dialog", () => expand(page, el("Viewers pane in \"Peptides settings\" dialog")));
      await session.step(33, "Then \"Sequence space\" checkbox in \"Peptides settings\" dialog should be unchecked", () => shouldBe(page, el("\"Sequence space\" checkbox in \"Peptides settings\" dialog"), "unchecked"));
      await session.step(34, "When user checks \"Sequence space\" checkbox in \"Peptides settings\" dialog", () => check(page, el("\"Sequence space\" checkbox in \"Peptides settings\" dialog")));
      await session.step(35, "Given user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
      await session.step(36, "When user clicks on OK button in \"Peptides settings\" dialog", () => clickOn(page, el("OK button in \"Peptides settings\" dialog")));
      await session.step(37, "Then \"Peptides settings\" dialog should be hidden", () => shouldBe(page, el("\"Peptides settings\" dialog"), "hidden"));
      await session.step(38, "And the SAR analysis should be ready", () => sarReady(page));
      await session.step(39, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(40, "And the SAR setting \"showSequenceSpace\" should be \"true\"", () => sarSetting(page, "showSequenceSpace", "true"));
      await session.step(41, "And scatter plot viewer should be added to the open tableview", () => viewerAdded(page, "scatter plot"));
      await session.step(42, "And the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
      await session.step(43, "And the table should have a column \"Embed_X_1\"", () => hasColumn(page, "Embed_X_1"));
      await session.step(44, "And the table should have a column \"Embed_Y_1\"", () => hasColumn(page, "Embed_Y_1"));
      await session.step(45, "And the table should have a column \"Cluster (DBSCAN)\"", () => hasColumn(page, "Cluster (DBSCAN)"));
      await session.step(46, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
      await session.step(47, "And the \"rows shown\" reading of scatter plot viewer should be 647", () => readingIs(page, "rows shown", el("scatter plot viewer"), 647));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Unchecking Sequence space removes the scatter plot and its columns", async () => {
      await session.step(51, "When user clicks on \"Peptides analysis settings\" icon", () => clickOn(page, el("\"Peptides analysis settings\" icon")));
      await session.step(52, "And user expands Viewers pane in \"Peptides settings\" dialog", () => expand(page, el("Viewers pane in \"Peptides settings\" dialog")));
      await session.step(53, "Then \"Sequence space\" checkbox in \"Peptides settings\" dialog should be checked", () => shouldBe(page, el("\"Sequence space\" checkbox in \"Peptides settings\" dialog"), "checked"));
      await session.step(54, "When user unchecks \"Sequence space\" checkbox in \"Peptides settings\" dialog", () => uncheck(page, el("\"Sequence space\" checkbox in \"Peptides settings\" dialog")));
      await session.step(55, "Given user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
      await session.step(56, "When user clicks on OK button in \"Peptides settings\" dialog", () => clickOn(page, el("OK button in \"Peptides settings\" dialog")));
      await session.step(57, "Then \"Peptides settings\" dialog should be hidden", () => shouldBe(page, el("\"Peptides settings\" dialog"), "hidden"));
      await session.step(58, "And the SAR analysis should be ready", () => sarReady(page));
      await session.step(59, "And the SAR setting \"showSequenceSpace\" should be \"false\"", () => sarSetting(page, "showSequenceSpace", "false"));
      await session.step(60, "And the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
      await session.step(61, "And the table should not have a column \"Embed_X_1\"", () => hasNoColumn(page, "Embed_X_1"));
      await session.step(62, "And the table should not have a column \"Cluster (DBSCAN)\"", () => hasNoColumn(page, "Cluster (DBSCAN)"));
      await session.step(63, "And the open tableview should have 1 MCL viewer", () => viewerCount(page, 1, "MCL"));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
      await session.step(65, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Checking Sequence space again adds one scatter plot", async () => {
      await session.step(68, "When user clicks on \"Peptides analysis settings\" icon", () => clickOn(page, el("\"Peptides analysis settings\" icon")));
      await session.step(69, "And user expands Viewers pane in \"Peptides settings\" dialog", () => expand(page, el("Viewers pane in \"Peptides settings\" dialog")));
      await session.step(70, "And user checks \"Sequence space\" checkbox in \"Peptides settings\" dialog", () => check(page, el("\"Sequence space\" checkbox in \"Peptides settings\" dialog")));
      await session.step(71, "Given user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
      await session.step(72, "When user clicks on OK button in \"Peptides settings\" dialog", () => clickOn(page, el("OK button in \"Peptides settings\" dialog")));
      await session.step(73, "Then the SAR analysis should be ready", () => sarReady(page));
      await session.step(74, "And the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
      await session.step(75, "And the table should have a column \"Embed_X_1\"", () => hasColumn(page, "Embed_X_1"));
      await session.step(76, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
      await session.step(77, "And no errors should have been logged", () => noErrors(page));
      await session.step(78, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The Sequence space pane shows the inputs of the chosen distance, and a new distance embeds the peptides again", async () => {
      await session.step(81, "When user clicks on \"Peptides analysis settings\" icon", () => clickOn(page, el("\"Peptides analysis settings\" icon")));
      await session.step(82, "And user expands \"Sequence space\" pane in \"Peptides settings\" dialog", () => expand(page, el("\"Sequence space\" pane in \"Peptides settings\" dialog")));
      await session.step(83, "Then \"Distance Function\" input in \"Sequence space\" pane in \"Peptides settings\" dialog should have value \"Needlemann-Wunsch\"", () => shouldHaveValue(page, el("\"Distance Function\" input in \"Sequence space\" pane in \"Peptides settings\" dialog"), "Needlemann-Wunsch"));
      await session.step(84, "And \"Gap Open Penalty\" input in \"Sequence space\" pane in \"Peptides settings\" dialog should be visible", () => shouldBe(page, el("\"Gap Open Penalty\" input in \"Sequence space\" pane in \"Peptides settings\" dialog"), "visible"));
      await session.step(85, "And \"Fingerprint Type\" input in \"Sequence space\" pane in \"Peptides settings\" dialog should be visible", () => shouldBe(page, el("\"Fingerprint Type\" input in \"Sequence space\" pane in \"Peptides settings\" dialog"), "visible"));
      await session.step(86, "When user selects \"Hamming\" in \"Distance Function\" input in \"Sequence space\" pane in \"Peptides settings\" dialog", () => selectIn(page, "Hamming", el("\"Distance Function\" input in \"Sequence space\" pane in \"Peptides settings\" dialog")));
      await session.step(87, "Then \"Gap Open Penalty\" input in \"Sequence space\" pane in \"Peptides settings\" dialog should be hidden", () => shouldBe(page, el("\"Gap Open Penalty\" input in \"Sequence space\" pane in \"Peptides settings\" dialog"), "hidden"));
      await session.step(88, "And \"Fingerprint Type\" input in \"Sequence space\" pane in \"Peptides settings\" dialog should be hidden", () => shouldBe(page, el("\"Fingerprint Type\" input in \"Sequence space\" pane in \"Peptides settings\" dialog"), "hidden"));
      await session.step(89, "Given user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
      await session.step(90, "When user clicks on OK button in \"Peptides settings\" dialog", () => clickOn(page, el("OK button in \"Peptides settings\" dialog")));
      await session.step(91, "Then \"Peptides settings\" dialog should be hidden", () => shouldBe(page, el("\"Peptides settings\" dialog"), "hidden"));
      await session.step(92, "And the SAR analysis should be ready", () => sarReady(page));
      await session.step(93, "And the SAR setting \"sequenceSpaceParams.distanceF\" should be \"Hamming\"", () => sarSetting(page, "sequenceSpaceParams.distanceF", "Hamming"));
      await session.step(94, "And the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
      await session.step(95, "And the table should have a column \"Embed_X_1\"", () => hasColumn(page, "Embed_X_1"));
      await session.step(96, "And the \"rows shown\" reading of scatter plot viewer should be 647", () => readingIs(page, "rows shown", el("scatter plot viewer"), 647));
      await session.step(97, "And no errors should have been logged", () => noErrors(page));
      await session.step(98, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
