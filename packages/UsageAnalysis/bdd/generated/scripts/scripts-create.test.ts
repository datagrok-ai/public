/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/scripts/scripts-create.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.scripts]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {ribbonReady, scriptResult} from '../../bindings/scripts.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {alertShown, check, clickOn, enterInto, expand, recordAlerts, selectIn, shouldBe, shouldContainText, shouldHaveText, typeInto, visibleCount} from '@datagrok-libraries/bdd/bindings/common/steps';
import {tableOpen, tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeCurrentView, dialogCloses, noScriptOnServer, saveScript, scriptsOnServer, switchView, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, menuLists, noBalloons, noErrors, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Creating a script", () => {
  const session = feature(test, "features/scripts/scripts-create.feature", import.meta.url);
  test("Creating a script", {tag: ["@journey", "@serial", "@realizes:views.scripts"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(31, "Given user is logged in", () => loggedIn(page));
    await session.step(32, "And no script named \"BddScriptCreate{time}\" is on the server", () => noScriptOnServer(page, session.text("BddScriptCreate{time}")));
    await session.step(33, "And the browse panel is open", () => browsePanelOpen(page));
    await run.scenario("The Scripts view opens from Platform > Functions", async () => {
      await session.step(36, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
      await session.step(37, "And user expands \"Platform > Functions\" tree node inside browse tree", () => expand(page, el("\"Platform > Functions\" tree node inside browse tree")));
      await session.step(38, "And user clicks on \"Platform > Functions > Scripts\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Functions > Scripts\" tree node inside browse tree")));
      await session.step(39, "Then the \"Scripts\" view should be current", () => viewIsCurrent(page, "Scripts"));
      await session.step(40, "And gallery should be visible", () => shouldBe(page, el("gallery"), "visible"));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("New offers every language and opens a template unsaved", async () => {
      await session.step(44, "When user clicks on New button", () => clickOn(page, el("New button")));
      await session.step(45, "Then the open menu should list \"R Script...\"", () => menuLists(page, "R Script..."));
      await session.step(46, "And the open menu should list \"Python Script...\"", () => menuLists(page, "Python Script..."));
      await session.step(47, "And the open menu should list \"Octave Script...\"", () => menuLists(page, "Octave Script..."));
      await session.step(48, "And the open menu should list \"NodeJS Script...\"", () => menuLists(page, "NodeJS Script..."));
      await session.step(49, "And the open menu should list \"Julia Script...\"", () => menuLists(page, "Julia Script..."));
      await session.step(50, "And the open menu should list \"JavaScript Script...\"", () => menuLists(page, "JavaScript Script..."));
      await session.step(51, "And the open menu should list \"Grok Script...\"", () => menuLists(page, "Grok Script..."));
      await session.step(52, "And the open menu should list \"Pyodide Script...\"", () => menuLists(page, "Pyodide Script..."));
      await session.step(53, "When user picks \"Grok Script...\" from the open menu", () => pickFromOpenMenu(page, "Grok Script..."));
      await session.step(54, "Then the \"Template\" view should be current", () => viewIsCurrent(page, "Template"));
      await session.step(55, "And code editor should contain the text \"#language: grok\"", () => shouldContainText(page, el("code editor"), "#language: grok"));
      await session.step(56, "And code editor should be visible", () => shouldBe(page, el("code editor"), "visible"));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The sample icon opens the sample table", async () => {
      await session.step(60, "When user clicks on \"Open script sample table\" icon", () => clickOn(page, el("\"Open script sample table\" icon")));
      await session.step(61, "Then table \"cars\" should be open", () => tableOpen(page, "cars"));
      await session.step(62, "And table \"cars\" should have 30 rows", () => tableRows(page, "cars", 30));
      await session.step(63, "And the \"Template\" view should be current", () => viewIsCurrent(page, "Template"));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Signature editor names the script and adds a parameter to its header", async () => {
      await session.step(67, "Given user switches to the \"Template\" view", () => switchView(page, "Template"));
      await session.step(71, "And the ribbon of the current view should be ready", () => ribbonReady(page));
      await session.step(72, "When user clicks on \"Open Signature Editor\" icon", () => clickOn(page, el("\"Open Signature Editor\" icon")));
      await session.step(73, "Then PARAMETERS tab should be visible", () => shouldBe(page, el("PARAMETERS tab"), "visible"));
      await session.step(74, "When user enters \"BddScriptCreate{time}\" into Name input", () => enterInto(page, session.text("BddScriptCreate{time}"), el("Name input")));
      await session.step(75, "And user clicks on PARAMETERS tab", () => clickOn(page, el("PARAMETERS tab")));
      await session.step(76, "Then there should be 2 visible \"Add the param\" icon", () => visibleCount(page, 2, el("\"Add the param\" icon")));
      await session.step(77, "When user clicks on first \"Add the param\" icon", () => clickOn(page, el("first \"Add the param\" icon")));
      await session.step(78, "Then there should be 3 visible \"Add the param\" icon", () => visibleCount(page, 3, el("\"Add the param\" icon")));
      await session.step(79, "Then \"Open function editor\" icon should be visible", () => shouldBe(page, el("\"Open function editor\" icon"), "visible"));
      await session.step(80, "When user clicks on \"Open function editor\" icon", () => clickOn(page, el("\"Open function editor\" icon")));
      await session.step(84, "Then \"Open function editor\" icon should be hidden", () => shouldBe(page, el("\"Open function editor\" icon"), "hidden"));
      await session.step(85, "And the \"BddScriptCreate{time}\" view should be current", () => viewIsCurrent(page, session.text("BddScriptCreate{time}")));
      await session.step(86, "And code editor should be visible", () => shouldBe(page, el("code editor"), "visible"));
      await session.step(93, "And the ribbon of the current view should be ready", () => ribbonReady(page));
      await session.step(94, "And code editor should contain the text \"#name: BddScriptCreate{time}\"", () => shouldContainText(page, el("code editor"), session.text("#name: BddScriptCreate{time}")));
      await session.step(95, "And code editor should contain the text \"#input: bool newParam\"", () => shouldContainText(page, el("code editor"), "#input: bool newParam"));
      await session.step(96, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Running with cars answers the number of cells", async () => {
      await session.step(100, "Given user switches to the \"BddScriptCreate{time}\" view", () => switchView(page, session.text("BddScriptCreate{time}")));
      await session.step(101, "When user clicks on \"Run script (F5)\" icon", () => clickOn(page, el("\"Run script (F5)\" icon")));
      await session.step(102, "Then \"BddScriptCreate{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"BddScriptCreate{time}\" dialog")), "visible"));
      await session.step(103, "When user selects \"cars\" in Table input in \"BddScriptCreate{time}\" dialog", () => selectIn(page, "cars", el(session.text("Table input in \"BddScriptCreate{time}\" dialog"))));
      await session.step(106, "And user checks \"New Param\" input in \"BddScriptCreate{time}\" dialog", () => check(page, el(session.text("\"New Param\" input in \"BddScriptCreate{time}\" dialog"))));
      await session.step(107, "And user clicks on OK button in \"BddScriptCreate{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"BddScriptCreate{time}\" dialog"))));
      await session.step(108, "Then the \"BddScriptCreate{time}\" dialog should close", () => dialogCloses(page, session.text("BddScriptCreate{time}")));
      await session.step(109, "And the script results should show \"count\" as \"510\"", () => scriptResult(page, "count", "510"));
      await session.step(110, "And no errors should have been logged", () => noErrors(page));
      await session.step(111, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Save stores the script", async () => {
      await session.step(114, "When user saves the script", () => saveScript(page));
      await session.step(115, "Then an info balloon containing \"Script saved.\" should have been shown", () => infoBalloonText(page, "Script saved."));
      await session.step(116, "And 1 script named \"BddScriptCreate{time}\" should be on the server", () => scriptsOnServer(page, 1, session.text("BddScriptCreate{time}")));
      await session.step(117, "And the \"BddScriptCreate{time}\" view should be current", () => viewIsCurrent(page, session.text("BddScriptCreate{time}")));
      await session.step(118, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Closing the editor returns to Scripts, where the new script is found", async () => {
      await session.step(121, "When user closes the current view", () => closeCurrentView(page));
      await session.step(122, "Then the \"Scripts\" view should be current", () => viewIsCurrent(page, "Scripts"));
      await session.step(123, "When user types \"BddScriptCreate{time}\" into gallery search", () => typeInto(page, session.text("BddScriptCreate{time}"), el("gallery search")));
      await session.step(124, "Then gallery counter should have text \"1\"", () => shouldHaveText(page, el("gallery counter"), "1"));
      await session.step(125, "And \"BddScriptCreate{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BddScriptCreate{time}\" link in gallery")), "visible"));
      await session.step(126, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The JavaScript template raises its alert", async () => {
      await session.step(129, "Given browser alerts are recorded", () => recordAlerts(page));
      await session.step(130, "When user clicks on New button", () => clickOn(page, el("New button")));
      await session.step(131, "And user picks \"JavaScript Script...\" from the open menu", () => pickFromOpenMenu(page, "JavaScript Script..."));
      await session.step(132, "Then the \"Template\" view should be current", () => viewIsCurrent(page, "Template"));
      await session.step(133, "And code editor should contain the text \"alert('Hello World!')\"", () => shouldContainText(page, el("code editor"), "alert('Hello World!')"));
      await session.step(134, "When user clicks on \"Run script (F5)\" icon", () => clickOn(page, el("\"Run script (F5)\" icon")));
      await session.step(135, "Then the browser should have shown the alert \"Hello World!\"", () => alertShown(page, "Hello World!"));
      await session.step(136, "And no errors should have been logged", () => noErrors(page));
      await session.step(137, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
