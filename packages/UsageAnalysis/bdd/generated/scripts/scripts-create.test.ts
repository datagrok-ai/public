/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/scripts/scripts-create.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.scripts]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {saveScript, scriptResult} from '../../bindings/scripts.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, expand, selectIn, shouldBe, shouldContainText, typeInto, visibleCount} from '@datagrok-libraries/bdd/bindings/common/steps';
import {tableOpen, tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeCurrentView, dialogCloses, galleryCountLower, noScriptOnServer, rememberGalleryCount, scriptHasParam, scriptsOnServer, switchView, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, menuLists, noBalloons, noErrors, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Creating a script", () => {
  const session = feature(test, "features/scripts/scripts-create.feature", import.meta.url);
  test("Creating a script", {tag: ["@journey", "@full-stand", "@serial", "@realizes:views.scripts"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And no script named \"BddScriptCreate{time}\" is on the server", () => noScriptOnServer(page, session.text("BddScriptCreate{time}")));
    await session.step(26, "And the browse panel is open", () => browsePanelOpen(page));
    await run.scenario("The Scripts view opens from Platform > Functions", async () => {
      await session.step(29, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
      await session.step(30, "And user expands \"Platform > Functions\" tree node inside browse tree", () => expand(page, el("\"Platform > Functions\" tree node inside browse tree")));
      await session.step(31, "And user clicks on \"Platform > Functions > Scripts\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Functions > Scripts\" tree node inside browse tree")));
      await session.step(32, "Then the \"Scripts\" view should be current", () => viewIsCurrent(page, "Scripts"));
      await session.step(33, "And gallery should be visible", () => shouldBe(page, el("gallery"), "visible"));
      await session.step(34, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("New offers every language and opens the R template unsaved", async () => {
      await session.step(37, "When user clicks on New button", () => clickOn(page, el("New button")));
      await session.step(38, "Then the open menu should list \"R Script...\"", () => menuLists(page, "R Script..."));
      await session.step(39, "And the open menu should list \"Python Script...\"", () => menuLists(page, "Python Script..."));
      await session.step(40, "And the open menu should list \"Octave Script...\"", () => menuLists(page, "Octave Script..."));
      await session.step(41, "And the open menu should list \"NodeJS Script...\"", () => menuLists(page, "NodeJS Script..."));
      await session.step(42, "And the open menu should list \"Julia Script...\"", () => menuLists(page, "Julia Script..."));
      await session.step(43, "And the open menu should list \"JavaScript Script...\"", () => menuLists(page, "JavaScript Script..."));
      await session.step(44, "And the open menu should list \"Grok Script...\"", () => menuLists(page, "Grok Script..."));
      await session.step(45, "And the open menu should list \"Pyodide Script...\"", () => menuLists(page, "Pyodide Script..."));
      await session.step(46, "When user picks \"R Script...\" from the open menu", () => pickFromOpenMenu(page, "R Script..."));
      await session.step(47, "Then the \"Template\" view should be current", () => viewIsCurrent(page, "Template"));
      await session.step(48, "And code editor should contain the text \"#language: r\"", () => shouldContainText(page, el("code editor"), "#language: r"));
      await session.step(49, "And code editor should be visible", () => shouldBe(page, el("code editor"), "visible"));
      await session.step(50, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The sample icon opens the sample table", async () => {
      await session.step(53, "When user clicks on \"Open script sample table\" icon", () => clickOn(page, el("\"Open script sample table\" icon")));
      await session.step(54, "Then table \"cars\" should be open", () => tableOpen(page, "cars"));
      await session.step(55, "And table \"cars\" should have 30 rows", () => tableRows(page, "cars", 30));
      await session.step(56, "And the \"Template\" view should be current", () => viewIsCurrent(page, "Template"));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Signature editor names the script and adds a parameter to its header", async () => {
      await session.step(60, "Given user switches to the \"Template\" view", () => switchView(page, "Template"));
      await session.step(61, "When user clicks on \"Open Signature Editor\" icon", () => clickOn(page, el("\"Open Signature Editor\" icon")));
      await session.step(62, "Then PARAMETERS tab should be visible", () => shouldBe(page, el("PARAMETERS tab"), "visible"));
      await session.step(63, "When user enters \"BddScriptCreate{time}\" into Name input", () => enterInto(page, session.text("BddScriptCreate{time}"), el("Name input")));
      await session.step(64, "And user clicks on PARAMETERS tab", () => clickOn(page, el("PARAMETERS tab")));
      await session.step(65, "Then there should be 2 visible \"Add the param\" icon", () => visibleCount(page, 2, el("\"Add the param\" icon")));
      await session.step(66, "When user clicks on first \"Add the param\" icon", () => clickOn(page, el("first \"Add the param\" icon")));
      await session.step(67, "Then there should be 3 visible \"Add the param\" icon", () => visibleCount(page, 3, el("\"Add the param\" icon")));
      await session.step(68, "Then \"Open function editor\" icon should be visible", () => shouldBe(page, el("\"Open function editor\" icon"), "visible"));
      await session.step(69, "When user clicks on \"Open function editor\" icon", () => clickOn(page, el("\"Open function editor\" icon")));
      await session.step(70, "Then \"Run script (F5)\" icon should be visible", () => shouldBe(page, el("\"Run script (F5)\" icon"), "visible"));
      await session.step(71, "And code editor should contain the text \"#name: BddScriptCreate{time}\"", () => shouldContainText(page, el("code editor"), session.text("#name: BddScriptCreate{time}")));
      await session.step(72, "And code editor should contain the text \"#input: bool newParam\"", () => shouldContainText(page, el("code editor"), "#input: bool newParam"));
      await session.step(73, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Running with cars answers the number of cells", async () => {
      await session.step(77, "Given user switches to the \"Template\" view", () => switchView(page, "Template"));
      await session.step(78, "When user clicks on \"Run script (F5)\" icon", () => clickOn(page, el("\"Run script (F5)\" icon")));
      await session.step(79, "Then \"Template\" dialog should be visible", () => shouldBe(page, el("\"Template\" dialog"), "visible"));
      await session.step(80, "When user selects \"cars\" in Table input in \"Template\" dialog", () => selectIn(page, "cars", el("Table input in \"Template\" dialog")));
      await session.step(81, "And user clicks on OK button in \"Template\" dialog", () => clickOn(page, el("OK button in \"Template\" dialog")));
      await session.step(82, "Then the \"Template\" dialog should close", () => dialogCloses(page, "Template"));
      await session.step(83, "And the script results should show \"count\" as \"510\"", () => scriptResult(page, "count", "510"));
      await session.step(84, "And no errors should have been logged", () => noErrors(page));
      await session.step(85, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Save stores the script with the parameters of its header", async () => {
      await session.step(88, "When user saves the script", () => saveScript(page));
      await session.step(89, "Then an info balloon containing \"Script saved.\" should have been shown", () => infoBalloonText(page, "Script saved."));
      await session.step(90, "And 1 script named \"BddScriptCreate{time}\" should be on the server", () => scriptsOnServer(page, 1, session.text("BddScriptCreate{time}")));
      await session.step(91, "And the script \"BddScriptCreate{time}\" on the server should have an output \"count\" of type \"int\"", () => scriptHasParam(page, session.text("BddScriptCreate{time}"), "output", "count", "int"));
      await session.step(92, "And the script \"BddScriptCreate{time}\" on the server should have an input \"newParam\" of type \"bool\"", () => scriptHasParam(page, session.text("BddScriptCreate{time}"), "input", "newParam", "bool"));
      await session.step(93, "And the script \"BddScriptCreate{time}\" on the server should have an input \"table\" of type \"dataframe\"", () => scriptHasParam(page, session.text("BddScriptCreate{time}"), "input", "table", "dataframe"));
      await session.step(94, "And the \"BddScriptCreate{time}\" view should be current", () => viewIsCurrent(page, session.text("BddScriptCreate{time}")));
      await session.step(95, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Closing the editor returns to Scripts, where the new script is found", async () => {
      await session.step(98, "When user closes the current view", () => closeCurrentView(page));
      await session.step(99, "Then the \"Scripts\" view should be current", () => viewIsCurrent(page, "Scripts"));
      await session.step(100, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(101, "And user types \"BddScriptCreate{time}\" into gallery search", () => typeInto(page, session.text("BddScriptCreate{time}"), el("gallery search")));
      await session.step(102, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(103, "And \"BddScriptCreate{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BddScriptCreate{time}\" link in gallery")), "visible"));
      await session.step(104, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
