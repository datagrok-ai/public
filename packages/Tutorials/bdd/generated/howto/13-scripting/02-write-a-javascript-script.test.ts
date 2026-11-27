/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/13-scripting/02-write-a-javascript-script.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-data.js';
import '../../../bindings/howto-databases.js';
import '../../../bindings/howto-filters.js';
import '../../../bindings/howto-ml.js';
import '../../../bindings/howto-molecules.js';
import '../../../bindings/howto-more-viewers.js';
import '../../../bindings/howto-sequences.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {scriptResultReads, scriptsSearchPutBack} from '../../../bindings/howto-scripting.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, expand, insertLineAfter, replaceLine, shouldBe, shouldContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, contextPanelClosed, dialogCloses, noScriptOnServer, saveScript, scriptsOnServer, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Scripting", () => {
  const session = feature(test, "features/howto/13-scripting/02-write-a-javascript-script.feature", import.meta.url);
  test("Write a JavaScript script", {tag: ["@guide", "@tutorial", "@viewport:1280x720"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And no script named \"Body mass index\" is on the server", () => noScriptOnServer(page, "Body mass index"));
    await session.step(19, "And the Scripts gallery search is emptied at feature end", () => scriptsSearchPutBack(page));
    await session.step(20, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(21, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(22, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(24, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")), undefined, "Scripts live in Browse, under Platform > Functions");
    await session.step(25, "And user expands \"Platform > Functions\" tree node inside browse tree", () => expand(page, el("\"Platform > Functions\" tree node inside browse tree")));
    await session.step(26, "And user clicks on \"Platform > Functions > Scripts\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Functions > Scripts\" tree node inside browse tree")));
    await session.step(27, "Then the \"Scripts\" view should be current", () => viewIsCurrent(page, "Scripts"));
    await session.step(29, "When user clicks on New button", () => clickOn(page, el("New button")), undefined, "New > JavaScript Script: it runs in your browser");
    await session.step(30, "And user picks \"JavaScript Script...\" from the open menu", () => pickFromOpenMenu(page, "JavaScript Script..."));
    await session.step(31, "Then the \"Template\" view should be current", () => viewIsCurrent(page, "Template"));
    await session.step(32, "And code editor should contain the text \"alert('Hello World!')\"", () => shouldContainText(page, el("code editor"), "alert('Hello World!')"));
    await session.step(34, "When user replaces the line starting with \"//name:\" in code editor with \"//name: Body mass index\"", () => replaceLine(page, "//name:", el("code editor"), "//name: Body mass index"), undefined, "The comment lines at the top describe the script");
    await session.step(35, "And user replaces the line starting with \"//description:\" in code editor with \"//description: BMI from height and weight\"", () => replaceLine(page, "//description:", el("code editor"), "//description: BMI from height and weight"));
    await session.step(37, "And user puts \"//input: double height = 175 {caption: Height, cm}\" on a new line after the line starting with \"//language:\" in code editor", () => insertLineAfter(page, "//input: double height = 175 {caption: Height, cm}", "//language:", el("code editor")), undefined, "Inputs: type, name, default value and a caption");
    await session.step(38, "And user puts \"//input: double weight = 70 {caption: Weight, kg}\" on a new line after the line starting with \"//input: double height\" in code editor", () => insertLineAfter(page, "//input: double weight = 70 {caption: Weight, kg}", "//input: double height", el("code editor")));
    await session.step(40, "And user puts \"//output: double bmi\" on a new line after the line starting with \"//input: double weight\" in code editor", () => insertLineAfter(page, "//output: double bmi", "//input: double weight", el("code editor")), undefined, "The output is a variable the code sets");
    await session.step(41, "And user replaces the line starting with \"alert(\" in code editor with \"bmi = weight / (height / 100) ** 2;\"", () => replaceLine(page, "alert(", el("code editor"), "bmi = weight / (height / 100) ** 2;"));
    await session.step(43, "And user clicks on \"Run script (F5)\" icon", () => clickOn(page, el("\"Run script (F5)\" icon")), undefined, "Run: the inputs come as a form");
    await session.step(44, "Then \"Body mass index\" dialog should be visible", () => shouldBe(page, el("\"Body mass index\" dialog"), "visible"));
    await session.step(45, "When user enters \"180\" into \"Height, cm\" input in \"Body mass index\" dialog", () => enterInto(page, "180", el("\"Height, cm\" input in \"Body mass index\" dialog")));
    await session.step(46, "And user clicks on OK button in \"Body mass index\" dialog", () => clickOn(page, el("OK button in \"Body mass index\" dialog")));
    await session.step(47, "Then the \"Body mass index\" dialog should close", () => dialogCloses(page, "Body mass index"));
    await session.step(49, "And the script result \"bmi\" should read \"21.60\"", () => scriptResultReads(page, "bmi", "21.60"), undefined, "The result shows under the code");
    await session.step(51, "When user saves the script", () => saveScript(page), undefined, "Save it to run it again later");
    await session.step(52, "Then 1 script named \"Body mass index\" should be on the server", () => scriptsOnServer(page, 1, "Body mass index"));
    await session.step(54, "When user clicks on \"Scripts\" view tab", () => clickOn(page, el("\"Scripts\" view tab")), undefined, "It is listed in Scripts with the others");
    await session.step(55, "And user types \"Body mass\" into gallery search", () => typeInto(page, "Body mass", el("gallery search")));
    await session.step(56, "Then \"Body mass index\" link in gallery should be visible", () => shouldBe(page, el("\"Body mass index\" link in gallery"), "visible"));
    await session.step(57, "And no errors should have been logged", () => noErrors(page));
  });
});
