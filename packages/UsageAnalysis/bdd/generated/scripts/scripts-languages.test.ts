/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/scripts/scripts-languages.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.scripts]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {scriptResult} from '../../bindings/scripts.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {alertShown, clickOn, recordAlerts, selectIn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {tableOpen} from '@datagrok-libraries/bdd/bindings/platform/data';
import {dialogCloses, scriptsView, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A script in every language runs its template", () => {
  const session = feature(test, "features/scripts/scripts-languages.feature", import.meta.url);
  test("The R template runs with cars [language=R]", {tag: ["@serial", "@realizes:views.scripts", "@full-stand"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And user opens the Scripts view", () => scriptsView(page));
    await session.step(25, "When user clicks on New button", () => clickOn(page, el("New button")));
    await session.step(26, "And user picks \"R Script...\" from the open menu", () => pickFromOpenMenu(page, "R Script..."));
    await session.step(27, "Then the \"Template\" view should be current", () => viewIsCurrent(page, "Template"));
    await session.step(28, "And code editor should contain the text \"sample: cars.csv\"", () => shouldContainText(page, el("code editor"), "sample: cars.csv"));
    await session.step(29, "When user clicks on \"Open script sample table\" icon", () => clickOn(page, el("\"Open script sample table\" icon")));
    await session.step(30, "Then table \"cars\" should be open", () => tableOpen(page, "cars"));
    await session.step(31, "When user clicks on \"Run script (F5)\" icon", () => clickOn(page, el("\"Run script (F5)\" icon")));
    await session.step(32, "Then \"Template\" dialog should be visible", () => shouldBe(page, el("\"Template\" dialog"), "visible"));
    await session.step(33, "When user selects \"cars\" in Table input in \"Template\" dialog", () => selectIn(page, "cars", el("Table input in \"Template\" dialog")));
    await session.step(34, "And user clicks on OK button in \"Template\" dialog", () => clickOn(page, el("OK button in \"Template\" dialog")));
    await session.step(35, "Then the \"Template\" dialog should close", () => dialogCloses(page, "Template"));
    await session.step(36, "And the script results should show \"count\" as \"510\"", () => scriptResult(page, "count", "510"));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
    await session.step(38, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Python template runs with cars [language=Python]", {tag: ["@serial", "@realizes:views.scripts", "@full-stand"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And user opens the Scripts view", () => scriptsView(page));
    await session.step(25, "When user clicks on New button", () => clickOn(page, el("New button")));
    await session.step(26, "And user picks \"Python Script...\" from the open menu", () => pickFromOpenMenu(page, "Python Script..."));
    await session.step(27, "Then the \"Template\" view should be current", () => viewIsCurrent(page, "Template"));
    await session.step(28, "And code editor should contain the text \"sample: cars.csv\"", () => shouldContainText(page, el("code editor"), "sample: cars.csv"));
    await session.step(29, "When user clicks on \"Open script sample table\" icon", () => clickOn(page, el("\"Open script sample table\" icon")));
    await session.step(30, "Then table \"cars\" should be open", () => tableOpen(page, "cars"));
    await session.step(31, "When user clicks on \"Run script (F5)\" icon", () => clickOn(page, el("\"Run script (F5)\" icon")));
    await session.step(32, "Then \"Template\" dialog should be visible", () => shouldBe(page, el("\"Template\" dialog"), "visible"));
    await session.step(33, "When user selects \"cars\" in Table input in \"Template\" dialog", () => selectIn(page, "cars", el("Table input in \"Template\" dialog")));
    await session.step(34, "And user clicks on OK button in \"Template\" dialog", () => clickOn(page, el("OK button in \"Template\" dialog")));
    await session.step(35, "Then the \"Template\" dialog should close", () => dialogCloses(page, "Template"));
    await session.step(36, "And the script results should show \"count\" as \"510\"", () => scriptResult(page, "count", "510"));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
    await session.step(38, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The NodeJS template runs with cars [language=NodeJS]", {tag: ["@serial", "@realizes:views.scripts", "@full-stand"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And user opens the Scripts view", () => scriptsView(page));
    await session.step(25, "When user clicks on New button", () => clickOn(page, el("New button")));
    await session.step(26, "And user picks \"NodeJS Script...\" from the open menu", () => pickFromOpenMenu(page, "NodeJS Script..."));
    await session.step(27, "Then the \"Template\" view should be current", () => viewIsCurrent(page, "Template"));
    await session.step(28, "And code editor should contain the text \"sample: cars.csv\"", () => shouldContainText(page, el("code editor"), "sample: cars.csv"));
    await session.step(29, "When user clicks on \"Open script sample table\" icon", () => clickOn(page, el("\"Open script sample table\" icon")));
    await session.step(30, "Then table \"cars\" should be open", () => tableOpen(page, "cars"));
    await session.step(31, "When user clicks on \"Run script (F5)\" icon", () => clickOn(page, el("\"Run script (F5)\" icon")));
    await session.step(32, "Then \"Template\" dialog should be visible", () => shouldBe(page, el("\"Template\" dialog"), "visible"));
    await session.step(33, "When user selects \"cars\" in Table input in \"Template\" dialog", () => selectIn(page, "cars", el("Table input in \"Template\" dialog")));
    await session.step(34, "And user clicks on OK button in \"Template\" dialog", () => clickOn(page, el("OK button in \"Template\" dialog")));
    await session.step(35, "Then the \"Template\" dialog should close", () => dialogCloses(page, "Template"));
    await session.step(36, "And the script results should show \"count\" as \"510\"", () => scriptResult(page, "count", "510"));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
    await session.step(38, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Julia template runs with cars [language=Julia]", {tag: ["@serial", "@realizes:views.scripts", "@full-stand"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And user opens the Scripts view", () => scriptsView(page));
    await session.step(25, "When user clicks on New button", () => clickOn(page, el("New button")));
    await session.step(26, "And user picks \"Julia Script...\" from the open menu", () => pickFromOpenMenu(page, "Julia Script..."));
    await session.step(27, "Then the \"Template\" view should be current", () => viewIsCurrent(page, "Template"));
    await session.step(28, "And code editor should contain the text \"sample: cars.csv\"", () => shouldContainText(page, el("code editor"), "sample: cars.csv"));
    await session.step(29, "When user clicks on \"Open script sample table\" icon", () => clickOn(page, el("\"Open script sample table\" icon")));
    await session.step(30, "Then table \"cars\" should be open", () => tableOpen(page, "cars"));
    await session.step(31, "When user clicks on \"Run script (F5)\" icon", () => clickOn(page, el("\"Run script (F5)\" icon")));
    await session.step(32, "Then \"Template\" dialog should be visible", () => shouldBe(page, el("\"Template\" dialog"), "visible"));
    await session.step(33, "When user selects \"cars\" in Table input in \"Template\" dialog", () => selectIn(page, "cars", el("Table input in \"Template\" dialog")));
    await session.step(34, "And user clicks on OK button in \"Template\" dialog", () => clickOn(page, el("OK button in \"Template\" dialog")));
    await session.step(35, "Then the \"Template\" dialog should close", () => dialogCloses(page, "Template"));
    await session.step(36, "And the script results should show \"count\" as \"510\"", () => scriptResult(page, "count", "510"));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
    await session.step(38, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Grok template runs with cars [language=Grok]", {tag: ["@serial", "@realizes:views.scripts"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And user opens the Scripts view", () => scriptsView(page));
    await session.step(25, "When user clicks on New button", () => clickOn(page, el("New button")));
    await session.step(26, "And user picks \"Grok Script...\" from the open menu", () => pickFromOpenMenu(page, "Grok Script..."));
    await session.step(27, "Then the \"Template\" view should be current", () => viewIsCurrent(page, "Template"));
    await session.step(28, "And code editor should contain the text \"sample: cars.csv\"", () => shouldContainText(page, el("code editor"), "sample: cars.csv"));
    await session.step(29, "When user clicks on \"Open script sample table\" icon", () => clickOn(page, el("\"Open script sample table\" icon")));
    await session.step(30, "Then table \"cars\" should be open", () => tableOpen(page, "cars"));
    await session.step(31, "When user clicks on \"Run script (F5)\" icon", () => clickOn(page, el("\"Run script (F5)\" icon")));
    await session.step(32, "Then \"Template\" dialog should be visible", () => shouldBe(page, el("\"Template\" dialog"), "visible"));
    await session.step(33, "When user selects \"cars\" in Table input in \"Template\" dialog", () => selectIn(page, "cars", el("Table input in \"Template\" dialog")));
    await session.step(34, "And user clicks on OK button in \"Template\" dialog", () => clickOn(page, el("OK button in \"Template\" dialog")));
    await session.step(35, "Then the \"Template\" dialog should close", () => dialogCloses(page, "Template"));
    await session.step(36, "And the script results should show \"count\" as \"510\"", () => scriptResult(page, "count", "510"));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
    await session.step(38, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Pyodide template runs with cars [language=Pyodide]", {tag: ["@serial", "@realizes:views.scripts"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And user opens the Scripts view", () => scriptsView(page));
    await session.step(25, "When user clicks on New button", () => clickOn(page, el("New button")));
    await session.step(26, "And user picks \"Pyodide Script...\" from the open menu", () => pickFromOpenMenu(page, "Pyodide Script..."));
    await session.step(27, "Then the \"Template\" view should be current", () => viewIsCurrent(page, "Template"));
    await session.step(28, "And code editor should contain the text \"sample: cars.csv\"", () => shouldContainText(page, el("code editor"), "sample: cars.csv"));
    await session.step(29, "When user clicks on \"Open script sample table\" icon", () => clickOn(page, el("\"Open script sample table\" icon")));
    await session.step(30, "Then table \"cars\" should be open", () => tableOpen(page, "cars"));
    await session.step(31, "When user clicks on \"Run script (F5)\" icon", () => clickOn(page, el("\"Run script (F5)\" icon")));
    await session.step(32, "Then \"Template\" dialog should be visible", () => shouldBe(page, el("\"Template\" dialog"), "visible"));
    await session.step(33, "When user selects \"cars\" in Table input in \"Template\" dialog", () => selectIn(page, "cars", el("Table input in \"Template\" dialog")));
    await session.step(34, "And user clicks on OK button in \"Template\" dialog", () => clickOn(page, el("OK button in \"Template\" dialog")));
    await session.step(35, "Then the \"Template\" dialog should close", () => dialogCloses(page, "Template"));
    await session.step(36, "And the script results should show \"count\" as \"510\"", () => scriptResult(page, "count", "510"));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
    await session.step(38, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The JavaScript template raises its alert", {tag: ["@serial", "@realizes:views.scripts"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And user opens the Scripts view", () => scriptsView(page));
    await session.step(54, "Given browser alerts are recorded", () => recordAlerts(page));
    await session.step(55, "When user clicks on New button", () => clickOn(page, el("New button")));
    await session.step(56, "And user picks \"JavaScript Script...\" from the open menu", () => pickFromOpenMenu(page, "JavaScript Script..."));
    await session.step(57, "Then the \"Template\" view should be current", () => viewIsCurrent(page, "Template"));
    await session.step(58, "And code editor should contain the text \"alert('Hello World!')\"", () => shouldContainText(page, el("code editor"), "alert('Hello World!')"));
    await session.step(59, "When user clicks on \"Run script (F5)\" icon", () => clickOn(page, el("\"Run script (F5)\" icon")));
    await session.step(60, "Then the browser should have shown the alert \"Hello World!\"", () => alertShown(page, "Hello World!"));
    await session.step(61, "And no errors should have been logged", () => noErrors(page));
    await session.step(62, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
