/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/scripting.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [diffstudio.app.diff-studio]
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {hubListsScript, openLibraryModel, openModelHub, openSavedScript, saveScript} from '../bindings/diff-studio.js';
import {lookedDifferent, takePicture} from '@datagrok-libraries/bdd/bindings/common/pixels';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, insertLine, shouldBe, shouldContainText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The equations behind a model, and the script they become", () => {
  const session = feature(test, "features/scripting.feature", import.meta.url);
  test("The equations behind a model, and the script they become", {tag: ["@journey", "@diffstudio", "@realizes:diffstudio.app.diff-studio"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens the \"Bioreactor\" model of the Diff Studio library", () => openLibraryModel(page, "Bioreactor"));
    await run.scenario("The model shows its inputs before anything is edited", async () => {
      await session.step(21, "Then the \"Bioreactor\" view should be current", () => viewIsCurrent(page, "Bioreactor"));
      await session.step(22, "And \"Process mode\" input should be visible", () => shouldBe(page, el("\"Process mode\" input"), "visible"));
      await session.step(23, "And code editor should be absent", () => shouldBe(page, el("code editor"), "absent"));
    });
    await run.scenario("Edit opens the equations editor in place of the form", async () => {
      await session.step(26, "When user clicks on Edit ribbon item", () => clickOn(page, el("Edit ribbon item")));
      await session.step(27, "Then code editor should be visible", () => shouldBe(page, el("code editor"), "visible"));
      await session.step(28, "And \"Process mode\" input should be absent", () => shouldBe(page, el("\"Process mode\" input"), "absent"));
    });
    await run.scenario("The angle brackets turn the model into a script", async () => {
      await session.step(31, "When user clicks on \"</>\" ribbon item", () => clickOn(page, el("\"</>\" ribbon item")));
      await session.step(32, "Then Sensitivity ribbon item should be absent", () => shouldBe(page, el("Sensitivity ribbon item"), "absent"));
      await session.step(33, "And \"Run script (F5)\" icon should be visible", () => shouldBe(page, el("\"Run script (F5)\" icon"), "visible"));
    });
    await run.scenario("The script is tagged as a model and saved", async () => {
      await session.step(36, "When user puts \"//tags: model\" on the first line of code editor", () => insertLine(page, "//tags: model", el("code editor")));
      await session.step(37, "Then code editor should contain the text \"//tags: model\"", () => shouldContainText(page, el("code editor"), "//tags: model"));
      await session.step(38, "When user saves the script", () => saveScript(page));
    });
    await run.scenario("The script runs, and its chart follows the input it exposes", async () => {
      await session.step(41, "When user clicks on \"Run script (F5)\" icon", () => clickOn(page, el("\"Run script (F5)\" icon")));
      await session.step(42, "Then Final input should be visible", () => shouldBe(page, el("Final input"), "visible"));
      await session.step(43, "When user takes a picture of viewer", () => takePicture(page, el("viewer")));
      await session.step(44, "And user enters \"500\" into Final input", () => enterInto(page, "500", el("Final input")));
      await session.step(45, "Then Final input should have value \"500\"", () => shouldHaveValue(page, el("Final input"), "500"));
      await session.step(46, "And viewer should look different", () => lookedDifferent(page, el("viewer")));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Model Hub lists the saved script, and it answers to its inputs there", async () => {
      await session.step(50, "Given user opens the Model Hub", () => openModelHub(page));
      await session.step(51, "Then the Model Hub should list the saved script", () => hubListsScript(page));
      await session.step(52, "When user opens the saved script from the Model Hub", () => openSavedScript(page));
      await session.step(53, "Then Final input should be visible", () => shouldBe(page, el("Final input"), "visible"));
      await session.step(54, "When user takes a picture of viewer", () => takePicture(page, el("viewer")));
      await session.step(55, "And user enters \"800\" into Final input", () => enterInto(page, "800", el("Final input")));
      await session.step(56, "Then Final input should have value \"800\"", () => shouldHaveValue(page, el("Final input"), "800"));
      await session.step(57, "And viewer should look different", () => lookedDifferent(page, el("viewer")));
      await session.step(58, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
