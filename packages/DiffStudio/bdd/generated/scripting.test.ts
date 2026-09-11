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
import {deleteSavedScript, hubDoesNotListScript, hubListsScript, openLibraryModel, openModelHub, openSavedScript, saveScript} from '../bindings/diff-studio.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, insertLine, shouldBe, shouldContainText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, readingReads, repainted, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The equations behind a model, and the script they become", () => {
  const session = feature(test, "features/scripting.feature", import.meta.url);
  test("The equations behind a model, and the script they become", {tag: ["@journey", "@diffstudio", "@realizes:diffstudio.app.diff-studio"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens the \"Bioreactor\" model of the Diff Studio library", () => openLibraryModel(page, "Bioreactor"));
    await run.scenario("The model shows its inputs before anything is edited", async () => {
      await session.step(22, "Then the \"Bioreactor\" view should be current", () => viewIsCurrent(page, "Bioreactor"));
      await session.step(23, "And \"Process mode\" input should be visible", () => shouldBe(page, el("\"Process mode\" input"), "visible"));
      await session.step(24, "And code editor should be absent", () => shouldBe(page, el("code editor"), "absent"));
    });
    await run.scenario("Edit opens the equations editor in place of the form", async () => {
      await session.step(27, "When user clicks on Edit ribbon item", () => clickOn(page, el("Edit ribbon item")));
      await session.step(28, "Then code editor should be visible", () => shouldBe(page, el("code editor"), "visible"));
      await session.step(29, "And \"Process mode\" input should be absent", () => shouldBe(page, el("\"Process mode\" input"), "absent"));
    });
    await run.scenario("The angle brackets turn the model into a script", async () => {
      await session.step(32, "When user clicks on \"</>\" ribbon item", () => clickOn(page, el("\"</>\" ribbon item")));
      await session.step(33, "Then Sensitivity ribbon item should be absent", () => shouldBe(page, el("Sensitivity ribbon item"), "absent"));
      await session.step(34, "And \"Run script (F5)\" icon should be visible", () => shouldBe(page, el("\"Run script (F5)\" icon"), "visible"));
    });
    await run.scenario("The script is tagged as a model and saved", async () => {
      await session.step(37, "When user puts \"//tags: model\" on the first line of code editor", () => insertLine(page, "//tags: model", el("code editor")));
      await session.step(38, "Then code editor should contain the text \"//tags: model\"", () => shouldContainText(page, el("code editor"), "//tags: model"));
      await session.step(39, "When user saves the script", () => saveScript(page));
    });
    await run.scenario("The script runs, and its table and chart follow the input it exposes", async () => {
      await session.step(42, "When user clicks on \"Run script (F5)\" icon", () => clickOn(page, el("\"Run script (F5)\" icon")));
      await session.step(43, "Then Final input should be visible", () => shouldBe(page, el("Final input"), "visible"));
      await session.step(44, "And \"Bioreactor / Grid\" tab should be selected", () => shouldBe(page, el("\"Bioreactor / Grid\" tab"), "selected"));
      await session.step(45, "And the \"rows\" reading of grid viewer should be \"1001\"", () => readingReads(page, "rows", el("grid viewer"), "1001"));
      await session.step(46, "When user clicks on \"Bioreactor / DiffStudio Facet\" tab", () => clickOn(page, el("\"Bioreactor / DiffStudio Facet\" tab")));
      await session.step(47, "And user takes a snapshot of line chart viewer", () => takeSnapshot(page, el("line chart viewer")));
      await session.step(48, "And user enters \"500\" into Final input", () => enterInto(page, "500", el("Final input")));
      await session.step(49, "Then Final input should have value \"500\"", () => shouldHaveValue(page, el("Final input"), "500"));
      await session.step(50, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(51, "When user clicks on \"Bioreactor / Grid\" tab", () => clickOn(page, el("\"Bioreactor / Grid\" tab")));
      await session.step(52, "Then the \"rows\" reading of grid viewer should be \"501\"", () => readingReads(page, "rows", el("grid viewer"), "501"));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Model Hub lists the saved script, and it answers to its inputs there", async () => {
      await session.step(56, "Given user opens the Model Hub", () => openModelHub(page));
      await session.step(57, "Then the Model Hub should list the saved script", () => hubListsScript(page));
      await session.step(58, "When user opens the saved script from the Model Hub", () => openSavedScript(page));
      await session.step(59, "Then Final input should be visible", () => shouldBe(page, el("Final input"), "visible"));
      await session.step(60, "And \"Bioreactor / Grid\" tab should be selected", () => shouldBe(page, el("\"Bioreactor / Grid\" tab"), "selected"));
      await session.step(61, "When user clicks on \"Bioreactor / DiffStudio Facet\" tab", () => clickOn(page, el("\"Bioreactor / DiffStudio Facet\" tab")));
      await session.step(62, "And user takes a snapshot of line chart viewer", () => takeSnapshot(page, el("line chart viewer")));
      await session.step(63, "And user enters \"800\" into Final input", () => enterInto(page, "800", el("Final input")));
      await session.step(64, "Then Final input should have value \"800\"", () => shouldHaveValue(page, el("Final input"), "800"));
      await session.step(65, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(66, "When user clicks on \"Bioreactor / Grid\" tab", () => clickOn(page, el("\"Bioreactor / Grid\" tab")));
      await session.step(67, "Then the \"rows\" reading of grid viewer should be \"801\"", () => readingReads(page, "rows", el("grid viewer"), "801"));
      await session.step(68, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Refresh re-fetches the catalog, so a model removed behind its back disappears", async () => {
      await session.step(71, "When the saved script is deleted on the server", () => deleteSavedScript(page));
      await session.step(72, "And user clicks on model hub refresh icon", () => clickOn(page, el("model hub refresh icon")));
      await session.step(73, "Then the Model Hub should not list the saved script", () => hubDoesNotListScript(page));
      await session.step(74, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
