/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/catalog.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [diffstudio.app.diff-studio]
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openLibraryModel, openModelHub, saveToLibrary} from '../bindings/diff-studio.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {customFired, listenCustom} from '@datagrok-libraries/bdd/bindings/platform/events';
import {viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, readingReads, repainted, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A model saved to the library", () => {
  const session = feature(test, "features/catalog.feature", import.meta.url);
  test("A model saved to the library", {tag: ["@journey", "@diffstudio", "@realizes:diffstudio.app.diff-studio"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens the \"PK-PD\" model of the Diff Studio library", () => openLibraryModel(page, "PK-PD"));
    await run.scenario("The model is the one the library serves", async () => {
      await session.step(24, "Then the \"PK-PD\" view should be current", () => viewIsCurrent(page, "PK-PD"));
      await session.step(25, "And dose input should be visible", () => shouldBe(page, el("dose input"), "visible"));
    });
    await run.scenario("The model answers to its inputs", async () => {
      await session.step(28, "When user clicks on Multiaxis tab", () => clickOn(page, el("Multiaxis tab")));
      await session.step(29, "And user takes a snapshot of line chart viewer", () => takeSnapshot(page, el("line chart viewer")));
      await session.step(30, "And user enters \"5000\" into dose input", () => enterInto(page, "5000", el("dose input")));
      await session.step(31, "Then dose input should have value \"5000\"", () => shouldHaveValue(page, el("dose input"), "5000"));
      await session.step(32, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
    });
    await run.scenario("Saving it to the library is announced", async () => {
      await session.step(35, "Given user listens for \"diff-studio:library-changed\" custom event", () => listenCustom(page, "diff-studio:library-changed"));
      await session.step(36, "When user saves the model to the Diff Studio library", () => saveToLibrary(page));
      await session.step(37, "Then the \"diff-studio:library-changed\" custom event should have fired", () => customFired(page, "diff-studio:library-changed"));
    });
    await run.scenario("The Model Hub lists the saved model", async () => {
      await session.step(40, "Given user opens the Model Hub", () => openModelHub(page));
      await session.step(41, "Then PK-PD link in gallery should be visible", () => shouldBe(page, el("PK-PD link in gallery"), "visible"));
    });
    await run.scenario("The model runs from the catalog and its chart follows its inputs", async () => {
      await session.step(44, "When user double-clicks on PK-PD link in gallery", () => doubleClickOn(page, el("PK-PD link in gallery")));
      await session.step(45, "Then dose input should be visible", () => shouldBe(page, el("dose input"), "visible"));
      await session.step(46, "And \"PK-PD / DiffStudio Facet\" tab should be selected", () => shouldBe(page, el("\"PK-PD / DiffStudio Facet\" tab"), "selected"));
      await session.step(47, "When user takes a snapshot of line chart viewer", () => takeSnapshot(page, el("line chart viewer")));
      await session.step(48, "And user enters \"5000\" into dose input", () => enterInto(page, "5000", el("dose input")));
      await session.step(49, "Then dose input should have value \"5000\"", () => shouldHaveValue(page, el("dose input"), "5000"));
      await session.step(50, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(51, "When user clicks on \"PK-PD / Grid\" tab", () => clickOn(page, el("\"PK-PD / Grid\" tab")));
      await session.step(52, "Then \"PK-PD / Grid\" tab should be selected", () => shouldBe(page, el("\"PK-PD / Grid\" tab"), "selected"));
      await session.step(53, "And the \"rows\" reading of grid viewer should be \"610\"", () => readingReads(page, "rows", el("grid viewer"), "610"));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
