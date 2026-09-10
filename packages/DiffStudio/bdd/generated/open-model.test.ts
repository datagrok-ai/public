/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/open-model.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [diffstudio.app.diff-studio, diffstudio.model.bioreactor]
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openLibraryModel} from '../bindings/diff-studio.js';
import {canvasColors} from '@datagrok-libraries/bdd/bindings/common/pixels';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, dragSliderTo, enterInto, selectIn, shouldBe, shouldHaveValue, shouldHaveValueBetween, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, repainted, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Opening a model from the library", () => {
  const session = feature(test, "features/open-model.feature", import.meta.url);
  test("Opening a model from the library", {tag: ["@journey", "@diffstudio", "@realizes:diffstudio.app.diff-studio", "@realizes:diffstudio.model.bioreactor"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens the \"Bioreactor\" model of the Diff Studio library", () => openLibraryModel(page, "Bioreactor"));
    await run.scenario("The model arrives with its table and its inputs", async () => {
      await session.step(21, "Then the \"Bioreactor\" view should be current", () => viewIsCurrent(page, "Bioreactor"));
      await session.step(22, "And grid should be visible", () => shouldBe(page, el("grid"), "visible"));
      await session.step(23, "And \"Process mode\" input should be visible", () => shouldBe(page, el("\"Process mode\" input"), "visible"));
      await session.step(24, "And \"switch at\" input should have value \"135\"", () => shouldHaveValue(page, el("\"switch at\" input"), "135"));
      await session.step(25, "And FFox input should have value \"0.20\"", () => shouldHaveValue(page, el("FFox input"), "0.20"));
    });
    await run.scenario("The line chart offers the Multiaxis and Facet tabs", async () => {
      await session.step(28, "Then Multiaxis tab should be visible", () => shouldBe(page, el("Multiaxis tab"), "visible"));
      await session.step(29, "And Facet tab should be visible", () => shouldBe(page, el("Facet tab"), "visible"));
      await session.step(30, "And Grid tab should be present", () => shouldBe(page, el("Grid tab"), "present"));
    });
    await run.scenario("The Facet tab draws the model as small multiples", async () => {
      await session.step(33, "When user clicks on Facet tab", () => clickOn(page, el("Facet tab")));
      await session.step(34, "Then Facet tab should be selected", () => shouldBe(page, el("Facet tab"), "selected"));
      await session.step(35, "And Multiaxis tab should not be selected", () => shouldNotBe(page, el("Multiaxis tab"), "selected"));
      await session.step(36, "And the canvases of open tableview should be painted in at least 10 colors", () => canvasColors(page, el("open tableview"), 10));
    });
    await run.scenario("Changing \"switch at\" redraws the table and the chart", async () => {
      await session.step(39, "When user clicks on Multiaxis tab", () => clickOn(page, el("Multiaxis tab")));
      await session.step(40, "And user takes a snapshot of line chart viewer", () => takeSnapshot(page, el("line chart viewer")));
      await session.step(41, "And user enters \"150\" into \"switch at\" input", () => enterInto(page, "150", el("\"switch at\" input")));
      await session.step(42, "Then \"switch at\" input should have value \"150\"", () => shouldHaveValue(page, el("\"switch at\" input"), "150"));
      await session.step(43, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
    });
    await run.scenario("The slider moves \"switch at\" and the chart follows", async () => {
      await session.step(46, "When user takes a snapshot of line chart viewer", () => takeSnapshot(page, el("line chart viewer")));
      await session.step(47, "And user drags the slider of \"switch at\" input to 100", () => dragSliderTo(page, el("\"switch at\" input"), 100));
      await session.step(48, "Then \"switch at\" input should have a value between 95 and 105", () => shouldHaveValueBetween(page, el("\"switch at\" input"), 95, 105));
      await session.step(49, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
    });
    await run.scenario("Process mode cascades into the parameters below it", async () => {
      await session.step(52, "When user takes a snapshot of line chart viewer", () => takeSnapshot(page, el("line chart viewer")));
      await session.step(53, "And user selects \"Mode 1\" in \"Process mode\" input", () => selectIn(page, "Mode 1", el("\"Process mode\" input")));
      await session.step(54, "Then \"Process mode\" input should have value \"Mode 1\"", () => shouldHaveValue(page, el("\"Process mode\" input"), "Mode 1"));
      await session.step(55, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
