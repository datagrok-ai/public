/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/layouts/layouts-known-failures.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.layouts, core#1252]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/file-formats.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, isExpanded, pressKey, shouldBe, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {closeAllViews, contextPanelOpen, contextPanelShows, layoutsDeleted, openDatasetRowsAs, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("The known failures of layouts", () => {
  const session = feature(test, "features/layouts/layouts-known-failures.feature", import.meta.url);
  test("A layout renamed to an empty name in the context panel keeps its name", {tag: ["@realizes:views.layouts", "@edge-cases", "@known-failure:36", "@realizes:core#1252"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the context panel is open", () => contextPanelOpen(page));
    await knownFailure(async () => {
      await session.step(20, "Given the layouts named \"bdd-edge-name-{time}\" are deleted when the feature ends", () => layoutsDeleted(page, session.text("bdd-edge-name-{time}")));
      await session.step(21, "And user opens demog-1000 dataset keeping the first 1000 rows as \"bdd-edge-name-{time}\"", () => openDatasetRowsAs(page, ds("demog-1000"), 1000, session.text("bdd-edge-name-{time}")));
      await session.step(22, "And the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(23, "And Layouts accordion header in toolbox is expanded", () => isExpanded(page, el("Layouts accordion header in toolbox")));
      await session.step(24, "When user clicks on the \"header AGE\" area of grid", () => clickArea(page, "header AGE", el("grid")));
      await session.step(25, "Then the context panel should show \"AGE\"", () => contextPanelShows(page, "AGE"));
      await session.step(26, "When user clicks on Save button in layouts pane", () => clickOn(page, el("Save button in layouts pane")));
      await session.step(27, "Then \"bdd-edge-name-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"bdd-edge-name-{time}\" layout card")), "visible"));
      await session.step(28, "And the context panel should show \"bdd-edge-name-{time}\"", () => contextPanelShows(page, session.text("bdd-edge-name-{time}")));
      await session.step(29, "When user clears Name details field in context panel", () => clearField(page, el("Name details field in context panel")));
      await session.step(30, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(31, "Then value of Name details field in context panel should have text \"bdd-edge-name-{time}\"", () => shouldHaveText(page, el("value of Name details field in context panel"), session.text("bdd-edge-name-{time}")));
      await session.step(32, "When user closes all views", () => closeAllViews(page));
      await session.step(33, "Given user opens demog-1000 dataset keeping the first 1000 rows as \"bdd-edge-name-{time}\"", () => openDatasetRowsAs(page, ds("demog-1000"), 1000, session.text("bdd-edge-name-{time}")));
      await session.step(34, "And the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(35, "And Layouts accordion header in toolbox is expanded", () => isExpanded(page, el("Layouts accordion header in toolbox")));
      await session.step(36, "Then \"bdd-edge-name-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"bdd-edge-name-{time}\" layout card")), "visible"));
    }, 36);
  });
});
