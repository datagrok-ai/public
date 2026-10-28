/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/layouts/layout-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.layouts, core#1179]
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
import {addViewerWith, clickArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Layout edge cases: the Layouts pane suggests a layout only for a table that has its columns", () => {
  const session = feature(test, "features/layouts/layout-edge-cases.feature", import.meta.url);
  test("A layout plotting HEIGHT and WEIGHT is not suggested for a table without those columns", {tag: ["@realizes:views.layouts", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the layouts named \"bdd-edge-layout-{time}\" are deleted when the feature ends", () => layoutsDeleted(page, session.text("bdd-edge-layout-{time}")));
    await session.step(23, "And the layouts named \"bdd-edge-cars-{time}\" are deleted when the feature ends", () => layoutsDeleted(page, session.text("bdd-edge-cars-{time}")));
    await session.step(24, "And the layouts named \"bdd-edge-grid-{time}\" are deleted when the feature ends", () => layoutsDeleted(page, session.text("bdd-edge-grid-{time}")));
    await session.step(25, "And the layouts named \"bdd-edge-carsb-{time}\" are deleted when the feature ends", () => layoutsDeleted(page, session.text("bdd-edge-carsb-{time}")));
    await session.step(26, "And the layouts named \"bdd-edge-name-{time}\" are deleted when the feature ends", () => layoutsDeleted(page, session.text("bdd-edge-name-{time}")));
    await session.step(27, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(28, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(32, "Given user opens demog dataset keeping the first 100 rows as \"bdd-edge-layout-{time}\"", () => openDatasetRowsAs(page, ds("demog"), 100, session.text("bdd-edge-layout-{time}")));
    await session.step(33, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","HEIGHT"],["yColumnName","WEIGHT"]]), [["xColumnName","HEIGHT"],["yColumnName","WEIGHT"]]);
    await session.step(36, "And Layouts accordion header in toolbox is expanded", () => isExpanded(page, el("Layouts accordion header in toolbox")));
    await session.step(37, "When user clicks on the \"header AGE\" area of grid", () => clickArea(page, "header AGE", el("grid")));
    await session.step(38, "Then the context panel should show \"AGE\"", () => contextPanelShows(page, "AGE"));
    await session.step(39, "When user clicks on Save button in layouts pane", () => clickOn(page, el("Save button in layouts pane")));
    await session.step(40, "Then \"bdd-edge-layout-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"bdd-edge-layout-{time}\" layout card")), "visible"));
    await session.step(41, "When user closes all views", () => closeAllViews(page));
    await session.step(42, "Given user opens cars dataset keeping the first 30 rows as \"bdd-edge-cars-{time}\"", () => openDatasetRowsAs(page, ds("cars"), 30, session.text("bdd-edge-cars-{time}")));
    await session.step(43, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","peak.rpm"],["yColumnName","city.mpg"]]), [["xColumnName","peak.rpm"],["yColumnName","city.mpg"]]);
    await session.step(46, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(47, "And Layouts accordion header in toolbox is expanded", () => isExpanded(page, el("Layouts accordion header in toolbox")));
    await session.step(48, "When user clicks on the \"header model\" area of grid", () => clickArea(page, "header model", el("grid")));
    await session.step(49, "Then the context panel should show \"model\"", () => contextPanelShows(page, "model"));
    await session.step(50, "When user clicks on Save button in layouts pane", () => clickOn(page, el("Save button in layouts pane")));
    await session.step(51, "Then \"bdd-edge-cars-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"bdd-edge-cars-{time}\" layout card")), "visible"));
    await session.step(52, "When user closes all views", () => closeAllViews(page));
    await session.step(53, "Given user opens cars dataset keeping the first 30 rows as \"bdd-edge-cars-{time}\"", () => openDatasetRowsAs(page, ds("cars"), 30, session.text("bdd-edge-cars-{time}")));
    await session.step(54, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(55, "And Layouts accordion header in toolbox is expanded", () => isExpanded(page, el("Layouts accordion header in toolbox")));
    await session.step(56, "Then \"bdd-edge-cars-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"bdd-edge-cars-{time}\" layout card")), "visible"));
    await session.step(57, "And \"bdd-edge-layout-{time}\" layout card should be absent", () => shouldBe(page, el(session.text("\"bdd-edge-layout-{time}\" layout card")), "absent"));
    await session.step(58, "When user closes all views", () => closeAllViews(page));
    await session.step(59, "Given user opens demog dataset keeping the first 100 rows as \"bdd-edge-layout-{time}\"", () => openDatasetRowsAs(page, ds("demog"), 100, session.text("bdd-edge-layout-{time}")));
    await session.step(60, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(61, "And Layouts accordion header in toolbox is expanded", () => isExpanded(page, el("Layouts accordion header in toolbox")));
    await session.step(62, "Then \"bdd-edge-layout-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"bdd-edge-layout-{time}\" layout card")), "visible"));
    await session.step(63, "And \"bdd-edge-cars-{time}\" layout card should be absent", () => shouldBe(page, el(session.text("\"bdd-edge-cars-{time}\" layout card")), "absent"));
    await session.step(64, "And no errors should have been logged", () => noErrors(page));
  });
  test("A layout of a grid alone is suggested only for a table that has the grid's columns", {tag: ["@realizes:views.layouts", "@edge-cases", "@known-failure:85", "@realizes:core#1179"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the layouts named \"bdd-edge-layout-{time}\" are deleted when the feature ends", () => layoutsDeleted(page, session.text("bdd-edge-layout-{time}")));
    await session.step(23, "And the layouts named \"bdd-edge-cars-{time}\" are deleted when the feature ends", () => layoutsDeleted(page, session.text("bdd-edge-cars-{time}")));
    await session.step(24, "And the layouts named \"bdd-edge-grid-{time}\" are deleted when the feature ends", () => layoutsDeleted(page, session.text("bdd-edge-grid-{time}")));
    await session.step(25, "And the layouts named \"bdd-edge-carsb-{time}\" are deleted when the feature ends", () => layoutsDeleted(page, session.text("bdd-edge-carsb-{time}")));
    await session.step(26, "And the layouts named \"bdd-edge-name-{time}\" are deleted when the feature ends", () => layoutsDeleted(page, session.text("bdd-edge-name-{time}")));
    await session.step(27, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(28, "And the context panel is open", () => contextPanelOpen(page));
    await knownFailure(async () => {
      await session.step(71, "Given user opens demog dataset keeping the first 100 rows as \"bdd-edge-grid-{time}\"", () => openDatasetRowsAs(page, ds("demog"), 100, session.text("bdd-edge-grid-{time}")));
      await session.step(72, "And Layouts accordion header in toolbox is expanded", () => isExpanded(page, el("Layouts accordion header in toolbox")));
      await session.step(73, "When user clicks on the \"header AGE\" area of grid", () => clickArea(page, "header AGE", el("grid")));
      await session.step(74, "Then the context panel should show \"AGE\"", () => contextPanelShows(page, "AGE"));
      await session.step(75, "When user clicks on Save button in layouts pane", () => clickOn(page, el("Save button in layouts pane")));
      await session.step(76, "Then \"bdd-edge-grid-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"bdd-edge-grid-{time}\" layout card")), "visible"));
      await session.step(77, "When user closes all views", () => closeAllViews(page));
      await session.step(78, "Given user opens cars dataset keeping the first 30 rows as \"bdd-edge-carsb-{time}\"", () => openDatasetRowsAs(page, ds("cars"), 30, session.text("bdd-edge-carsb-{time}")));
      await session.step(79, "And the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(80, "And Layouts accordion header in toolbox is expanded", () => isExpanded(page, el("Layouts accordion header in toolbox")));
      await session.step(81, "When user clicks on the \"header model\" area of grid", () => clickArea(page, "header model", el("grid")));
      await session.step(82, "Then the context panel should show \"model\"", () => contextPanelShows(page, "model"));
      await session.step(83, "When user clicks on Save button in layouts pane", () => clickOn(page, el("Save button in layouts pane")));
      await session.step(84, "Then \"bdd-edge-carsb-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"bdd-edge-carsb-{time}\" layout card")), "visible"));
      await session.step(85, "And \"bdd-edge-grid-{time}\" layout card should be absent", () => shouldBe(page, el(session.text("\"bdd-edge-grid-{time}\" layout card")), "absent"));
      await session.step(86, "When user closes all views", () => closeAllViews(page));
      await session.step(87, "Given user opens demog dataset keeping the first 100 rows as \"bdd-edge-grid-{time}\"", () => openDatasetRowsAs(page, ds("demog"), 100, session.text("bdd-edge-grid-{time}")));
      await session.step(88, "And the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(89, "And Layouts accordion header in toolbox is expanded", () => isExpanded(page, el("Layouts accordion header in toolbox")));
      await session.step(90, "Then \"bdd-edge-grid-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"bdd-edge-grid-{time}\" layout card")), "visible"));
      await session.step(91, "And no errors should have been logged", () => noErrors(page));
    }, 85);
  });
  test("A layout renamed to an empty name in the context panel keeps its name", {tag: ["@realizes:views.layouts", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the layouts named \"bdd-edge-layout-{time}\" are deleted when the feature ends", () => layoutsDeleted(page, session.text("bdd-edge-layout-{time}")));
    await session.step(23, "And the layouts named \"bdd-edge-cars-{time}\" are deleted when the feature ends", () => layoutsDeleted(page, session.text("bdd-edge-cars-{time}")));
    await session.step(24, "And the layouts named \"bdd-edge-grid-{time}\" are deleted when the feature ends", () => layoutsDeleted(page, session.text("bdd-edge-grid-{time}")));
    await session.step(25, "And the layouts named \"bdd-edge-carsb-{time}\" are deleted when the feature ends", () => layoutsDeleted(page, session.text("bdd-edge-carsb-{time}")));
    await session.step(26, "And the layouts named \"bdd-edge-name-{time}\" are deleted when the feature ends", () => layoutsDeleted(page, session.text("bdd-edge-name-{time}")));
    await session.step(27, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(28, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(95, "Given user opens demog dataset keeping the first 100 rows as \"bdd-edge-name-{time}\"", () => openDatasetRowsAs(page, ds("demog"), 100, session.text("bdd-edge-name-{time}")));
    await session.step(96, "And Layouts accordion header in toolbox is expanded", () => isExpanded(page, el("Layouts accordion header in toolbox")));
    await session.step(97, "When user clicks on the \"header AGE\" area of grid", () => clickArea(page, "header AGE", el("grid")));
    await session.step(98, "Then the context panel should show \"AGE\"", () => contextPanelShows(page, "AGE"));
    await session.step(99, "When user clicks on Save button in layouts pane", () => clickOn(page, el("Save button in layouts pane")));
    await session.step(100, "Then \"bdd-edge-name-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"bdd-edge-name-{time}\" layout card")), "visible"));
    await session.step(101, "And the context panel should show \"bdd-edge-name-{time}\"", () => contextPanelShows(page, session.text("bdd-edge-name-{time}")));
    await session.step(102, "When user clears Name details field in context panel", () => clearField(page, el("Name details field in context panel")));
    await session.step(103, "And user presses Enter", () => pressKey(page, "Enter"));
    await session.step(104, "Then \"bdd-edge-name-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"bdd-edge-name-{time}\" layout card")), "visible"));
    await session.step(105, "And value of Name details field in context panel should have text \"bdd-edge-name-{time}\"", () => shouldHaveText(page, el("value of Name details field in context panel"), session.text("bdd-edge-name-{time}")));
    await session.step(106, "And no errors should have been logged", () => noErrors(page));
  });
});
