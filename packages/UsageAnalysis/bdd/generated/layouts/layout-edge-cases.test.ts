/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/layouts/layout-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.layouts]
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
import {clickOn, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {closeAllViews, contextPanelOpen, contextPanelShows, layoutsDeleted, openDatasetRowsAs, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Layout edge cases: the Layouts pane suggests a layout only for a table that has its columns", () => {
  const session = feature(test, "features/layouts/layout-edge-cases.feature", import.meta.url);
  test("A layout plotting HEIGHT and WEIGHT is not suggested for a table without those columns", {tag: ["@realizes:views.layouts", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the layouts named \"bdd-edge-layout-{time}\" are deleted when the feature ends", () => layoutsDeleted(page, session.text("bdd-edge-layout-{time}")));
    await session.step(21, "And the layouts named \"bdd-edge-cars-{time}\" are deleted when the feature ends", () => layoutsDeleted(page, session.text("bdd-edge-cars-{time}")));
    await session.step(22, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(23, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(27, "Given user opens demog dataset keeping the first 100 rows as \"bdd-edge-layout-{time}\"", () => openDatasetRowsAs(page, ds("demog"), 100, session.text("bdd-edge-layout-{time}")));
    await session.step(28, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","HEIGHT"],["yColumnName","WEIGHT"]]), [["xColumnName","HEIGHT"],["yColumnName","WEIGHT"]]);
    await session.step(31, "And Layouts accordion header in toolbox is expanded", () => isExpanded(page, el("Layouts accordion header in toolbox")));
    await session.step(32, "When user clicks on the \"header AGE\" area of grid", () => clickArea(page, "header AGE", el("grid")));
    await session.step(33, "Then the context panel should show \"AGE\"", () => contextPanelShows(page, "AGE"));
    await session.step(34, "When user clicks on Save button in layouts pane", () => clickOn(page, el("Save button in layouts pane")));
    await session.step(35, "Then \"bdd-edge-layout-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"bdd-edge-layout-{time}\" layout card")), "visible"));
    await session.step(36, "When user closes all views", () => closeAllViews(page));
    await session.step(37, "Given user opens cars dataset keeping the first 30 rows as \"bdd-edge-cars-{time}\"", () => openDatasetRowsAs(page, ds("cars"), 30, session.text("bdd-edge-cars-{time}")));
    await session.step(38, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","peak.rpm"],["yColumnName","city.mpg"]]), [["xColumnName","peak.rpm"],["yColumnName","city.mpg"]]);
    await session.step(41, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(42, "And Layouts accordion header in toolbox is expanded", () => isExpanded(page, el("Layouts accordion header in toolbox")));
    await session.step(43, "When user clicks on the \"header model\" area of grid", () => clickArea(page, "header model", el("grid")));
    await session.step(44, "Then the context panel should show \"model\"", () => contextPanelShows(page, "model"));
    await session.step(45, "When user clicks on Save button in layouts pane", () => clickOn(page, el("Save button in layouts pane")));
    await session.step(46, "Then \"bdd-edge-cars-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"bdd-edge-cars-{time}\" layout card")), "visible"));
    await session.step(47, "When user closes all views", () => closeAllViews(page));
    await session.step(48, "Given user opens cars dataset keeping the first 30 rows as \"bdd-edge-cars-{time}\"", () => openDatasetRowsAs(page, ds("cars"), 30, session.text("bdd-edge-cars-{time}")));
    await session.step(49, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(50, "And Layouts accordion header in toolbox is expanded", () => isExpanded(page, el("Layouts accordion header in toolbox")));
    await session.step(51, "Then \"bdd-edge-cars-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"bdd-edge-cars-{time}\" layout card")), "visible"));
    await session.step(52, "And \"bdd-edge-layout-{time}\" layout card should be absent", () => shouldBe(page, el(session.text("\"bdd-edge-layout-{time}\" layout card")), "absent"));
    await session.step(53, "When user closes all views", () => closeAllViews(page));
    await session.step(54, "Given user opens demog dataset keeping the first 100 rows as \"bdd-edge-layout-{time}\"", () => openDatasetRowsAs(page, ds("demog"), 100, session.text("bdd-edge-layout-{time}")));
    await session.step(55, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(56, "And Layouts accordion header in toolbox is expanded", () => isExpanded(page, el("Layouts accordion header in toolbox")));
    await session.step(57, "Then \"bdd-edge-layout-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"bdd-edge-layout-{time}\" layout card")), "visible"));
    await session.step(58, "And \"bdd-edge-cars-{time}\" layout card should be absent", () => shouldBe(page, el(session.text("\"bdd-edge-cars-{time}\" layout card")), "absent"));
    await session.step(59, "And no errors should have been logged", () => noErrors(page));
  });
});
