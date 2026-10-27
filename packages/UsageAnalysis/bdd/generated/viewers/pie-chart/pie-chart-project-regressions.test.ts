/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pie-chart/pie-chart-project-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pie-chart, public#482]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/file-formats.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {closeAllViews, openDataset, openProjectWithTable, ownProjectGone, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, addViewerWith, legendLists, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {stackViewer, switchTab, viewersTabbed} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Pie chart regressions — the legend of a pie chart stacked as a second tab", () => {
  const session = feature(test, "features/viewers/pie-chart/pie-chart-project-regressions.feature", import.meta.url);
  test("A pie chart stacked as the second tab shows its legend after the project is reopened (public#482)", {tag: ["@viewers", "@realizes:viewers.pie-chart", "@realizes:public#482"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "Given the user's own project \"bdd-pie-tab-{time}\" is removed now and at feature end", () => ownProjectGone(page, session.text("bdd-pie-tab-{time}")));
    await session.step(16, "And user adds a bar chart viewer", () => addViewer(page, "bar chart"));
    await session.step(17, "And user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"],["Legend Visibility","Always"]]), [["Category","RACE"],["Legend Visibility","Always"]]);
    await session.step(20, "When user stacks pie chart viewer onto bar chart viewer as a tab", () => stackViewer(page, el("pie chart viewer"), el("bar chart viewer")));
    await session.step(21, "Then the viewers \"Bar chart, Pie chart\" should be the tabs of one panel", () => viewersTabbed(page, "Bar chart, Pie chart"));
    await session.step(22, "When user switches the tabbed panel to the \"Pie chart\" tab", () => switchTab(page, "Pie chart"));
    await session.step(23, "Then the legend of pie chart viewer should list 4 items", () => legendLists(page, el("pie chart viewer"), 4));
    await session.step(24, "When user saves the current view as project \"bdd-pie-tab-{time}\"", () => saveAsProject(page, session.text("bdd-pie-tab-{time}")));
    await session.step(25, "And user closes all views", () => closeAllViews(page));
    await session.step(26, "And user opens the \"bdd-pie-tab-{time}\" project and waits for its table", () => openProjectWithTable(page, session.text("bdd-pie-tab-{time}")));
    await session.step(27, "Then the viewers \"Bar chart, Pie chart\" should be the tabs of one panel", () => viewersTabbed(page, "Bar chart, Pie chart"));
    await session.step(28, "When user switches the tabbed panel to the \"Pie chart\" tab", () => switchTab(page, "Pie chart"));
    await session.step(29, "Then the legend of pie chart viewer should list 4 items", () => legendLists(page, el("pie chart viewer"), 4));
    await session.step(30, "And legend of pie chart viewer should be visible", () => shouldBe(page, el("legend of pie chart viewer"), "visible"));
    await session.step(31, "And no errors should have been logged", () => noErrors(page));
  });
});
