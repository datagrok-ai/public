/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/viewers-docking.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.chrome]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/nx.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {closeAllViews, openDataset, viewHoldsViewers} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, loadLayout, noErrors, saveLayoutToServer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {dockBesideViewer, dockToViewEdge, dockedAtViewEdge, dockedBeside, dockedInCorner, notDockedAtViewEdge, notDockedBeside} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Docking viewers by their title bars", () => {
  const session = feature(test, "features/viewers/viewers-docking.feature", import.meta.url);
  test("A viewer dropped at the right edge of the view docks along the whole edge", {tag: ["@viewers", "@realizes:viewers.chrome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(25, "And user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(26, "And user adds a histogram viewer", () => addViewer(page, "histogram"));
    await session.step(27, "And user adds a bar chart viewer", () => addViewer(page, "bar chart"));
    await session.step(28, "Then the current view should hold at least 3 viewers", () => viewHoldsViewers(page, 3));
    await session.step(31, "Then histogram viewer should not be docked along the right edge of the view", () => notDockedAtViewEdge(page, el("histogram viewer"), "right"));
    await session.step(32, "When user docks histogram viewer to the right edge of the view", () => dockToViewEdge(page, el("histogram viewer"), "right"));
    await session.step(33, "Then histogram viewer should be docked along the right edge of the view", () => dockedAtViewEdge(page, el("histogram viewer"), "right"));
    await session.step(34, "And bar chart viewer should be visible", () => shouldBe(page, el("bar chart viewer"), "visible"));
    await session.step(35, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(36, "And no errors should have been logged", () => noErrors(page));
  });
  test("A viewer dropped at the bottom edge after the right edge was taken runs under all of them", {tag: ["@viewers", "@realizes:viewers.chrome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(25, "And user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(26, "And user adds a histogram viewer", () => addViewer(page, "histogram"));
    await session.step(27, "And user adds a bar chart viewer", () => addViewer(page, "bar chart"));
    await session.step(28, "Then the current view should hold at least 3 viewers", () => viewHoldsViewers(page, 3));
    await session.step(39, "When user docks histogram viewer to the right edge of the view", () => dockToViewEdge(page, el("histogram viewer"), "right"));
    await session.step(40, "Then histogram viewer should be docked along the right edge of the view", () => dockedAtViewEdge(page, el("histogram viewer"), "right"));
    await session.step(41, "And bar chart viewer should not be docked along the bottom edge of the view", () => notDockedAtViewEdge(page, el("bar chart viewer"), "bottom"));
    await session.step(42, "When user docks bar chart viewer to the bottom edge of the view", () => dockToViewEdge(page, el("bar chart viewer"), "bottom"));
    await session.step(43, "Then bar chart viewer should be docked along the bottom edge of the view", () => dockedAtViewEdge(page, el("bar chart viewer"), "bottom"));
    await session.step(44, "And histogram viewer should not be docked along the right edge of the view", () => notDockedAtViewEdge(page, el("histogram viewer"), "right"));
    await session.step(45, "And histogram viewer should be docked in the top right corner of the view", () => dockedInCorner(page, el("histogram viewer"), "top", "right"));
    await session.step(46, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(47, "And no errors should have been logged", () => noErrors(page));
  });
  test("A viewer dropped on a side of another viewer docks next to that viewer only", {tag: ["@viewers", "@realizes:viewers.chrome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(25, "And user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(26, "And user adds a histogram viewer", () => addViewer(page, "histogram"));
    await session.step(27, "And user adds a bar chart viewer", () => addViewer(page, "bar chart"));
    await session.step(28, "Then the current view should hold at least 3 viewers", () => viewHoldsViewers(page, 3));
    await session.step(50, "Then histogram viewer should not be docked above scatter plot viewer", () => notDockedBeside(page, el("histogram viewer"), "above", el("scatter plot viewer")));
    await session.step(51, "When user docks histogram viewer to the top side of scatter plot viewer", () => dockBesideViewer(page, el("histogram viewer"), "top", el("scatter plot viewer")));
    await session.step(52, "Then histogram viewer should be docked above scatter plot viewer", () => dockedBeside(page, el("histogram viewer"), "above", el("scatter plot viewer")));
    await session.step(53, "And bar chart viewer should not be docked left-of scatter plot viewer", () => notDockedBeside(page, el("bar chart viewer"), "left-of", el("scatter plot viewer")));
    await session.step(54, "When user docks bar chart viewer to the left side of scatter plot viewer", () => dockBesideViewer(page, el("bar chart viewer"), "left", el("scatter plot viewer")));
    await session.step(55, "Then bar chart viewer should be docked left-of scatter plot viewer", () => dockedBeside(page, el("bar chart viewer"), "left-of", el("scatter plot viewer")));
    await session.step(56, "And no errors should have been logged", () => noErrors(page));
  });
  test("A docking arrangement comes back from a layout applied to the table opened anew", {tag: ["@viewers", "@realizes:viewers.chrome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(25, "And user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(26, "And user adds a histogram viewer", () => addViewer(page, "histogram"));
    await session.step(27, "And user adds a bar chart viewer", () => addViewer(page, "bar chart"));
    await session.step(28, "Then the current view should hold at least 3 viewers", () => viewHoldsViewers(page, 3));
    await session.step(59, "When user docks histogram viewer to the right edge of the view", () => dockToViewEdge(page, el("histogram viewer"), "right"));
    await session.step(60, "Then bar chart viewer should not be docked below histogram viewer", () => notDockedBeside(page, el("bar chart viewer"), "below", el("histogram viewer")));
    await session.step(61, "When user docks bar chart viewer to the bottom side of histogram viewer", () => dockBesideViewer(page, el("bar chart viewer"), "bottom", el("histogram viewer")));
    await session.step(62, "Then bar chart viewer should be docked below histogram viewer", () => dockedBeside(page, el("bar chart viewer"), "below", el("histogram viewer")));
    await session.step(63, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
    await session.step(64, "And user closes all views", () => closeAllViews(page));
    await session.step(65, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(66, "And user loads the saved layout", () => loadLayout(page));
    await session.step(67, "Then scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(68, "And bar chart viewer should be docked below histogram viewer", () => dockedBeside(page, el("bar chart viewer"), "below", el("histogram viewer")));
    await session.step(69, "And no errors should have been logged", () => noErrors(page));
  });
});
