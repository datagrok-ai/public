/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tooltips/tooltip-context-menu.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tooltips]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/grid.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, closeContextMenu, menuLists, noErrors, openContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Tooltip group of a viewer's context menu", () => {
  const session = feature(test, "features/viewers/tooltips/tooltip-context-menu.feature", import.meta.url);
  test("The Tooltip group of a viewer's context menu", {tag: ["@journey", "@viewers", "@realizes:viewers.tooltips"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a histogram viewer", () => addViewer(page, "histogram"));
    await session.step(16, "And user adds a line chart viewer", () => addViewer(page, "line chart"));
    await session.step(17, "And user adds a bar chart viewer", () => addViewer(page, "bar chart"));
    await session.step(18, "And user adds a trellis plot viewer", () => addViewer(page, "trellis plot"));
    await run.scenario("The context menu of the grid offers the four Tooltip actions [viewer=grid]", async () => {
      await session.step(21, "When user opens the context menu of grid", () => openContextMenu(page, el("grid")));
      await session.step(22, "Then the open menu should list \"Tooltip > Hide\"", () => menuLists(page, "Tooltip > Hide"));
      await session.step(23, "And the open menu should list \"Tooltip > Edit...\"", () => menuLists(page, "Tooltip > Edit..."));
      await session.step(24, "And the open menu should list \"Tooltip > Use as Group Tooltip\"", () => menuLists(page, "Tooltip > Use as Group Tooltip"));
      await session.step(25, "And the open menu should list \"Tooltip > Remove Group Tooltip\"", () => menuLists(page, "Tooltip > Remove Group Tooltip"));
      await session.step(26, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(27, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The context menu of the histogram viewer offers the four Tooltip actions [viewer=histogram viewer]", async () => {
      await session.step(21, "When user opens the context menu of histogram viewer", () => openContextMenu(page, el("histogram viewer")));
      await session.step(22, "Then the open menu should list \"Tooltip > Hide\"", () => menuLists(page, "Tooltip > Hide"));
      await session.step(23, "And the open menu should list \"Tooltip > Edit...\"", () => menuLists(page, "Tooltip > Edit..."));
      await session.step(24, "And the open menu should list \"Tooltip > Use as Group Tooltip\"", () => menuLists(page, "Tooltip > Use as Group Tooltip"));
      await session.step(25, "And the open menu should list \"Tooltip > Remove Group Tooltip\"", () => menuLists(page, "Tooltip > Remove Group Tooltip"));
      await session.step(26, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(27, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The context menu of the line chart viewer offers the four Tooltip actions [viewer=line chart viewer]", async () => {
      await session.step(21, "When user opens the context menu of line chart viewer", () => openContextMenu(page, el("line chart viewer")));
      await session.step(22, "Then the open menu should list \"Tooltip > Hide\"", () => menuLists(page, "Tooltip > Hide"));
      await session.step(23, "And the open menu should list \"Tooltip > Edit...\"", () => menuLists(page, "Tooltip > Edit..."));
      await session.step(24, "And the open menu should list \"Tooltip > Use as Group Tooltip\"", () => menuLists(page, "Tooltip > Use as Group Tooltip"));
      await session.step(25, "And the open menu should list \"Tooltip > Remove Group Tooltip\"", () => menuLists(page, "Tooltip > Remove Group Tooltip"));
      await session.step(26, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(27, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The context menu of the bar chart viewer offers the four Tooltip actions [viewer=bar chart viewer]", async () => {
      await session.step(21, "When user opens the context menu of bar chart viewer", () => openContextMenu(page, el("bar chart viewer")));
      await session.step(22, "Then the open menu should list \"Tooltip > Hide\"", () => menuLists(page, "Tooltip > Hide"));
      await session.step(23, "And the open menu should list \"Tooltip > Edit...\"", () => menuLists(page, "Tooltip > Edit..."));
      await session.step(24, "And the open menu should list \"Tooltip > Use as Group Tooltip\"", () => menuLists(page, "Tooltip > Use as Group Tooltip"));
      await session.step(25, "And the open menu should list \"Tooltip > Remove Group Tooltip\"", () => menuLists(page, "Tooltip > Remove Group Tooltip"));
      await session.step(26, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(27, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The context menu of the trellis plot viewer offers the four Tooltip actions [viewer=trellis plot viewer]", async () => {
      await session.step(21, "When user opens the context menu of trellis plot viewer", () => openContextMenu(page, el("trellis plot viewer")));
      await session.step(22, "Then the open menu should list \"Tooltip > Hide\"", () => menuLists(page, "Tooltip > Hide"));
      await session.step(23, "And the open menu should list \"Tooltip > Edit...\"", () => menuLists(page, "Tooltip > Edit..."));
      await session.step(24, "And the open menu should list \"Tooltip > Use as Group Tooltip\"", () => menuLists(page, "Tooltip > Use as Group Tooltip"));
      await session.step(25, "And the open menu should list \"Tooltip > Remove Group Tooltip\"", () => menuLists(page, "Tooltip > Remove Group Tooltip"));
      await session.step(26, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(27, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
