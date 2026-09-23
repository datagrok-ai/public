/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/browse-apps-and-dashboards.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.browse]
--- */
import {test} from '@playwright/test';
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
import {clickOn, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, closeAllViews, openDataset, openProject, saveAsProject, urlShouldContain, viewHoldsViewers, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Apps and Dashboards sections of the Browse tree", () => {
  const session = feature(test, "features/browse/browse-apps-and-dashboards.feature", import.meta.url);
  test("The Apps section lists the installed applications", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(34, "Given user is logged in", () => loggedIn(page));
    await session.step(35, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(38, "Given Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(39, "Then Tutorials tree node inside browse tree should be visible", () => shouldBe(page, el("Tutorials tree node inside browse tree"), "visible"));
    await session.step(40, "And Misc tree node inside browse tree should be visible", () => shouldBe(page, el("Misc tree node inside browse tree"), "visible"));
    await session.step(41, "And no errors should have been logged", () => noErrors(page));
    await session.step(42, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("An application opens from the tree and becomes the current object", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(34, "Given user is logged in", () => loggedIn(page));
    await session.step(35, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(45, "Given Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(46, "When user clicks on Tutorials tree node inside browse tree", () => clickOn(page, el("Tutorials tree node inside browse tree")));
    await session.step(47, "Then the \"Tutorials\" view should be current", () => viewIsCurrent(page, "Tutorials"));
    await session.step(48, "And no errors should have been logged", () => noErrors(page));
    await session.step(49, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Model Hub opens from the Compute group", {tag: ["@browse", "@realizes:views.browse", "@compute"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(34, "Given user is logged in", () => loggedIn(page));
    await session.step(35, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(55, "Given Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(56, "And Apps---Compute tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Compute tree node inside browse tree")));
    await session.step(57, "When user clicks on Apps---Compute---Model-Hub tree node inside browse tree", () => clickOn(page, el("Apps---Compute---Model-Hub tree node inside browse tree")));
    await session.step(58, "Then the \"Model Hub\" view should be current", () => viewIsCurrent(page, "Model Hub"));
    await session.step(59, "And no errors should have been logged", () => noErrors(page));
    await session.step(60, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Dashboards node opens the list of projects", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(34, "Given user is logged in", () => loggedIn(page));
    await session.step(35, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(63, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(64, "Then the \"Projects\" view should be current", () => viewIsCurrent(page, "Projects"));
    await session.step(65, "And the page address should contain \"/projects\"", () => urlShouldContain(page, "/projects"));
    await session.step(68, "And card of gallery should be visible", () => shouldBe(page, el("card of gallery"), "visible"));
    await session.step(69, "And no errors should have been logged", () => noErrors(page));
    await session.step(70, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A dashboard saved from a table view opens again with its viewers", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(34, "Given user is logged in", () => loggedIn(page));
    await session.step(35, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(73, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(74, "And user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(75, "When user saves the current view as project \"BDD-Browse-Dash-{run}\"", () => saveAsProject(page, session.text("BDD-Browse-Dash-{run}")));
    await session.step(76, "And user closes all views", () => closeAllViews(page));
    await session.step(77, "And user opens the \"BDD-Browse-Dash-{run}\" project", () => openProject(page, session.text("BDD-Browse-Dash-{run}")));
    await session.step(78, "Then the current view should hold at least 2 viewers", () => viewHoldsViewers(page, 2));
    await session.step(79, "And no errors should have been logged", () => noErrors(page));
    await session.step(80, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
