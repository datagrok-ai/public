/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/browse-apps-matrix.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.browse]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, urlShouldContain, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Every application a stand runs in the browser opens from Browse > Apps with its own content", () => {
  const session = feature(test, "features/browse/browse-apps-matrix.feature", import.meta.url);
  test("Apps---Chem---Reactions---Reaction-Enumerator opens from the tree with its own content [group=Chem, parent=Apps---Chem---Reactions, node=Apps---Chem---Reactions---Reaction-Enumerator, view=Reaction Enumerator, address=/apps/Chem/ReactionEnumerator, content=\"Number of steps\" input, more=\"Next: Reactions\" button]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(24, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(27, "Given Apps---Chem tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Chem tree node inside browse tree")));
    await session.step(28, "And Apps---Chem---Reactions tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Chem---Reactions tree node inside browse tree")));
    await session.step(29, "When user clicks on Apps---Chem---Reactions---Reaction-Enumerator tree node inside browse tree", () => clickOn(page, el("Apps---Chem---Reactions---Reaction-Enumerator tree node inside browse tree")));
    await session.step(30, "Then the \"Reaction Enumerator\" view should be current", () => viewIsCurrent(page, "Reaction Enumerator"));
    await session.step(31, "And the page address should contain \"/apps/Chem/ReactionEnumerator\"", () => urlShouldContain(page, "/apps/Chem/ReactionEnumerator"));
    await session.step(32, "And \"Number of steps\" input should be visible", () => shouldBe(page, el("\"Number of steps\" input"), "visible"));
    await session.step(33, "And \"Next: Reactions\" button should be visible", () => shouldBe(page, el("\"Next: Reactions\" button"), "visible"));
    await session.step(34, "And no errors should have been logged", () => noErrors(page));
    await session.step(35, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Apps---Chem---Reactions---Transformation-Reactions opens from the tree with its own content [group=Chem, parent=Apps---Chem---Reactions, node=Apps---Chem---Reactions---Transformation-Reactions, view=Transformation Reactions, address=/apps/Chem/TransformationReactions, content=\"Molecules\" input, more=\"Run Reaction\" button]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(24, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(27, "Given Apps---Chem tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Chem tree node inside browse tree")));
    await session.step(28, "And Apps---Chem---Reactions tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Chem---Reactions tree node inside browse tree")));
    await session.step(29, "When user clicks on Apps---Chem---Reactions---Transformation-Reactions tree node inside browse tree", () => clickOn(page, el("Apps---Chem---Reactions---Transformation-Reactions tree node inside browse tree")));
    await session.step(30, "Then the \"Transformation Reactions\" view should be current", () => viewIsCurrent(page, "Transformation Reactions"));
    await session.step(31, "And the page address should contain \"/apps/Chem/TransformationReactions\"", () => urlShouldContain(page, "/apps/Chem/TransformationReactions"));
    await session.step(32, "And \"Molecules\" input should be visible", () => shouldBe(page, el("\"Molecules\" input"), "visible"));
    await session.step(33, "And \"Run Reaction\" button should be visible", () => shouldBe(page, el("\"Run Reaction\" button"), "visible"));
    await session.step(34, "And no errors should have been logged", () => noErrors(page));
    await session.step(35, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Apps---Chem---Reactions---Two-Component-Reactions opens from the tree with its own content [group=Chem, parent=Apps---Chem---Reactions, node=Apps---Chem---Reactions---Two-Component-Reactions, view=Two-Component Reactions, address=/apps/Chem/TwoComponentReactions, content=\"Reactant 1\" input, more=\"Run Reaction\" button]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(24, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(27, "Given Apps---Chem tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Chem tree node inside browse tree")));
    await session.step(28, "And Apps---Chem---Reactions tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Chem---Reactions tree node inside browse tree")));
    await session.step(29, "When user clicks on Apps---Chem---Reactions---Two-Component-Reactions tree node inside browse tree", () => clickOn(page, el("Apps---Chem---Reactions---Two-Component-Reactions tree node inside browse tree")));
    await session.step(30, "Then the \"Two-Component Reactions\" view should be current", () => viewIsCurrent(page, "Two-Component Reactions"));
    await session.step(31, "And the page address should contain \"/apps/Chem/TwoComponentReactions\"", () => urlShouldContain(page, "/apps/Chem/TwoComponentReactions"));
    await session.step(32, "And \"Reactant 1\" input should be visible", () => shouldBe(page, el("\"Reactant 1\" input"), "visible"));
    await session.step(33, "And \"Run Reaction\" button should be visible", () => shouldBe(page, el("\"Run Reaction\" button"), "visible"));
    await session.step(34, "And no errors should have been logged", () => noErrors(page));
    await session.step(35, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Apps---Dev---Reports-Browser opens from the tree with its own content [group=Dev, parent=Apps---Dev, node=Apps---Dev---Reports-Browser, view=Reports, address=/apps/U2demo/ReportsBrowser, content=first item in list, more=first row actions]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(24, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(27, "Given Apps---Dev tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Dev tree node inside browse tree")));
    await session.step(28, "And Apps---Dev tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Dev tree node inside browse tree")));
    await session.step(29, "When user clicks on Apps---Dev---Reports-Browser tree node inside browse tree", () => clickOn(page, el("Apps---Dev---Reports-Browser tree node inside browse tree")));
    await session.step(30, "Then the \"Reports\" view should be current", () => viewIsCurrent(page, "Reports"));
    await session.step(31, "And the page address should contain \"/apps/U2demo/ReportsBrowser\"", () => urlShouldContain(page, "/apps/U2demo/ReportsBrowser"));
    await session.step(32, "And first item in list should be visible", () => shouldBe(page, el("first item in list"), "visible"));
    await session.step(33, "And first row actions should be visible", () => shouldBe(page, el("first row actions"), "visible"));
    await session.step(34, "And no errors should have been logged", () => noErrors(page));
    await session.step(35, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Apps---Dev---U2-Designer opens from the tree with its own content [group=Dev, parent=Apps---Dev, node=Apps---Dev---U2-Designer, view=U2 Designer, address=/apps/U2demo/U2Designer, content=\"nameInput\" element, more=\"Save\" button]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(24, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(27, "Given Apps---Dev tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Dev tree node inside browse tree")));
    await session.step(28, "And Apps---Dev tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Dev tree node inside browse tree")));
    await session.step(29, "When user clicks on Apps---Dev---U2-Designer tree node inside browse tree", () => clickOn(page, el("Apps---Dev---U2-Designer tree node inside browse tree")));
    await session.step(30, "Then the \"U2 Designer\" view should be current", () => viewIsCurrent(page, "U2 Designer"));
    await session.step(31, "And the page address should contain \"/apps/U2demo/U2Designer\"", () => urlShouldContain(page, "/apps/U2demo/U2Designer"));
    await session.step(32, "And \"nameInput\" element should be visible", () => shouldBe(page, el("\"nameInput\" element"), "visible"));
    await session.step(33, "And \"Save\" button should be visible", () => shouldBe(page, el("\"Save\" button"), "visible"));
    await session.step(34, "And no errors should have been logged", () => noErrors(page));
    await session.step(35, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A double click keeps Apps---Chem---Reactions---Transformation-Reactions open through the next click in the tree [group=Chem, parent=Apps---Chem---Reactions, node=Apps---Chem---Reactions---Transformation-Reactions, view=Transformation Reactions]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(24, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(46, "Given Apps---Chem tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Chem tree node inside browse tree")));
    await session.step(47, "And Apps---Chem---Reactions tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Chem---Reactions tree node inside browse tree")));
    await session.step(48, "When user double-clicks on Apps---Chem---Reactions---Transformation-Reactions tree node inside browse tree", () => doubleClickOn(page, el("Apps---Chem---Reactions---Transformation-Reactions tree node inside browse tree")));
    await session.step(49, "Then the \"Transformation Reactions\" view should be current", () => viewIsCurrent(page, "Transformation Reactions"));
    await session.step(50, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(51, "Then Projects view should be visible", () => shouldBe(page, el("Projects view"), "visible"));
    await session.step(52, "And \"Transformation Reactions\" view should be present", () => shouldBe(page, el("\"Transformation Reactions\" view"), "present"));
    await session.step(53, "And no errors should have been logged", () => noErrors(page));
    await session.step(54, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A double click keeps Apps---Dev---U2-Designer open through the next click in the tree [group=Dev, parent=Apps---Dev, node=Apps---Dev---U2-Designer, view=U2 Designer]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(24, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(46, "Given Apps---Dev tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Dev tree node inside browse tree")));
    await session.step(47, "And Apps---Dev tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Dev tree node inside browse tree")));
    await session.step(48, "When user double-clicks on Apps---Dev---U2-Designer tree node inside browse tree", () => doubleClickOn(page, el("Apps---Dev---U2-Designer tree node inside browse tree")));
    await session.step(49, "Then the \"U2 Designer\" view should be current", () => viewIsCurrent(page, "U2 Designer"));
    await session.step(50, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(51, "Then Projects view should be visible", () => shouldBe(page, el("Projects view"), "visible"));
    await session.step(52, "And \"U2 Designer\" view should be present", () => shouldBe(page, el("\"U2 Designer\" view"), "present"));
    await session.step(53, "And no errors should have been logged", () => noErrors(page));
    await session.step(54, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
