/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/browse-platform-and-databases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.browse]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, collapse, followingShouldBe, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Platform and Databases sections of the Browse tree", () => {
  const session = feature(test, "features/browse/browse-platform-and-databases.feature", import.meta.url);
  test("The Platform section lists what an administrator manages", {tag: ["@browse", "@realizes:views.browse", "@full-stand"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(29, "Given Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(30, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["Platform---Plugins tree node inside browse tree"],["Platform---Credentials tree node inside browse tree"],["Platform---Functions tree node inside browse tree"],["Platform---Users tree node inside browse tree"],["Platform---Groups tree node inside browse tree"],["Platform---Roles tree node inside browse tree"],["Platform---Predictive-models tree node inside browse tree"],["Platform---Dockers tree node inside browse tree"]]), [["Platform---Plugins tree node inside browse tree"],["Platform---Credentials tree node inside browse tree"],["Platform---Functions tree node inside browse tree"],["Platform---Users tree node inside browse tree"],["Platform---Groups tree node inside browse tree"],["Platform---Roles tree node inside browse tree"],["Platform---Predictive-models tree node inside browse tree"],["Platform---Dockers tree node inside browse tree"]]);
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
    await session.step(40, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A Platform node opens the view it is named after", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(43, "Given Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(44, "When user clicks on Platform---Users tree node inside browse tree", () => clickOn(page, el("Platform---Users tree node inside browse tree")));
    await session.step(45, "Then the \"Users\" view should be current", () => viewIsCurrent(page, "Users"));
    await session.step(46, "When user clicks on Platform---Groups tree node inside browse tree", () => clickOn(page, el("Platform---Groups tree node inside browse tree")));
    await session.step(47, "Then the \"Groups\" view should be current", () => viewIsCurrent(page, "Groups"));
    await session.step(48, "And no errors should have been logged", () => noErrors(page));
    await session.step(49, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Databases section lists the connected providers", {tag: ["@browse", "@realizes:views.browse", "@full-stand"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(53, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(54, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["Databases---Postgres tree node inside browse tree"],["Databases---MySQL tree node inside browse tree"],["Databases---Oracle tree node inside browse tree"],["Databases---MariaDB tree node inside browse tree"],["Databases---MS-SQL tree node inside browse tree"]]), [["Databases---Postgres tree node inside browse tree"],["Databases---MySQL tree node inside browse tree"],["Databases---Oracle tree node inside browse tree"],["Databases---MariaDB tree node inside browse tree"],["Databases---MS-SQL tree node inside browse tree"]]);
    await session.step(60, "And no errors should have been logged", () => noErrors(page));
    await session.step(61, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A provider opens down to its connections", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(64, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(65, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(66, "Then Databases---Postgres---Datagrok tree node inside browse tree should be visible", () => shouldBe(page, el("Databases---Postgres---Datagrok tree node inside browse tree"), "visible"));
    await session.step(67, "And Databases---Postgres---CHEMBL tree node inside browse tree should be visible", () => shouldBe(page, el("Databases---Postgres---CHEMBL tree node inside browse tree"), "visible"));
    await session.step(68, "When user collapses Databases---Postgres tree node inside browse tree", () => collapse(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(69, "Then Databases---Postgres---Datagrok tree node inside browse tree should be hidden", () => shouldBe(page, el("Databases---Postgres---Datagrok tree node inside browse tree"), "hidden"));
    await session.step(70, "And Databases---Postgres tree node inside browse tree should be visible", () => shouldBe(page, el("Databases---Postgres tree node inside browse tree"), "visible"));
    await session.step(71, "And no errors should have been logged", () => noErrors(page));
    await session.step(72, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
