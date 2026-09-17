/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/browse-node-matrix.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.browse]
--- */
import {test} from '@playwright/test';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, collapse, expand, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, contextPanelOpen} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Every section of the Browse tree opens without an error", () => {
  const session = feature(test, "features/browse/browse-node-matrix.feature", import.meta.url);
  test("Clicking My stuff logs no error [node=My stuff]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(22, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(28, "When user clicks on first My stuff tree node inside browse tree", () => clickOn(page, el("first My stuff tree node inside browse tree")));
    await session.step(29, "Then first My stuff tree node inside browse tree should be selected", () => shouldBe(page, el("first My stuff tree node inside browse tree"), "selected"));
    await session.step(30, "And no errors should have been logged", () => noErrors(page));
    await session.step(31, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Clicking Spaces logs no error [node=Spaces]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(22, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(28, "When user clicks on first Spaces tree node inside browse tree", () => clickOn(page, el("first Spaces tree node inside browse tree")));
    await session.step(29, "Then first Spaces tree node inside browse tree should be selected", () => shouldBe(page, el("first Spaces tree node inside browse tree"), "selected"));
    await session.step(30, "And no errors should have been logged", () => noErrors(page));
    await session.step(31, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Clicking Apps logs no error [node=Apps]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(22, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(28, "When user clicks on first Apps tree node inside browse tree", () => clickOn(page, el("first Apps tree node inside browse tree")));
    await session.step(29, "Then first Apps tree node inside browse tree should be selected", () => shouldBe(page, el("first Apps tree node inside browse tree"), "selected"));
    await session.step(30, "And no errors should have been logged", () => noErrors(page));
    await session.step(31, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Clicking Files logs no error [node=Files]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(22, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(28, "When user clicks on first Files tree node inside browse tree", () => clickOn(page, el("first Files tree node inside browse tree")));
    await session.step(29, "Then first Files tree node inside browse tree should be selected", () => shouldBe(page, el("first Files tree node inside browse tree"), "selected"));
    await session.step(30, "And no errors should have been logged", () => noErrors(page));
    await session.step(31, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Clicking Dashboards logs no error [node=Dashboards]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(22, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(28, "When user clicks on first Dashboards tree node inside browse tree", () => clickOn(page, el("first Dashboards tree node inside browse tree")));
    await session.step(29, "Then first Dashboards tree node inside browse tree should be selected", () => shouldBe(page, el("first Dashboards tree node inside browse tree"), "selected"));
    await session.step(30, "And no errors should have been logged", () => noErrors(page));
    await session.step(31, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Clicking Databases logs no error [node=Databases]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(22, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(28, "When user clicks on first Databases tree node inside browse tree", () => clickOn(page, el("first Databases tree node inside browse tree")));
    await session.step(29, "Then first Databases tree node inside browse tree should be selected", () => shouldBe(page, el("first Databases tree node inside browse tree"), "selected"));
    await session.step(30, "And no errors should have been logged", () => noErrors(page));
    await session.step(31, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Clicking Platform logs no error [node=Platform]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(22, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(28, "When user clicks on first Platform tree node inside browse tree", () => clickOn(page, el("first Platform tree node inside browse tree")));
    await session.step(29, "Then first Platform tree node inside browse tree should be selected", () => shouldBe(page, el("first Platform tree node inside browse tree"), "selected"));
    await session.step(30, "And no errors should have been logged", () => noErrors(page));
    await session.step(31, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Opening the My stuff section logs no error [section=My stuff, child=My-stuff---Recent]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(22, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(49, "Given user collapses first My stuff tree node inside browse tree", () => collapse(page, el("first My stuff tree node inside browse tree")));
    await session.step(50, "And first My stuff tree node inside browse tree should be collapsed", () => shouldBe(page, el("first My stuff tree node inside browse tree"), "collapsed"));
    await session.step(51, "When user expands first My stuff tree node inside browse tree", () => expand(page, el("first My stuff tree node inside browse tree")));
    await session.step(52, "Then first My stuff tree node inside browse tree should be expanded", () => shouldBe(page, el("first My stuff tree node inside browse tree"), "expanded"));
    await session.step(53, "And My-stuff---Recent tree node inside browse tree should be visible", () => shouldBe(page, el("My-stuff---Recent tree node inside browse tree"), "visible"));
    await session.step(54, "And no errors should have been logged", () => noErrors(page));
    await session.step(55, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Opening the Apps section logs no error [section=Apps, child=Apps---Compute]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(22, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(49, "Given user collapses first Apps tree node inside browse tree", () => collapse(page, el("first Apps tree node inside browse tree")));
    await session.step(50, "And first Apps tree node inside browse tree should be collapsed", () => shouldBe(page, el("first Apps tree node inside browse tree"), "collapsed"));
    await session.step(51, "When user expands first Apps tree node inside browse tree", () => expand(page, el("first Apps tree node inside browse tree")));
    await session.step(52, "Then first Apps tree node inside browse tree should be expanded", () => shouldBe(page, el("first Apps tree node inside browse tree"), "expanded"));
    await session.step(53, "And Apps---Compute tree node inside browse tree should be visible", () => shouldBe(page, el("Apps---Compute tree node inside browse tree"), "visible"));
    await session.step(54, "And no errors should have been logged", () => noErrors(page));
    await session.step(55, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Opening the Files section logs no error [section=Files, child=Files---Demo]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(22, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(49, "Given user collapses first Files tree node inside browse tree", () => collapse(page, el("first Files tree node inside browse tree")));
    await session.step(50, "And first Files tree node inside browse tree should be collapsed", () => shouldBe(page, el("first Files tree node inside browse tree"), "collapsed"));
    await session.step(51, "When user expands first Files tree node inside browse tree", () => expand(page, el("first Files tree node inside browse tree")));
    await session.step(52, "Then first Files tree node inside browse tree should be expanded", () => shouldBe(page, el("first Files tree node inside browse tree"), "expanded"));
    await session.step(53, "And Files---Demo tree node inside browse tree should be visible", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "visible"));
    await session.step(54, "And no errors should have been logged", () => noErrors(page));
    await session.step(55, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Opening the Databases section logs no error [section=Databases, child=Databases---Postgres]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(22, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(49, "Given user collapses first Databases tree node inside browse tree", () => collapse(page, el("first Databases tree node inside browse tree")));
    await session.step(50, "And first Databases tree node inside browse tree should be collapsed", () => shouldBe(page, el("first Databases tree node inside browse tree"), "collapsed"));
    await session.step(51, "When user expands first Databases tree node inside browse tree", () => expand(page, el("first Databases tree node inside browse tree")));
    await session.step(52, "Then first Databases tree node inside browse tree should be expanded", () => shouldBe(page, el("first Databases tree node inside browse tree"), "expanded"));
    await session.step(53, "And Databases---Postgres tree node inside browse tree should be visible", () => shouldBe(page, el("Databases---Postgres tree node inside browse tree"), "visible"));
    await session.step(54, "And no errors should have been logged", () => noErrors(page));
    await session.step(55, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Opening the Platform section logs no error [section=Platform, child=Platform---Users]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(22, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(49, "Given user collapses first Platform tree node inside browse tree", () => collapse(page, el("first Platform tree node inside browse tree")));
    await session.step(50, "And first Platform tree node inside browse tree should be collapsed", () => shouldBe(page, el("first Platform tree node inside browse tree"), "collapsed"));
    await session.step(51, "When user expands first Platform tree node inside browse tree", () => expand(page, el("first Platform tree node inside browse tree")));
    await session.step(52, "Then first Platform tree node inside browse tree should be expanded", () => shouldBe(page, el("first Platform tree node inside browse tree"), "expanded"));
    await session.step(53, "And Platform---Users tree node inside browse tree should be visible", () => shouldBe(page, el("Platform---Users tree node inside browse tree"), "visible"));
    await session.step(54, "And no errors should have been logged", () => noErrors(page));
    await session.step(55, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
