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
    await session.step(21, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(24, "When user clicks on My stuff tree node inside browse tree", () => clickOn(page, el("My stuff tree node inside browse tree")));
    await session.step(25, "Then My stuff tree node inside browse tree should be selected", () => shouldBe(page, el("My stuff tree node inside browse tree"), "selected"));
    await session.step(26, "And no errors should have been logged", () => noErrors(page));
    await session.step(27, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Clicking Spaces logs no error [node=Spaces]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(24, "When user clicks on Spaces tree node inside browse tree", () => clickOn(page, el("Spaces tree node inside browse tree")));
    await session.step(25, "Then Spaces tree node inside browse tree should be selected", () => shouldBe(page, el("Spaces tree node inside browse tree"), "selected"));
    await session.step(26, "And no errors should have been logged", () => noErrors(page));
    await session.step(27, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Clicking Apps logs no error [node=Apps]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(24, "When user clicks on Apps tree node inside browse tree", () => clickOn(page, el("Apps tree node inside browse tree")));
    await session.step(25, "Then Apps tree node inside browse tree should be selected", () => shouldBe(page, el("Apps tree node inside browse tree"), "selected"));
    await session.step(26, "And no errors should have been logged", () => noErrors(page));
    await session.step(27, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Clicking Files logs no error [node=Files]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(24, "When user clicks on Files tree node inside browse tree", () => clickOn(page, el("Files tree node inside browse tree")));
    await session.step(25, "Then Files tree node inside browse tree should be selected", () => shouldBe(page, el("Files tree node inside browse tree"), "selected"));
    await session.step(26, "And no errors should have been logged", () => noErrors(page));
    await session.step(27, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Clicking Dashboards logs no error [node=Dashboards]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(24, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(25, "Then Dashboards tree node inside browse tree should be selected", () => shouldBe(page, el("Dashboards tree node inside browse tree"), "selected"));
    await session.step(26, "And no errors should have been logged", () => noErrors(page));
    await session.step(27, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Clicking Databases logs no error [node=Databases]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(24, "When user clicks on Databases tree node inside browse tree", () => clickOn(page, el("Databases tree node inside browse tree")));
    await session.step(25, "Then Databases tree node inside browse tree should be selected", () => shouldBe(page, el("Databases tree node inside browse tree"), "selected"));
    await session.step(26, "And no errors should have been logged", () => noErrors(page));
    await session.step(27, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Clicking Platform logs no error [node=Platform]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(24, "When user clicks on Platform tree node inside browse tree", () => clickOn(page, el("Platform tree node inside browse tree")));
    await session.step(25, "Then Platform tree node inside browse tree should be selected", () => shouldBe(page, el("Platform tree node inside browse tree"), "selected"));
    await session.step(26, "And no errors should have been logged", () => noErrors(page));
    await session.step(27, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Opening the My stuff section logs no error [section=My stuff]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(42, "Given user collapses My stuff tree node inside browse tree", () => collapse(page, el("My stuff tree node inside browse tree")));
    await session.step(43, "And My stuff tree node inside browse tree should be collapsed", () => shouldBe(page, el("My stuff tree node inside browse tree"), "collapsed"));
    await session.step(44, "When user expands My stuff tree node inside browse tree", () => expand(page, el("My stuff tree node inside browse tree")));
    await session.step(45, "Then My stuff tree node inside browse tree should be expanded", () => shouldBe(page, el("My stuff tree node inside browse tree"), "expanded"));
    await session.step(46, "And no errors should have been logged", () => noErrors(page));
    await session.step(47, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Opening the Apps section logs no error [section=Apps]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(42, "Given user collapses Apps tree node inside browse tree", () => collapse(page, el("Apps tree node inside browse tree")));
    await session.step(43, "And Apps tree node inside browse tree should be collapsed", () => shouldBe(page, el("Apps tree node inside browse tree"), "collapsed"));
    await session.step(44, "When user expands Apps tree node inside browse tree", () => expand(page, el("Apps tree node inside browse tree")));
    await session.step(45, "Then Apps tree node inside browse tree should be expanded", () => shouldBe(page, el("Apps tree node inside browse tree"), "expanded"));
    await session.step(46, "And no errors should have been logged", () => noErrors(page));
    await session.step(47, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Opening the Files section logs no error [section=Files]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(42, "Given user collapses Files tree node inside browse tree", () => collapse(page, el("Files tree node inside browse tree")));
    await session.step(43, "And Files tree node inside browse tree should be collapsed", () => shouldBe(page, el("Files tree node inside browse tree"), "collapsed"));
    await session.step(44, "When user expands Files tree node inside browse tree", () => expand(page, el("Files tree node inside browse tree")));
    await session.step(45, "Then Files tree node inside browse tree should be expanded", () => shouldBe(page, el("Files tree node inside browse tree"), "expanded"));
    await session.step(46, "And no errors should have been logged", () => noErrors(page));
    await session.step(47, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Opening the Databases section logs no error [section=Databases]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(42, "Given user collapses Databases tree node inside browse tree", () => collapse(page, el("Databases tree node inside browse tree")));
    await session.step(43, "And Databases tree node inside browse tree should be collapsed", () => shouldBe(page, el("Databases tree node inside browse tree"), "collapsed"));
    await session.step(44, "When user expands Databases tree node inside browse tree", () => expand(page, el("Databases tree node inside browse tree")));
    await session.step(45, "Then Databases tree node inside browse tree should be expanded", () => shouldBe(page, el("Databases tree node inside browse tree"), "expanded"));
    await session.step(46, "And no errors should have been logged", () => noErrors(page));
    await session.step(47, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Opening the Platform section logs no error [section=Platform]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(42, "Given user collapses Platform tree node inside browse tree", () => collapse(page, el("Platform tree node inside browse tree")));
    await session.step(43, "And Platform tree node inside browse tree should be collapsed", () => shouldBe(page, el("Platform tree node inside browse tree"), "collapsed"));
    await session.step(44, "When user expands Platform tree node inside browse tree", () => expand(page, el("Platform tree node inside browse tree")));
    await session.step(45, "Then Platform tree node inside browse tree should be expanded", () => shouldBe(page, el("Platform tree node inside browse tree"), "expanded"));
    await session.step(46, "And no errors should have been logged", () => noErrors(page));
    await session.step(47, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
