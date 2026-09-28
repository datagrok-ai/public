/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/browse-gallery-filters.feature
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
import {clickOn, isExpanded, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, galleryCountNotLower, rememberGalleryCount, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The quick filters of a Browse gallery", () => {
  const session = feature(test, "features/browse/browse-gallery-filters.feature", import.meta.url);
  test("Created recently puts its query in the search and All takes it back", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(23, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(24, "When user clicks on Files---Demo tree node inside browse tree", () => clickOn(page, el("Files---Demo tree node inside browse tree")));
    await session.step(25, "Then the \"Demo\" view should be current", () => viewIsCurrent(page, "Demo"));
    await session.step(26, "When user clicks on \"Toggle filters\" icon", () => clickOn(page, el("\"Toggle filters\" icon")));
    await session.step(27, "And user remembers the gallery counter", () => rememberGalleryCount(page));
    await session.step(28, "And user clicks on \"Created recently\" tag", () => clickOn(page, el("\"Created recently\" tag")));
    await session.step(29, "Then gallery search should have value \"created < 1w\"", () => shouldHaveValue(page, el("gallery search"), "created < 1w"));
    await session.step(30, "When user clicks on \"All\" tag", () => clickOn(page, el("\"All\" tag")));
    await session.step(31, "Then gallery search should have value \"\"", () => shouldHaveValue(page, el("gallery search"), ""));
    await session.step(32, "And the gallery counter should not be lower than remembered", () => galleryCountNotLower(page));
    await session.step(33, "And no errors should have been logged", () => noErrors(page));
    await session.step(34, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Dockers gallery offers Used by me and Created recently", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(37, "Given Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(38, "When user clicks on Platform---Dockers tree node inside browse tree", () => clickOn(page, el("Platform---Dockers tree node inside browse tree")));
    await session.step(39, "Then the \"Dockers\" view should be current", () => viewIsCurrent(page, "Dockers"));
    await session.step(40, "When user clicks on \"Toggle filters\" icon", () => clickOn(page, el("\"Toggle filters\" icon")));
    await session.step(41, "Then \"Used by me\" tag should be visible", () => shouldBe(page, el("\"Used by me\" tag"), "visible"));
    await session.step(42, "And \"Created recently\" tag should be visible", () => shouldBe(page, el("\"Created recently\" tag"), "visible"));
    await session.step(43, "When user remembers the gallery counter", () => rememberGalleryCount(page));
    await session.step(44, "And user clicks on \"Used by me\" tag", () => clickOn(page, el("\"Used by me\" tag")));
    await session.step(45, "Then gallery search should have value \"usedBy = @current\"", () => shouldHaveValue(page, el("gallery search"), "usedBy = @current"));
    await session.step(46, "When user clicks on \"All\" tag", () => clickOn(page, el("\"All\" tag")));
    await session.step(47, "Then gallery search should have value \"\"", () => shouldHaveValue(page, el("gallery search"), ""));
    await session.step(48, "And the gallery counter should not be lower than remembered", () => galleryCountNotLower(page));
    await session.step(49, "And no errors should have been logged", () => noErrors(page));
    await session.step(50, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
