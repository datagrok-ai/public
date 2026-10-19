/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/browse-gallery-filters.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.browse]
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
import {clearField, clickOn, enterInto, isExpanded, shouldBe, shouldHaveText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {inFavorites} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {browsePanelOpen, galleryCountHigher, galleryCountLower, galleryCountNotLower, rememberGalleryCount, scriptOnServer, scriptsView, standRunsService, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The quick filters of a Browse gallery", () => {
  const session = feature(test, "features/browse/browse-gallery-filters.feature", import.meta.url);
  test("Created recently puts its query in the search and All takes it back", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(30, "Given user is logged in", () => loggedIn(page));
    await session.step(31, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(34, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(35, "When user clicks on Files---Demo tree node inside browse tree", () => clickOn(page, el("Files---Demo tree node inside browse tree")));
    await session.step(36, "Then the \"Demo\" view should be current", () => viewIsCurrent(page, "Demo"));
    await session.step(37, "When user clicks on \"Toggle filters\" icon", () => clickOn(page, el("\"Toggle filters\" icon")));
    await session.step(38, "And user remembers the gallery counter", () => rememberGalleryCount(page));
    await session.step(39, "And user clicks on \"Created recently\" tag", () => clickOn(page, el("\"Created recently\" tag")));
    await session.step(40, "Then gallery search should have value \"created < 1w\"", () => shouldHaveValue(page, el("gallery search"), "created < 1w"));
    await session.step(41, "When user clicks on \"All\" tag", () => clickOn(page, el("\"All\" tag")));
    await session.step(42, "Then gallery search should have value \"\"", () => shouldHaveValue(page, el("gallery search"), ""));
    await session.step(43, "And the gallery counter should not be lower than remembered", () => galleryCountNotLower(page));
    await session.step(44, "And no errors should have been logged", () => noErrors(page));
    await session.step(45, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Dockers gallery offers Used by me and Created recently", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(30, "Given user is logged in", () => loggedIn(page));
    await session.step(31, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(49, "Given the stand runs the \"Grok Spawner\" service", () => standRunsService(page, "Grok Spawner"));
    await session.step(50, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(51, "When user clicks on Platform---Dockers tree node inside browse tree", () => clickOn(page, el("Platform---Dockers tree node inside browse tree")));
    await session.step(52, "Then the \"Dockers\" view should be current", () => viewIsCurrent(page, "Dockers"));
    await session.step(53, "When user clicks on \"Toggle filters\" icon", () => clickOn(page, el("\"Toggle filters\" icon")));
    await session.step(54, "Then \"Used by me\" tag should be visible", () => shouldBe(page, el("\"Used by me\" tag"), "visible"));
    await session.step(55, "And \"Created recently\" tag should be visible", () => shouldBe(page, el("\"Created recently\" tag"), "visible"));
    await session.step(56, "When user remembers the gallery counter", () => rememberGalleryCount(page));
    await session.step(57, "And user clicks on \"Used by me\" tag", () => clickOn(page, el("\"Used by me\" tag")));
    await session.step(58, "Then gallery search should have value \"usedBy = @current\"", () => shouldHaveValue(page, el("gallery search"), "usedBy = @current"));
    await session.step(59, "When user clicks on \"All\" tag", () => clickOn(page, el("\"All\" tag")));
    await session.step(60, "Then gallery search should have value \"\"", () => shouldHaveValue(page, el("gallery search"), ""));
    await session.step(61, "And the gallery counter should not be lower than remembered", () => galleryCountNotLower(page));
    await session.step(62, "And no errors should have been logged", () => noErrors(page));
    await session.step(63, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A quick filter picked over a search keeps only what both let through", {tag: ["@browse", "@realizes:views.browse", "@serial", "@GROK-19690"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(30, "Given user is logged in", () => loggedIn(page));
    await session.step(31, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(67, "Given a script \"BddGallery{time}Starred\" is on the server:", () => scriptOnServer(page, session.text("BddGallery{time}Starred"), "//language: javascript\n//output: int x\nx = 1;"));
    await session.step(73, "And a script \"BddGallery{time}Other\" is on the server:", () => scriptOnServer(page, session.text("BddGallery{time}Other"), "//language: javascript\n//output: int x\nx = 2;"));
    await session.step(79, "And the script \"BddGallery{time}Starred\" is in favorites", () => inFavorites(page, "script", session.text("BddGallery{time}Starred")));
    await session.step(80, "And user opens the Scripts view", () => scriptsView(page));
    await session.step(81, "When user enters \"BddGallery{time}\" into gallery search", () => enterInto(page, session.text("BddGallery{time}"), el("gallery search")));
    await session.step(82, "Then gallery counter should have text \"2\"", () => shouldHaveText(page, el("gallery counter"), "2"));
    await session.step(83, "When user remembers the gallery counter", () => rememberGalleryCount(page));
    await session.step(84, "And user clicks on \"Toggle filters\" icon", () => clickOn(page, el("\"Toggle filters\" icon")));
    await session.step(85, "And user clicks on \"Favorites\" tag", () => clickOn(page, el("\"Favorites\" tag")));
    await session.step(86, "Then gallery counter should have text \"1\"", () => shouldHaveText(page, el("gallery counter"), "1"));
    await session.step(87, "And \"BddGallery{time}Starred\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BddGallery{time}Starred\" link in gallery")), "visible"));
    await session.step(88, "And \"BddGallery{time}Other\" link in gallery should be absent", () => shouldBe(page, el(session.text("\"BddGallery{time}Other\" link in gallery")), "absent"));
    await session.step(89, "And gallery search should have value \"BddGallery{time}\"", () => shouldHaveValue(page, el("gallery search"), session.text("BddGallery{time}")));
    await session.step(90, "When user clicks on \"All\" tag", () => clickOn(page, el("\"All\" tag")));
    await session.step(91, "Then gallery search should have value \"\"", () => shouldHaveValue(page, el("gallery search"), ""));
    await session.step(92, "And the gallery counter should be higher than remembered", () => galleryCountHigher(page));
    await session.step(93, "When user enters \"BddGallery{time}\" into gallery search", () => enterInto(page, session.text("BddGallery{time}"), el("gallery search")));
    await session.step(94, "Then gallery counter should have text \"2\"", () => shouldHaveText(page, el("gallery counter"), "2"));
    await session.step(95, "When user clears gallery search", () => clearField(page, el("gallery search")));
    await session.step(96, "Then gallery search should have value \"\"", () => shouldHaveValue(page, el("gallery search"), ""));
    await session.step(97, "And no errors should have been logged", () => noErrors(page));
    await session.step(98, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A search typed over a quick filter narrows what the filter lets through", {tag: ["@browse", "@realizes:views.browse", "@serial", "@GROK-19690"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(30, "Given user is logged in", () => loggedIn(page));
    await session.step(31, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(102, "Given a script \"BddGallery{time}Starred\" is on the server:", () => scriptOnServer(page, session.text("BddGallery{time}Starred"), "//language: javascript\n//output: int x\nx = 1;"));
    await session.step(108, "And a script \"BddGallery{time}Other\" is on the server:", () => scriptOnServer(page, session.text("BddGallery{time}Other"), "//language: javascript\n//output: int x\nx = 2;"));
    await session.step(114, "And a script \"BddFavoriteKeptAside\" is on the server:", () => scriptOnServer(page, "BddFavoriteKeptAside", "//language: javascript\n//output: int x\nx = 3;"));
    await session.step(120, "And the script \"BddGallery{time}Starred\" is in favorites", () => inFavorites(page, "script", session.text("BddGallery{time}Starred")));
    await session.step(121, "And the script \"BddFavoriteKeptAside\" is in favorites", () => inFavorites(page, "script", "BddFavoriteKeptAside"));
    await session.step(122, "And user opens the Scripts view", () => scriptsView(page));
    await session.step(123, "When user clicks on \"Toggle filters\" icon", () => clickOn(page, el("\"Toggle filters\" icon")));
    await session.step(124, "And user remembers the gallery counter", () => rememberGalleryCount(page));
    await session.step(125, "And user clicks on \"Favorites\" tag", () => clickOn(page, el("\"Favorites\" tag")));
    await session.step(126, "Then gallery search should have value \"starredBy = @current\"", () => shouldHaveValue(page, el("gallery search"), "starredBy = @current"));
    await session.step(127, "And the gallery counter should be lower than remembered", () => galleryCountLower(page));
    await session.step(128, "And \"BddGallery{time}Starred\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BddGallery{time}Starred\" link in gallery")), "visible"));
    await session.step(129, "And \"BddFavoriteKeptAside\" link in gallery should be visible", () => shouldBe(page, el("\"BddFavoriteKeptAside\" link in gallery"), "visible"));
    await session.step(130, "And \"BddGallery{time}Other\" link in gallery should be absent", () => shouldBe(page, el(session.text("\"BddGallery{time}Other\" link in gallery")), "absent"));
    await session.step(131, "When user enters \"BddGallery{time}\" into gallery search", () => enterInto(page, session.text("BddGallery{time}"), el("gallery search")));
    await session.step(132, "Then gallery counter should have text \"1\"", () => shouldHaveText(page, el("gallery counter"), "1"));
    await session.step(133, "And \"BddGallery{time}Starred\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BddGallery{time}Starred\" link in gallery")), "visible"));
    await session.step(134, "And \"BddFavoriteKeptAside\" link in gallery should be absent", () => shouldBe(page, el("\"BddFavoriteKeptAside\" link in gallery"), "absent"));
    await session.step(135, "And \"BddGallery{time}Other\" link in gallery should be absent", () => shouldBe(page, el(session.text("\"BddGallery{time}Other\" link in gallery")), "absent"));
    await session.step(136, "When user clicks on \"All\" tag", () => clickOn(page, el("\"All\" tag")));
    await session.step(137, "Then gallery search should have value \"\"", () => shouldHaveValue(page, el("gallery search"), ""));
    await session.step(138, "And the gallery counter should not be lower than remembered", () => galleryCountNotLower(page));
    await session.step(139, "And no errors should have been logged", () => noErrors(page));
    await session.step(140, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
