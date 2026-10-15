/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/search/power-search-results.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [powerpack.search.power-pack, powerpack.view.welcome]
--- */
import {test} from '@playwright/test';
import '../../bindings/add-new-column.js';
import '../../bindings/enrichment.js';
import '../../bindings/home.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {categoryLists, noResults, searchFinished} from '../../bindings/search.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldHaveValue, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, scriptOnServer, taggedProjectOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("What the Home search finds", () => {
  const session = feature(test, "features/search/power-search-results.feature", import.meta.url);
  test("The Home icon of the Browse toolbar brings back the Home page with its search box", {tag: ["@realizes:powerpack.search.power-pack", "@realizes:powerpack.view.welcome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "Given the browse panel is open", () => browsePanelOpen(page));
    await session.step(19, "And user clicks on \"Open text\" icon inside browse toolbar", () => clickOn(page, el("\"Open text\" icon inside browse toolbar")));
    await session.step(20, "And the \"Import text\" view should be current", () => viewIsCurrent(page, "Import text"));
    await session.step(21, "And home search should be hidden", () => shouldBe(page, el("home search"), "hidden"));
    await session.step(22, "When user clicks on \"Home\" icon inside browse toolbar", () => clickOn(page, el("\"Home\" icon inside browse toolbar")));
    await session.step(23, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
    await session.step(24, "And home search should be visible", () => shouldBe(page, el("home search"), "visible"));
    await session.step(25, "And no errors should have been logged", () => noErrors(page));
  });
  test("A script is found by its name", {tag: ["@realizes:powerpack.search.power-pack", "@realizes:powerpack.view.welcome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "Given a script \"BddPowerSearch{time}\" is on the server:", () => scriptOnServer(page, session.text("BddPowerSearch{time}"), "//language: javascript\n//output: int x\nx = 1;"));
    await session.step(34, "When user types \"BddPowerSearch{time}\" into home search", () => typeInto(page, session.text("BddPowerSearch{time}"), el("home search")));
    await session.step(35, "Then the search should have finished", () => searchFinished(page));
    await session.step(36, "And the \"Scripts\" search results should list \"BddPowerSearch{time}\"", () => categoryLists(page, "Scripts", session.text("BddPowerSearch{time}")));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
    await session.step(38, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A tag search finds the projects that carry the tag", {tag: ["@realizes:powerpack.search.power-pack", "@realizes:powerpack.view.welcome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "Given a project named \"bdd-power-search-{time}\" tagged \"bddsearch{time}\" is on the server", () => taggedProjectOnServer(page, session.text("bdd-power-search-{time}"), session.text("bddsearch{time}")));
    await session.step(42, "When user types \"#bddsearch{time}\" into home search", () => typeInto(page, session.text("#bddsearch{time}"), el("home search")));
    await session.step(43, "Then the search should have finished", () => searchFinished(page));
    await session.step(44, "And the \"Projects\" search results should list \"bdd-power-search-{time}\"", () => categoryLists(page, "Projects", session.text("bdd-power-search-{time}")));
    await session.step(45, "And no errors should have been logged", () => noErrors(page));
    await session.step(46, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A text nothing matches finds nothing", {tag: ["@realizes:powerpack.search.power-pack", "@realizes:powerpack.view.welcome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(49, "When user types \"zzqxnomatchqq\" into home search", () => typeInto(page, "zzqxnomatchqq", el("home search")));
    await session.step(50, "Then the search should have finished", () => searchFinished(page));
    await session.step(51, "And the search results should be empty", () => noResults(page));
    await session.step(52, "And home search should have value \"zzqxnomatchqq\"", () => shouldHaveValue(page, el("home search"), "zzqxnomatchqq"));
    await session.step(53, "And no errors should have been logged", () => noErrors(page));
    await session.step(54, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
