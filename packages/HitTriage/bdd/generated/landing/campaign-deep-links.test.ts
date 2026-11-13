/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/landing/campaign-deep-links.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {campaignsOnServer, listed, notListed, settingsDefaults, sortedBy} from '../../bindings/landing.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {openAddress} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, currentViewShows, packageInstalled, urlShouldContain, urlShouldNotContain, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {errorBalloonText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Campaign links: ?campaignId opens a campaign, and a link the user cannot follow lands on the landing", () => {
  const session = feature(test, "features/landing/campaign-deep-links.feature", import.meta.url);
  test("Campaign links: ?campaignId opens a campaign, and a link the user cannot follow lands on the landing", {tag: ["@journey", "@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the campaigns tables start from their defaults, put back at feature end", () => settingsDefaults(page));
    await session.step(18, "And the \"Hit Design\" campaigns are on the server:", () => campaignsOnServer(page, "Hit Design", [["code","name","created","author","can view","can edit","template"],[session.text("BDDDLK{time}-1"),session.text("BDD Link Open {time}"),"2099/1/13","me","everyone","everyone",""],[session.text("BDDDLK{time}-2"),session.text("BDD Link Hidden {time}"),"2099/1/12","sharing user","author","author",""],[session.text("BDDDLK{time}-3"),session.text("BDD Link Broken {time}"),"2099/1/11","me","everyone","everyone",session.text("gone: BDD Link Gone Template {time}")]]), [["code","name","created","author","can view","can edit","template"],[session.text("BDDDLK{time}-1"),session.text("BDD Link Open {time}"),"2099/1/13","me","everyone","everyone",""],[session.text("BDDDLK{time}-2"),session.text("BDD Link Hidden {time}"),"2099/1/12","sharing user","author","author",""],[session.text("BDDDLK{time}-3"),session.text("BDD Link Broken {time}"),"2099/1/11","me","everyone","everyone",session.text("gone: BDD Link Gone Template {time}")]]);
    await session.step(23, "And the \"PeptiHit\" campaigns are on the server:", () => campaignsOnServer(page, "PeptiHit", [["code","name","created","author","can view","can edit"],[session.text("BDDDLP{time}-1"),session.text("BDD Link PH Hidden {time}"),"2099/1/12","sharing user","author","author"]]), [["code","name","created","author","can view","can edit"],[session.text("BDDDLP{time}-1"),session.text("BDD Link PH Hidden {time}"),"2099/1/12","sharing user","author","author"]]);
    await session.step(26, "And the campaigns table of \"Hit Design\" is sorted by \"Created\" descending", () => sortedBy(page, "Hit Design", "Created", "descending"));
    await run.scenario("A link with a campaign's code opens the campaign", async () => {
      await session.step(29, "When user opens the address \"/apps/HitTriage/HitDesign?campaignId=BDDDLK{time}-1\"", () => openAddress(page, session.text("/apps/HitTriage/HitDesign?campaignId=BDDDLK{time}-1")));
      await session.step(30, "Then the \"BDD Link Open {time}\" view should be current", () => viewIsCurrent(page, session.text("BDD Link Open {time}")));
      await session.step(31, "And the page address should contain \"campaignId=BDDDLK{time}-1\"", () => urlShouldContain(page, session.text("campaignId=BDDDLK{time}-1")));
      await session.step(32, "And the table should have 1 row", () => rowCount(page, 1));
    });
    await run.scenario("A link to a campaign the user may not view lands on the landing and drops the code", async () => {
      await session.step(35, "When user closes all views", () => closeAllViews(page));
      await session.step(36, "And user opens the address \"/apps/HitTriage/HitDesign?campaignId=BDDDLK{time}-2\"", () => openAddress(page, session.text("/apps/HitTriage/HitDesign?campaignId=BDDDLK{time}-2")));
      await session.step(37, "Then the \"Hit Design\" view should be current", () => viewIsCurrent(page, "Hit Design"));
      await session.step(38, "And the campaigns table should list \"BDD Link Open {time}\"", () => listed(page, session.text("BDD Link Open {time}")));
      await session.step(39, "And the campaigns table should not list \"BDD Link Hidden {time}\"", () => notListed(page, session.text("BDD Link Hidden {time}")));
      await session.step(40, "And the page address should not contain \"campaignId\"", () => urlShouldNotContain(page, "campaignId"));
    });
    await run.scenario("A link with a code no campaign has lands on the landing", async () => {
      await session.step(43, "When user closes all views", () => closeAllViews(page));
      await session.step(44, "And user opens the address \"/apps/HitTriage/HitDesign?campaignId=BDDDLK{time}-99\"", () => openAddress(page, session.text("/apps/HitTriage/HitDesign?campaignId=BDDDLK{time}-99")));
      await session.step(45, "Then the \"Hit Design\" view should be current", () => viewIsCurrent(page, "Hit Design"));
      await session.step(46, "And the campaigns table should list \"BDD Link Open {time}\"", () => listed(page, session.text("BDD Link Open {time}")));
    });
    await run.scenario("A link to a campaign that cannot be opened still lands on a working landing", async () => {
      await session.step(49, "When user closes all views", () => closeAllViews(page));
      await session.step(50, "And user opens the address \"/apps/HitTriage/HitDesign?campaignId=BDDDLK{time}-3\"", () => openAddress(page, session.text("/apps/HitTriage/HitDesign?campaignId=BDDDLK{time}-3")));
      await session.step(51, "Then the \"Hit Design\" view should be current", () => viewIsCurrent(page, "Hit Design"));
      await session.step(52, "And the campaigns table should list \"BDD Link Open {time}\"", () => listed(page, session.text("BDD Link Open {time}")));
      await session.step(53, "And an error balloon containing \"Campaign BDDDLK{time}-3 could not be opened\" should have been shown", () => errorBalloonText(page, session.text("Campaign BDDDLK{time}-3 could not be opened")));
    });
    await run.scenario("In PeptiHit, a link to a campaign the user may not view lands on its landing and drops the code", async () => {
      await session.step(56, "Given the \"Bio\" package is installed", () => packageInstalled(page, "Bio"));
      await session.step(57, "When user closes all views", () => closeAllViews(page));
      await session.step(58, "And user opens the address \"/apps/HitTriage/PeptiHit?campaignId=BDDDLP{time}-1\"", () => openAddress(page, session.text("/apps/HitTriage/PeptiHit?campaignId=BDDDLP{time}-1")));
      await session.step(59, "Then the \"PeptiHit\" view should be current", () => viewIsCurrent(page, "PeptiHit"));
      await session.step(60, "And the current view should show the text \"Continue Campaign\"", () => currentViewShows(page, "Continue Campaign"));
      await session.step(61, "And the page address should not contain \"campaignId\"", () => urlShouldNotContain(page, "campaignId"));
    });
    run.finish();
  });
});
