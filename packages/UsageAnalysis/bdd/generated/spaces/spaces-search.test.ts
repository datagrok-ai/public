/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/spaces/spaces-search.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.space]
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
import {clearField, clickOn, doubleClickOn, enterInto, finishedUpdating, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, childSpaceOnServer, spaceOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Searching spaces", () => {
  const session = feature(test, "features/spaces/spaces-search.feature", import.meta.url);
  test("Searching spaces", {tag: ["@journey", "@spaces", "@realizes:views.space", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(18, "And a space named \"BDD-Find\" is on the server", () => spaceOnServer(page, "BDD-Find"));
    await session.step(19, "And a space named \"BDD-Miss\" is on the server", () => spaceOnServer(page, "BDD-Miss"));
    await session.step(20, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await run.scenario("The Spaces list shows both", async () => {
      await session.step(23, "When user clicks on Spaces tree node inside browse tree", () => clickOn(page, el("Spaces tree node inside browse tree")));
      await session.step(24, "And user clicks on \"Refresh\" icon", () => clickOn(page, el("\"Refresh\" icon")));
      await session.step(25, "Then BDD-Find link in gallery should be visible", () => shouldBe(page, el("BDD-Find link in gallery"), "visible"));
      await session.step(26, "And BDD-Miss link in gallery should be visible", () => shouldBe(page, el("BDD-Miss link in gallery"), "visible"));
    });
    await run.scenario("A whole name keeps only that space", async () => {
      await session.step(29, "When user enters \"BDD-Find\" into gallery search", () => enterInto(page, "BDD-Find", el("gallery search")));
      await session.step(30, "Then BDD-Find link in gallery should be visible", () => shouldBe(page, el("BDD-Find link in gallery"), "visible"));
      await session.step(31, "And BDD-Miss link in gallery should be absent", () => shouldBe(page, el("BDD-Miss link in gallery"), "absent"));
    });
    await run.scenario("Part of a name still matches", async () => {
      await session.step(34, "When user enters \"BDD-Fi\" into gallery search", () => enterInto(page, "BDD-Fi", el("gallery search")));
      await session.step(35, "Then BDD-Find link in gallery should be visible", () => shouldBe(page, el("BDD-Find link in gallery"), "visible"));
      await session.step(36, "And BDD-Miss link in gallery should be absent", () => shouldBe(page, el("BDD-Miss link in gallery"), "absent"));
    });
    await run.scenario("A name nothing carries empties the list", async () => {
      await session.step(39, "When user enters \"zzz-no-such-space\" into gallery search", () => enterInto(page, "zzz-no-such-space", el("gallery search")));
      await session.step(40, "Then BDD-Find link in gallery should be absent", () => shouldBe(page, el("BDD-Find link in gallery"), "absent"));
      await session.step(41, "And BDD-Miss link in gallery should be absent", () => shouldBe(page, el("BDD-Miss link in gallery"), "absent"));
    });
    await run.scenario("Clearing the search brings both back", async () => {
      await session.step(44, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(45, "Then BDD-Find link in gallery should be visible", () => shouldBe(page, el("BDD-Find link in gallery"), "visible"));
      await session.step(46, "And BDD-Miss link in gallery should be visible", () => shouldBe(page, el("BDD-Miss link in gallery"), "visible"));
    });
    await run.scenario("A child space is searchable inside its parent", async () => {
      await session.step(49, "Given a space named \"BDD-Find-Child\" under \"BDD-Find\" is on the server", () => childSpaceOnServer(page, "BDD-Find-Child", "BDD-Find"));
      await session.step(50, "When user double-clicks on Spaces---BDD-Find tree node inside browse tree", () => doubleClickOn(page, el("Spaces---BDD-Find tree node inside browse tree")));
      await session.step(51, "Then BDD-Find-Child link in gallery should be visible", () => shouldBe(page, el("BDD-Find-Child link in gallery"), "visible"));
      await session.step(52, "When user enters \"zzz-no-such-space\" into gallery search", () => enterInto(page, "zzz-no-such-space", el("gallery search")));
      await session.step(53, "Then BDD-Find-Child link in gallery should be absent", () => shouldBe(page, el("BDD-Find-Child link in gallery"), "absent"));
      await session.step(54, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(55, "Then BDD-Find-Child link in gallery should be visible", () => shouldBe(page, el("BDD-Find-Child link in gallery"), "visible"));
    });
    await run.scenario("Characters of the search syntax match no space, and clearing the search brings both back", async () => {
      await session.step(59, "When user clicks on Spaces tree node inside browse tree", () => clickOn(page, el("Spaces tree node inside browse tree")));
      await session.step(60, "And user clicks on \"Refresh\" icon", () => clickOn(page, el("\"Refresh\" icon")));
      await session.step(61, "Then BDD-Find link in gallery should be visible", () => shouldBe(page, el("BDD-Find link in gallery"), "visible"));
      await session.step(62, "And BDD-Miss link in gallery should be visible", () => shouldBe(page, el("BDD-Miss link in gallery"), "visible"));
      await session.step(63, "When user enters \"(\" into gallery search", () => enterInto(page, "(", el("gallery search")));
      await session.step(64, "Then gallery should have finished updating", () => finishedUpdating(page, el("gallery")));
      await session.step(65, "And BDD-Find link in gallery should be absent", () => shouldBe(page, el("BDD-Find link in gallery"), "absent"));
      await session.step(66, "And BDD-Miss link in gallery should be absent", () => shouldBe(page, el("BDD-Miss link in gallery"), "absent"));
      await session.step(67, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(68, "Then BDD-Find link in gallery should be visible", () => shouldBe(page, el("BDD-Find link in gallery"), "visible"));
      await session.step(69, "When user enters \"\\\"\" into gallery search", () => enterInto(page, "\"", el("gallery search")));
      await session.step(70, "Then gallery should have finished updating", () => finishedUpdating(page, el("gallery")));
      await session.step(71, "And BDD-Find link in gallery should be absent", () => shouldBe(page, el("BDD-Find link in gallery"), "absent"));
      await session.step(72, "And BDD-Miss link in gallery should be absent", () => shouldBe(page, el("BDD-Miss link in gallery"), "absent"));
      await session.step(73, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(74, "Then BDD-Find link in gallery should be visible", () => shouldBe(page, el("BDD-Find link in gallery"), "visible"));
      await session.step(75, "When user enters \"name =\" into gallery search", () => enterInto(page, "name =", el("gallery search")));
      await session.step(76, "Then gallery should have finished updating", () => finishedUpdating(page, el("gallery")));
      await session.step(77, "And BDD-Find link in gallery should be absent", () => shouldBe(page, el("BDD-Find link in gallery"), "absent"));
      await session.step(78, "And BDD-Miss link in gallery should be absent", () => shouldBe(page, el("BDD-Miss link in gallery"), "absent"));
      await session.step(79, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(80, "Then BDD-Find link in gallery should be visible", () => shouldBe(page, el("BDD-Find link in gallery"), "visible"));
      await session.step(81, "And BDD-Miss link in gallery should be visible", () => shouldBe(page, el("BDD-Miss link in gallery"), "visible"));
      await session.step(82, "And no errors should have been logged", () => noErrors(page));
      await session.step(83, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
