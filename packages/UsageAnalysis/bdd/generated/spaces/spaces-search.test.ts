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
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, doubleClickOn, enterInto, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, childSpaceOnServer, spaceOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Searching spaces", () => {
  const session = feature(test, "features/spaces/spaces-search.feature", import.meta.url);
  test("Searching spaces", {tag: ["@journey", "@spaces", "@realizes:views.space"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(14, "And a space named \"BDD-Find\" is on the server", () => spaceOnServer(page, "BDD-Find"));
    await session.step(15, "And a space named \"BDD-Miss\" is on the server", () => spaceOnServer(page, "BDD-Miss"));
    await session.step(16, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await run.scenario("The Spaces list shows both", async () => {
      await session.step(19, "When user clicks on Spaces tree node inside browse tree", () => clickOn(page, el("Spaces tree node inside browse tree")));
      await session.step(20, "And user clicks on \"Refresh\" icon", () => clickOn(page, el("\"Refresh\" icon")));
      await session.step(21, "Then BDD-Find link in gallery should be visible", () => shouldBe(page, el("BDD-Find link in gallery"), "visible"));
      await session.step(22, "And BDD-Miss link in gallery should be visible", () => shouldBe(page, el("BDD-Miss link in gallery"), "visible"));
    });
    await run.scenario("A whole name keeps only that space", async () => {
      await session.step(25, "When user enters \"BDD-Find\" into gallery search", () => enterInto(page, "BDD-Find", el("gallery search")));
      await session.step(26, "Then BDD-Find link in gallery should be visible", () => shouldBe(page, el("BDD-Find link in gallery"), "visible"));
      await session.step(27, "And BDD-Miss link in gallery should be absent", () => shouldBe(page, el("BDD-Miss link in gallery"), "absent"));
    });
    await run.scenario("Part of a name still matches", async () => {
      await session.step(30, "When user enters \"BDD-Fi\" into gallery search", () => enterInto(page, "BDD-Fi", el("gallery search")));
      await session.step(31, "Then BDD-Find link in gallery should be visible", () => shouldBe(page, el("BDD-Find link in gallery"), "visible"));
      await session.step(32, "And BDD-Miss link in gallery should be absent", () => shouldBe(page, el("BDD-Miss link in gallery"), "absent"));
    });
    await run.scenario("A name nothing carries empties the list", async () => {
      await session.step(35, "When user enters \"zzz-no-such-space\" into gallery search", () => enterInto(page, "zzz-no-such-space", el("gallery search")));
      await session.step(36, "Then BDD-Find link in gallery should be absent", () => shouldBe(page, el("BDD-Find link in gallery"), "absent"));
      await session.step(37, "And BDD-Miss link in gallery should be absent", () => shouldBe(page, el("BDD-Miss link in gallery"), "absent"));
    });
    await run.scenario("Clearing the search brings both back", async () => {
      await session.step(40, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(41, "Then BDD-Find link in gallery should be visible", () => shouldBe(page, el("BDD-Find link in gallery"), "visible"));
      await session.step(42, "And BDD-Miss link in gallery should be visible", () => shouldBe(page, el("BDD-Miss link in gallery"), "visible"));
    });
    await run.scenario("A child space is searchable inside its parent", async () => {
      await session.step(45, "Given a space named \"BDD-Find-Child\" under \"BDD-Find\" is on the server", () => childSpaceOnServer(page, "BDD-Find-Child", "BDD-Find"));
      await session.step(46, "When user double-clicks on Spaces---BDD-Find tree node inside browse tree", () => doubleClickOn(page, el("Spaces---BDD-Find tree node inside browse tree")));
      await session.step(47, "Then BDD-Find-Child link in gallery should be visible", () => shouldBe(page, el("BDD-Find-Child link in gallery"), "visible"));
      await session.step(48, "When user enters \"zzz-no-such-space\" into gallery search", () => enterInto(page, "zzz-no-such-space", el("gallery search")));
      await session.step(49, "Then BDD-Find-Child link in gallery should be absent", () => shouldBe(page, el("BDD-Find-Child link in gallery"), "absent"));
      await session.step(50, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(51, "Then BDD-Find-Child link in gallery should be visible", () => shouldBe(page, el("BDD-Find-Child link in gallery"), "visible"));
    });
    run.finish();
  });
});
