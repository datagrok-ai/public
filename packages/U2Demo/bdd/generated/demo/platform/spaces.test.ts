/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demo/platform/spaces.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.typeahead, u2.dg.entities]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {openDemoPage} from '../../../bindings/demo.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {selectIn, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {spaceOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Space picker", () => {
  const session = feature(test, "features/demo/platform/spaces.feature", import.meta.url);
  test("A space picked by its name shows in the preview", {tag: ["@demo", "@realizes:u2.typeahead", "@realizes:u2.dg.entities"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(7, "Given user is logged in", () => loggedIn(page));
    await session.step(8, "And a space named \"BDD-U2Demo-Space\" is on the server", () => spaceOnServer(page, "BDD-U2Demo-Space"));
    await session.step(9, "And user opens the \"Spaces\" demo page", () => openDemoPage(page, "Spaces"));
    enter(page, "U2 Demo");
    await session.step(10, "When user types \"BDD-U2Demo\" into Space typeahead", () => typeInto(page, "BDD-U2Demo", el("Space typeahead")));
    await session.step(11, "And user selects \"BDD-U2Demo-Space\" in Space typeahead", () => selectIn(page, "BDD-U2Demo-Space", el("Space typeahead")));
    await session.step(12, "Then \"BDD-U2Demo-Space\" entity card should be visible", () => shouldBe(page, el("\"BDD-U2Demo-Space\" entity card"), "visible"));
    await session.step(13, "And \"BDD-U2Demo-Space\" chip should be visible", () => shouldBe(page, el("\"BDD-U2Demo-Space\" chip"), "visible"));
  });
});
