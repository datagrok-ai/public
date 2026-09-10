/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/smoke.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/diff-studio.js';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Diff Studio smoke", () => {
  const session = feature(test, "features/smoke.feature", import.meta.url);
  test("The platform is up", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(6, "Given user is logged in", () => loggedIn(page));
    await session.step(7, "Then browse tab should be visible", () => shouldBe(page, el("browse tab"), "visible"));
  });
});
