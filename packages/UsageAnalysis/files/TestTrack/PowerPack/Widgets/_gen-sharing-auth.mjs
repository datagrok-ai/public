// One-off helper: sign in the non-admin sharing user and save e2e/.auth-sharing.json.
// Run from playwright-tests/: node e2e/Widgets/_gen-sharing-auth.mjs
import { chromium } from '@playwright/test';
import * as dotenv from 'dotenv';

dotenv.config();

const url = process.env.DATAGROK_URL;
const login = process.env.DATAGROK_SHARING_LOGIN;
const password = process.env.DATAGROK_SHARING_PASSWORD;

const browser = await chromium.launch();
const page = await browser.newPage();
await page.goto(url, { waitUntil: 'domcontentloaded', timeout: 60_000 });
const loginField = page.locator('#signup-login-fields input[placeholder="Login or Email"]');
await loginField.waitFor({ timeout: 60_000 });
await loginField.fill(login);
await page.locator('#signup-login-fields input[placeholder="Password"]').fill(password);
await page.locator('#signup-login-fields .signup-buttons button').click();
await page.waitForSelector('.d4-ribbon', { timeout: 90_000 });
await page.context().storageState({ path: 'e2e/.auth-sharing.json' });
await browser.close();
console.log('saved e2e/.auth-sharing.json for', login);
