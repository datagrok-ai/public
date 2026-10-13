/* A package function the UI offers no entry to, called the way the platform calls it — for what it then shows
   (a dialog a package opens only from code). What a function returns is a package test, not a feature. */
import {Page} from '@playwright/test';
import {When} from '../../src/registry.js';
import {callFunction} from '../../src/runtime/functions.js';

export const call = When('user calls {string} function', (page: Page, name: string) => callFunction(page, name),
  {tier: 'api', description: '"Package:function", called with no arguments; a rejected call fails the step with the message of the platform'});
