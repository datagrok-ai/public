'use strict';

const {name, version} = require('./package.json');

const rules = {
  'import-entry': require('./rules/import-entry'),
  'no-internal-selectors': require('./rules/no-internal-selectors'),
  'build-in-owner': require('./rules/build-in-owner'),
  'no-raw-subscribe': require('./rules/no-raw-subscribe'),
  'no-unowned-root': require('./rules/no-unowned-root'),
  'no-subclass': require('./rules/no-subclass'),
  'no-ui-in-class': require('./rules/no-ui-in-class'),
  'comments': require('./rules/comments'),
  'no-snapshot-in-live-slot': require('./rules/no-snapshot-in-live-slot'),
  'reason-text': require('./rules/reason-text'),
  'no-raw-column-list': require('./rules/no-raw-column-list'),
};

const plugin = {meta: {name, version}, rules, configs: {}};

// 'warn' is report mode (the fix lands with a later API; consumers-baseline.json ratchets the count), 'error' fails.
plugin.configs.consumer = {
  plugins: {u2: plugin},
  rules: {
    'u2/import-entry': 'warn',
    'u2/no-internal-selectors': 'warn',
    'u2/no-raw-subscribe': 'warn',
    'u2/no-unowned-root': 'warn',
    'u2/no-subclass': 'warn',
    'u2/no-ui-in-class': 'warn',
    'u2/build-in-owner': 'warn',
    'u2/comments': 'warn',
    'u2/no-snapshot-in-live-slot': 'warn',
    'u2/reason-text': 'warn',
    'u2/no-raw-column-list': 'warn',
  },
};

plugin.configs.library = {
  plugins: {u2: plugin},
  rules: {
    'u2/import-entry': 'warn',
    'u2/build-in-owner': 'error',
    'u2/no-unowned-root': 'warn',
    'u2/comments': 'warn',
    'u2/no-snapshot-in-live-slot': 'warn',
    'u2/reason-text': 'warn',
  },
};

module.exports = plugin;
