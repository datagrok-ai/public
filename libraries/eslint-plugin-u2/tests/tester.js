'use strict';

const {describe, it} = require('node:test');
const eslint = require('eslint');
const parser = require('@typescript-eslint/parser');

const RuleTester = require('eslint/use-at-your-own-risk').FlatRuleTester ?? eslint.RuleTester;
RuleTester.describe = describe;
RuleTester.it = it;
RuleTester.itOnly = it.only;

module.exports = new RuleTester({languageOptions: {parser, ecmaVersion: 'latest', sourceType: 'module'}});
