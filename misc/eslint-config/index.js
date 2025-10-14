module.exports = {
  "env": {
    "browser": true,
    "es2022": true
  },
  "extends": [
    "google"
  ],
  "parser": "@typescript-eslint/parser",
  "parserOptions": {
    "ecmaVersion": 12,
    "sourceType": "module",
    "project": "./tsconfig.json"
  },
  "plugins": [
    "@typescript-eslint"
  ],
  "rules": {
    "@typescript-eslint/no-unused-vars": ["warn", { "varsIgnorePattern": "^(_|ui$|grok$|DG$)", "argsIgnorePattern": "^_"}],
    "no-trailing-spaces": "off",
    "indent": [
      "error",
      2
    ],
    "max-len": [
      "error",
      120
    ],
    "spaced-comment": "off",
    "linebreak-style": "off",
    "@typescript-eslint/explicit-function-return-type": "off",
    "curly": [
      "error",
      "multi-or-nest"
    ],
    "brace-style": [
      "error",
      "1tbs",
      {
        "allowSingleLine": true
      }
    ],
    "block-spacing": 2
  },
  "overrides": [
    {
      "files": ["*.ts", "*.mts", "*.cts", "*.tsx"],
      "rules": {
        "@typescript-eslint/explicit-function-return-type": "error",
        "@typescript-eslint/no-unsafe-return": "error"
      },
    },
  ],
};
