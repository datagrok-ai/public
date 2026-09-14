// One bundle load shared by every unit test file. Run: npm test (node --test scripts/unit).
const {loadBundle} = require('../load-bundle.cjs');

module.exports = loadBundle().DG;
