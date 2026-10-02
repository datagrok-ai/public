// Libraries the Datagrok client serves at runtime are pinned once, in
// build-config/platform-deps.json, and enter the default catalog ("catalog:") from here.
const platform = require('./build-config/platform-deps.json');

module.exports = {
  hooks: {
    updateConfig(config) {
      const versions = Object.fromEntries(Object.entries(platform).map(([name, d]) => [name, d.version]));
      config.catalogs = {...config.catalogs, default: {...config.catalogs?.default, ...versions}};
      if (config.catalog) config.catalog = {...config.catalog, ...versions};
      return config;
    },
  },
};
