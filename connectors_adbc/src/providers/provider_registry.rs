use crate::protocol::{strip_engine_suffix, DataSourceDescriptor};
use crate::providers::{Advertised, DbProvider};
use std::collections::HashMap;
use std::sync::Arc;
use tracing::{error, info};

/// Registry of available providers.
pub struct ProviderRegistry {
    providers: HashMap<String, Arc<dyn DbProvider>>,
}

impl ProviderRegistry {
    pub fn new() -> Self {
        ProviderRegistry {
            providers: HashMap::new(),
        }
    }

    pub fn register(&mut self, provider: Arc<dyn DbProvider>) {
        self.providers.insert(provider.name().to_string(), provider);
    }

    /// Register a provider that asks to be advertised, but only after proving it can run.
    ///
    /// **Exits the process** when the driver won't load. Booting without it would leave
    /// `/conn` advertising a type we cannot execute, and datlas routes on that advertisement
    /// — the queries would reach us and fail one by one instead of falling to the Java
    /// default. A dead container is the louder, more honest signal. Startup-only for that
    /// reason: use [`ProviderRegistry::register`] anywhere a failure must not be fatal,
    /// tests included.
    ///
    /// A provider that isn't advertised is never probed: its driver is irrelevant to a
    /// service that won't serve it, and demanding one would break startup for everyone else.
    pub fn register_probed<P>(&mut self, provider: P)
    where
        P: Advertised + DbProvider + 'static,
    {
        if !P::ADVERTISED {
            info!(provider = provider.name(), "provider not advertised, skipping");
            return;
        }
        if let Err(e) = provider.probe() {
            error!(provider = provider.name(), error = ?e, "provider probe failed");
            std::process::exit(1);
        }
        info!(provider = provider.name(), "provider probed, driver loaded");
        self.register(Arc::new(provider));
    }

    /// Look a provider up by the type name a query arrived with, case-insensitively.
    ///
    /// Falls back to matching without the engine suffix: datlas sends back whatever name this
    /// service advertised, so a `/conn?stripEngineSuffix=true` endpoint is asked for
    /// `ClickHouse` while the registry is keyed `ClickHouseADBC`. Exact matches win, so a
    /// provider registered under the bare name is never shadowed by a suffixed one.
    pub fn get(&self, name: &str) -> Option<&Arc<dyn DbProvider>> {
        let name_lower = name.to_lowercase();
        for (k, v) in &self.providers {
            if k.to_lowercase() == name_lower {
                return Some(v);
            }
        }
        for (k, v) in &self.providers {
            if strip_engine_suffix(k).to_lowercase() == name_lower {
                return Some(v);
            }
        }
        None
    }

    pub fn descriptors(&self) -> Vec<DataSourceDescriptor> {
        self.providers.values().map(|p| p.descriptor()).collect()
    }

    pub fn all(&self) -> &HashMap<String, Arc<dyn DbProvider>> {
        &self.providers
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::providers::provider_stub::StubProvider;

    fn registry(names: &[&str]) -> ProviderRegistry {
        let mut registry = ProviderRegistry::new();
        for name in names {
            registry.register(Arc::new(StubProvider::new(name)));
        }
        registry
    }

    #[test]
    fn get_resolves_the_advertised_name_with_or_without_the_engine_suffix() {
        let registry = registry(&["AlphaADBC"]);
        for name in &["AlphaADBC", "alphaadbc", "Alpha", "alpha"] {
            assert_eq!(registry.get(name).map(|p| p.name()), Some("AlphaADBC"), "lookup of {}", name);
        }
        assert!(registry.get("Beta").is_none());
    }

    /// A provider registered under the bare name owns it — the suffix-stripped fallback
    /// must not steal a query meant for it.
    #[test]
    fn exact_match_wins_over_the_suffix_stripped_fallback() {
        let registry = registry(&["AlphaADBC", "Alpha"]);
        assert_eq!(registry.get("Alpha").map(|p| p.name()), Some("Alpha"));
        assert_eq!(registry.get("AlphaADBC").map(|p| p.name()), Some("AlphaADBC"));
    }
}
