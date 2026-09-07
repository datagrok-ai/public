use axum::extract::{Query, State};
use axum::response::Json;
use serde_json::{json, Value};
use crate::protocol::strip_engine_suffix;
use crate::routes::AdvertiseQuery;
use crate::state::AppState;

/// GET /health — simple health check.
pub async fn health() -> &'static str {
    "OK"
}

/// GET /info — service metadata.
///
/// `providers` is read from the registry, so it always matches what `/conn` advertises
/// and what the service can actually execute. Claiming a type this endpoint cannot serve
/// silently hijacks it from the Java default endpoint.
///
/// Takes the same `?stripEngineSuffix=true` option as `/conn`, so the two keep agreeing on
/// the names — see [`AdvertiseQuery`].
pub async fn info(State(state): State<AppState>, Query(opts): Query<AdvertiseQuery>) -> Json<Value> {
    let mut providers: Vec<&str> = state
        .registry
        .all()
        .keys()
        .map(|k| if opts.strip_engine_suffix { strip_engine_suffix(k) } else { k.as_str() })
        .collect();
    providers.sort();
    Json(json!({
        "name": "grok_connect_adbc",
        "version": env!("CARGO_PKG_VERSION"),
        "description": "Rust ADBC Connector Service (Arrow IPC)",
        "providers": providers
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::Config;
    use crate::providers::provider_stub::StubProvider;
    use crate::providers::ProviderRegistry;
    use crate::routes::conn::list_providers;
    use crate::state::ConnectionPoolManager;
    use std::sync::Arc;

    /// Stub providers, not real ones: what's under test is that the two endpoints report
    /// whatever the registry holds, which has nothing to do with any particular database.
    fn state_with_providers(names: &[&str]) -> AppState {
        let mut registry = ProviderRegistry::new();
        for name in names {
            registry.register(Arc::new(StubProvider::new(name)));
        }
        AppState {
            registry: Arc::new(registry),
            pool_manager: Arc::new(ConnectionPoolManager::new(1, 300)),
            config: Arc::new(Config::from_env()),
        }
    }

    /// The names `/info` and `/conn` advertise for the same registry, under the same options.
    async fn advertised_names(state: AppState, strip_engine_suffix: bool) -> (Vec<String>, Vec<String>) {
        let Json(info_body) = info(
            State(state.clone()),
            Query(AdvertiseQuery { strip_engine_suffix }),
        ).await;
        let mut from_info: Vec<String> = info_body["providers"]
            .as_array()
            .expect("providers must be an array")
            .iter()
            .map(|v| v.as_str().unwrap().to_string())
            .collect();
        from_info.sort();

        let Json(conn_body) = list_providers(
            State(state),
            Query(AdvertiseQuery { strip_engine_suffix }),
        ).await;
        let mut from_conn: Vec<String> = conn_body
            .as_array()
            .expect("/conn must return an array")
            .iter()
            .map(|d| d["type"].as_str().unwrap().to_string())
            .collect();
        from_conn.sort();
        (from_info, from_conn)
    }

    #[tokio::test]
    async fn info_and_conn_advertise_the_same_providers() {
        let state = state_with_providers(&["AlphaADBC", "BetaADBC", "GammaADBC"]);
        let (from_info, from_conn) = advertised_names(state, false).await;

        assert_eq!(from_info, from_conn);
        assert_eq!(from_info, vec!["AlphaADBC", "BetaADBC", "GammaADBC"]);
    }

    /// With the suffix stripped both endpoints advertise the base type; the engine is still
    /// on the descriptor, which is what tells the variants apart.
    #[tokio::test]
    async fn strip_engine_suffix_drops_the_suffix_from_both_endpoints() {
        let state = state_with_providers(&["AlphaADBC", "BetaADBC"]);
        let (from_info, from_conn) = advertised_names(state.clone(), true).await;

        assert_eq!(from_info, from_conn);
        assert_eq!(from_info, vec!["Alpha", "Beta"]);

        let Json(conn_body) = list_providers(State(state), Query(AdvertiseQuery { strip_engine_suffix: true })).await;
        for desc in conn_body.as_array().unwrap() {
            assert_eq!(desc["engine"].as_str(), Some("ADBC"));
        }
    }

    #[tokio::test]
    async fn info_reports_nothing_when_no_provider_is_registered() {
        let state = AppState {
            registry: Arc::new(ProviderRegistry::new()),
            pool_manager: Arc::new(ConnectionPoolManager::new(1, 300)),
            config: Arc::new(Config::from_env()),
        };
        let Json(body) = info(State(state), Query(AdvertiseQuery::default())).await;
        assert_eq!(body["providers"].as_array().unwrap().len(), 0);
    }
}
