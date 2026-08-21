use axum::extract::{Query, State};
use axum::response::Json;
use crate::protocol::strip_engine_suffix;
use crate::routes::AdvertiseQuery;
use crate::state::AppState;

/// GET /conn — return the descriptors of the providers this service actually serves.
///
/// The list comes from the provider registry built in `main.rs`, so it advertises exactly
/// what this service can natively execute — same set as `/info`. Datlas polls every
/// configured endpoint and merges capabilities: a routed type that this endpoint does not
/// advertise falls through to the default (Java) endpoint instead of being pinned here and
/// failing at query time.
///
/// Historically this returned the full Java descriptor list verbatim so a lone Datlas
/// pointed here would still see every data source — but that made the service lie about
/// what it serves, silently hijacking types (e.g. Postgres) it has no driver for.
///
/// `?stripEngineSuffix=true` advertises the base type names — see [`AdvertiseQuery`].
pub async fn list_providers(
    State(state): State<AppState>,
    Query(opts): Query<AdvertiseQuery>,
) -> Json<serde_json::Value> {
    let descriptors: Vec<serde_json::Value> = state
        .registry
        .descriptors()
        .into_iter()
        .map(|mut desc| {
            if opts.strip_engine_suffix {
                desc.source_type = strip_engine_suffix(&desc.source_type).to_string();
            }
            desc
        })
        .filter_map(|desc| serde_json::to_value(desc).ok())
        .collect();
    Json(serde_json::Value::Array(descriptors))
}
