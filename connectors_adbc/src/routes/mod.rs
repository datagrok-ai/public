pub mod health;
pub mod conn;
pub mod test_connection;
pub mod query_socket;
pub mod query;
pub mod schema;
pub mod cancel;

use serde::Deserialize;

/// Query options of the endpoints that advertise providers — `GET /conn` and `GET /info`.
#[derive(Debug, Default, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct AdvertiseQuery {
    /// `?stripEngineSuffix=true` advertises `ClickHouse` instead of `ClickHouseADBC`,
    /// leaving the engine to the descriptor's `engine` field alone.
    ///
    /// Off by default: a datlas that doesn't read `engine` tells the ADBC variant apart from
    /// the Java one by that suffix, and dropping it there would make this service shadow the
    /// Java provider instead of coexisting with it. Only the advertised spelling changes —
    /// queries arriving under either name resolve to the same provider.
    #[serde(default)]
    pub strip_engine_suffix: bool,
}
