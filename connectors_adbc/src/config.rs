use tracing::info;

/// Server configuration loaded from environment variables.
/// Vendor credentials are never configured here — they arrive per request in the
/// `FuncCall`. Test/bench credential loading lives in `crate::test_support`.
#[derive(Debug, Clone)]
pub struct Config {
    pub port: u16,
    pub pool_max_size: u32,
    pub pool_idle_timeout_secs: u64,
    pub pool_eviction_interval_secs: u64,
}

impl Config {
    pub fn from_env() -> Self {
        Config {
            port: env_or("GROK_CONNECT_PORT", 1234),
            pool_max_size: env_or("POOL_MAX_SIZE", 30),
            pool_idle_timeout_secs: env_or("POOL_IDLE_TIMEOUT_SECS", 300),
            pool_eviction_interval_secs: env_or("POOL_EVICTION_INTERVAL_SECS", 60),
        }
    }

    pub fn log_config(&self) {
        info!(port = self.port, "Server port");
        info!(max_size = self.pool_max_size, "Pool max size per target");
        info!(idle_timeout = self.pool_idle_timeout_secs, "Pool idle timeout (secs)");
        info!(eviction = self.pool_eviction_interval_secs, "Pool eviction interval (secs)");
    }
}

fn env_or<T: std::str::FromStr>(name: &str, default: T) -> T {
    std::env::var(name).ok().and_then(|s| s.parse().ok()).unwrap_or(default)
}
