fn main() {
    // The native Rust Databricks driver compiles as a normal Cargo dependency (git dep).
    // No special build script linking needed for the primary path.
    //
    // If using the dbc/FFI fallback for Databricks:
    //   - Set ADBC_DATABRICKS_DRIVER_PATH to the .so location
    //   - Uncomment: println!("cargo:rustc-link-search={}", driver_path);
    //
    // For adbc_snowflake (CGo), the crate's own build.rs handles Go compilation
    // when the "bundled" feature is enabled.

    // Rerun if the build script itself changes
    println!("cargo:rerun-if-changed=build.rs");
}
