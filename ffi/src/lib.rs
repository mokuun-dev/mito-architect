//! Safe synchronous ownership boundary for the native mtDNA engine.
mod engine;
mod error;
mod options;
mod raw;

pub use engine::MitoEngine;
pub use error::MitoError;
pub use options::AnalyzeOptions;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct EngineCapabilities {
    pub engine_version: String,
    pub schema_version: String,
    pub error_schema_version: String,
    pub htslib: bool,
}

#[cfg(test)]
mod tests;
