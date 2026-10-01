use axum::http::header;
use axum::http::HeaderValue;
use axum::http::Method;
use tower_http::cors::CorsLayer;

pub(super) const DEFAULT_MAX_UPLOAD_BYTES: usize = 64 * 1024 * 1024;

pub(super) fn positive_env_usize(name: &str) -> anyhow::Result<Option<usize>> {
    let value = match std::env::var(name) {
        Ok(value) => value,
        Err(std::env::VarError::NotPresent) => return Ok(None),
        Err(error) => return Err(error.into()),
    };
    let parsed = value
        .parse::<usize>()
        .map_err(|_| anyhow::anyhow!("{name} must be a positive integer"))?;
    if parsed == 0 {
        anyhow::bail!("{name} must be a positive integer");
    }
    Ok(Some(parsed))
}

pub(super) fn nonnegative_env_usize(name: &str) -> anyhow::Result<Option<usize>> {
    optional_env_integer(name, |value| value.parse::<usize>())
}

pub(super) fn env_u8(name: &str) -> anyhow::Result<Option<u8>> {
    optional_env_integer(name, |value| value.parse::<u8>())
}

pub(super) fn env_u16(name: &str) -> anyhow::Result<Option<u16>> {
    let value = match std::env::var(name) {
        Ok(value) => value,
        Err(std::env::VarError::NotPresent) => return Ok(None),
        Err(error) => return Err(error.into()),
    };
    let parsed = if let Some(hex) = value
        .strip_prefix("0x")
        .or_else(|| value.strip_prefix("0X"))
    {
        u16::from_str_radix(hex, 16)
    } else {
        value.parse::<u16>()
    }
    .map_err(|_| anyhow::anyhow!("{name} must be an unsigned 16-bit integer"))?;
    Ok(Some(parsed))
}

pub(super) fn env_u64(name: &str) -> anyhow::Result<Option<u64>> {
    let value = match std::env::var(name) {
        Ok(value) => value,
        Err(std::env::VarError::NotPresent) => return Ok(None),
        Err(error) => return Err(error.into()),
    };
    let parsed = if let Some(hex) = value
        .strip_prefix("0x")
        .or_else(|| value.strip_prefix("0X"))
    {
        u64::from_str_radix(hex, 16)
    } else {
        value.parse::<u64>()
    }
    .map_err(|_| anyhow::anyhow!("{name} must be an unsigned 64-bit integer"))?;
    Ok(Some(parsed))
}

pub(super) fn optional_env_number(name: &str) -> anyhow::Result<Option<f64>> {
    let value = match std::env::var(name) {
        Ok(value) => value,
        Err(std::env::VarError::NotPresent) => return Ok(None),
        Err(error) => return Err(error.into()),
    };
    let parsed = value
        .parse::<f64>()
        .map_err(|_| anyhow::anyhow!("{name} must be a finite number"))?;
    if !parsed.is_finite() {
        anyhow::bail!("{name} must be a finite number");
    }
    Ok(Some(parsed))
}

pub(super) fn unit_interval_env(name: &str, default: f64) -> anyhow::Result<f64> {
    let value = optional_env_number(name)?.unwrap_or(default);
    if !(0.0..=1.0).contains(&value) {
        anyhow::bail!("{name} must be between 0 and 1");
    }
    Ok(value)
}

pub(super) fn optional_env_integer<T, E>(
    name: &str,
    parse: impl FnOnce(&str) -> Result<T, E>,
) -> anyhow::Result<Option<T>> {
    let value = match std::env::var(name) {
        Ok(value) => value,
        Err(std::env::VarError::NotPresent) => return Ok(None),
        Err(error) => return Err(error.into()),
    };
    parse(&value)
        .map(Some)
        .map_err(|_| anyhow::anyhow!("{name} has an invalid integer value"))
}

pub(super) fn cors_layer() -> anyhow::Result<CorsLayer> {
    let configured = std::env::var("MITO_CORS_ORIGINS")
        .unwrap_or_else(|_| "http://127.0.0.1:5173,http://localhost:5173".to_string());
    let origins = configured
        .split(',')
        .map(str::trim)
        .filter(|origin| !origin.is_empty())
        .map(str::parse::<HeaderValue>)
        .collect::<Result<Vec<_>, _>>()?;
    if origins.is_empty() {
        anyhow::bail!("MITO_CORS_ORIGINS must contain at least one valid origin");
    }
    Ok(CorsLayer::new()
        .allow_origin(origins)
        .allow_methods([Method::GET, Method::POST])
        .allow_headers([header::CONTENT_TYPE, header::AUTHORIZATION]))
}
