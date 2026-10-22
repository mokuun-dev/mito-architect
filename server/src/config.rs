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

/// Reads an on/off switch; unset means off.
pub(super) fn env_flag(name: &str) -> anyhow::Result<bool> {
    match std::env::var(name) {
        Ok(value) => parse_flag(name, &value),
        Err(std::env::VarError::NotPresent) => Ok(false),
        Err(error) => Err(error.into()),
    }
}

pub(super) fn parse_flag(name: &str, value: &str) -> anyhow::Result<bool> {
    match value.trim().to_ascii_lowercase().as_str() {
        "" | "0" | "false" | "no" => Ok(false),
        "1" | "true" | "yes" => Ok(true),
        _ => anyhow::bail!("{name} must be 1, true, yes, 0, false, or no"),
    }
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

/// Browser origins allowed to call the API (MITO_CORS_ORIGINS). The same list
/// drives CORS and the server-side Origin check.
pub(super) fn cors_origins() -> anyhow::Result<Vec<HeaderValue>> {
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
    Ok(origins)
}

pub(super) fn cors_layer(origins: Vec<HeaderValue>) -> CorsLayer {
    CorsLayer::new()
        .allow_origin(origins)
        .allow_methods([Method::GET, Method::POST])
        .allow_headers([header::CONTENT_TYPE, header::AUTHORIZATION])
}

/// Host names accepted when the API runs without an API key: loopback names
/// plus MITO_ALLOWED_HOSTS. Rejecting other Host headers blocks DNS
/// rebinding from web pages against the unauthenticated local API.
pub(super) fn allowed_hosts() -> Vec<String> {
    let mut hosts = vec![
        "localhost".to_string(),
        "127.0.0.1".to_string(),
        "[::1]".to_string(),
    ];
    if let Ok(configured) = std::env::var("MITO_ALLOWED_HOSTS") {
        hosts.extend(
            configured
                .split(',')
                .map(|host| host.trim().to_ascii_lowercase())
                .filter(|host| !host.is_empty()),
        );
    }
    hosts
}
