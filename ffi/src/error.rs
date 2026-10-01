use std::fmt;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MitoError {
    pub code: String,
    pub message: String,
}

impl MitoError {
    pub fn new(code: impl Into<String>, message: impl Into<String>) -> Self {
        Self {
            code: code.into(),
            message: message.into(),
        }
    }
}

impl fmt::Display for MitoError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "[{}] {}", self.code, self.message)
    }
}

impl std::error::Error for MitoError {}
