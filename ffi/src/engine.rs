use crate::raw::*;
use crate::{AnalyzeOptions, EngineCapabilities, MitoError};
use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_void};
use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};

pub struct MitoEngine {
    raw: *mut c_void,
}

// SAFETY: AnalysisEngine owns no thread-affine state. The Rust wrapper is not Sync, so one engine
// cannot be called concurrently through safe Rust, and destruction occurs on the receiving thread.
unsafe impl Send for MitoEngine {}

impl MitoEngine {
    pub fn capabilities() -> EngineCapabilities {
        EngineCapabilities {
            engine_version: static_c_string(unsafe { mito_engine_version() }),
            schema_version: static_c_string(unsafe { mito_engine_schema_version() }),
            error_schema_version: static_c_string(unsafe { mito_engine_error_schema_version() }),
            htslib: unsafe { mito_engine_has_htslib() },
        }
    }

    pub fn new() -> Result<Self, MitoError> {
        let raw = unsafe { mito_engine_new() };
        if raw.is_null() {
            return Err(last_error());
        }
        Ok(Self { raw })
    }

    pub fn analyze(
        &self,
        input_path: &Path,
        reference_path: Option<&Path>,
    ) -> Result<String, MitoError> {
        self.analyze_with_options(input_path, reference_path, AnalyzeOptions::default())
    }

    pub fn analyze_with_options(
        &self,
        input_path: &Path,
        reference_path: Option<&Path>,
        options: AnalyzeOptions,
    ) -> Result<String, MitoError> {
        self.analyze_with_cancellation(input_path, reference_path, options, None)
    }

    pub fn analyze_with_cancel_flag(
        &self,
        input_path: &Path,
        reference_path: Option<&Path>,
        options: AnalyzeOptions,
        cancel_flag: &AtomicBool,
    ) -> Result<String, MitoError> {
        self.analyze_with_cancellation(input_path, reference_path, options, Some(cancel_flag))
    }

    fn analyze_with_cancellation(
        &self,
        input_path: &Path,
        reference_path: Option<&Path>,
        options: AnalyzeOptions,
        cancel_flag: Option<&AtomicBool>,
    ) -> Result<String, MitoError> {
        let input = path_to_cstring(input_path, "input")?;
        let reference = match reference_path {
            Some(path) => Some(path_to_cstring(path, "reference")?),
            None => None,
        };

        let reference_ptr = reference
            .as_ref()
            .map_or(std::ptr::null(), |value| value.as_ptr());
        let molecule_id_tag = protocol_tag_to_cstring(&options.molecule_id_tag)?;
        let umi_tag = protocol_tag_to_cstring(&options.umi_tag)?;
        let duplex_tag = protocol_tag_to_cstring(&options.duplex_tag)?;
        let (callback, user_data) = cancel_flag.map_or((None, std::ptr::null_mut()), |flag| {
            (
                Some(cancel_callback as unsafe extern "C" fn(*mut c_void) -> bool),
                flag as *const AtomicBool as *mut c_void,
            )
        });
        let native_options = MitoEngineAnalyzeOptionsV1 {
            struct_size: std::mem::size_of::<MitoEngineAnalyzeOptionsV1>() as u32,
            abi_version: 1,
            filter_numt: options.filter_numt,
            threads: options.threads.max(1),
            min_mapping_quality: options.min_mapping_quality,
            min_base_quality: options.min_base_quality,
            excluded_snp_flags: options.excluded_snp_flags,
            numt_threshold: options.numt_threshold,
            allow_development_tags: options.allow_development_tags,
            emit_evidence_graph: options.emit_evidence_graph,
            max_evidence_observations: options.max_evidence_observations,
            max_phase_links: options.max_phase_links,
            max_phase_work: options.max_phase_work,
            max_phase_molecule_references: options.max_phase_molecule_references,
            max_result_bytes: options.max_result_bytes,
            evidence_page_size: options.evidence_page_size,
            molecule_id_tag: molecule_id_tag.as_ptr(), umi_tag: umi_tag.as_ptr(),
            duplex_tag: duplex_tag.as_ptr(),
            min_architecture_molecules: options.min_architecture_molecules,
            max_candidate_architectures: options.max_candidate_architectures,
            architecture_max_distance: options.architecture_max_distance,
            architecture_ambiguity_margin: options.architecture_ambiguity_margin,
            architecture_min_overlap_fraction: options.architecture_min_overlap_fraction,
            architecture_consensus_fraction: options.architecture_consensus_fraction,
            architecture_optional_fraction: options.architecture_optional_fraction,
            architecture_stability_replicates: options.architecture_stability_replicates,
            architecture_seed: options.architecture_seed,
            should_cancel: callback, cancel_user_data: user_data,
        };
        let result = unsafe {
            mito_engine_analyze_with_options_v1(
                self.raw, input.as_ptr(), reference_ptr, &native_options)
        };

        if result.is_null() {
            return Err(last_error());
        }

        // Guard ownership before any fallible decoding or Rust allocation.
        let result = NativeString(result);
        result.to_utf8()
    }
}

fn static_c_string(value: *const c_char) -> String {
    if value.is_null() {
        return "unknown".to_owned();
    }
    unsafe { CStr::from_ptr(value) }
        .to_string_lossy()
        .into_owned()
}

fn protocol_tag_to_cstring(value: &str) -> Result<CString, MitoError> {
    CString::new(value).map_err(|_| {
        MitoError::new(
            "MITO-E1001",
            "protocol SAM tag contains an interior NUL byte",
        )
    })
}

fn path_to_cstring(path: &Path, label: &str) -> Result<CString, MitoError> {
    #[cfg(unix)]
    {
        use std::os::unix::ffi::OsStrExt;
        let bytes = path.as_os_str().as_bytes();
        std::str::from_utf8(bytes).map_err(|_| {
            MitoError::new("MITO-E1001", format!("{label} path is not valid UTF-8"))
        })?;
        CString::new(bytes).map_err(|_| {
            MitoError::new(
                "MITO-E1001",
                format!("{label} path contains an interior NUL byte"),
            )
        })
    }
    #[cfg(not(unix))]
    {
        let value = path.to_str().ok_or_else(|| {
            MitoError::new("MITO-E1001", format!("{label} path is not valid UTF-8"))
        })?;
        CString::new(value).map_err(|_| {
            MitoError::new(
                "MITO-E1001",
                format!("{label} path contains an interior NUL byte"),
            )
        })
    }
}

unsafe extern "C" fn cancel_callback(user_data: *mut c_void) -> bool {
    if user_data.is_null() {
        return false;
    }
    let flag = &*(user_data as *const AtomicBool);
    flag.load(Ordering::Relaxed)
}

impl Drop for MitoEngine {
    fn drop(&mut self) {
        if !self.raw.is_null() {
            unsafe { mito_engine_delete(self.raw) };
        }
    }
}

// This guard owns only strings returned by a successful analyze call.
struct NativeString(*const c_char);
impl NativeString {
    fn to_utf8(&self) -> Result<String, MitoError> {
        // SAFETY: native analyze returns a live NUL-terminated allocation.
        unsafe { CStr::from_ptr(self.0) }
            .to_str()
            .map(str::to_owned)
            .map_err(|_| {
                MitoError::new("MITO-E9001", "native analysis returned invalid UTF-8 JSON")
            })
    }
}
impl Drop for NativeString {
    fn drop(&mut self) {
        // SAFETY: uniquely owned native allocation; never borrowed error/version data.
        unsafe { mito_engine_free_string(self.0) };
    }
}

fn last_error() -> MitoError {
    let code = static_c_string(unsafe { mito_engine_get_last_error_code() });
    let message = static_c_string(unsafe { mito_engine_get_last_error() });
    MitoError::new(
        if code == "unknown" || code.is_empty() {
            "MITO-E9001".to_owned()
        } else {
            code
        },
        if message == "unknown" || message.is_empty() {
            "unknown mito engine error".to_owned()
        } else {
            message
        },
    )
}
