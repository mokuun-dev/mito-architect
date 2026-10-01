use serde_json::Value;

pub(super) fn render_html_report(value: &Value) -> String {
    let sample = value
        .pointer("/metadata/sample")
        .and_then(Value::as_str)
        .unwrap_or("mtDNA");
    let pretty = serde_json::to_string_pretty(value).unwrap_or_else(|_| "{}".to_string());
    format!(
        r#"<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Mito-Architect {sample}</title>
<style>
body {{ margin: 0; font-family: ui-sans-serif, system-ui, sans-serif; background: #08111f; color: #e5edf8; }}
header {{ padding: 24px 28px; border-bottom: 1px solid #263244; }}
main {{ padding: 24px 28px; display: grid; gap: 18px; }}
h1 {{ margin: 0; font-size: 32px; letter-spacing: 0; }}
.grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(190px, 1fr)); gap: 12px; }}
.metric {{ border: 1px solid #263244; border-radius: 8px; padding: 16px; background: #101a2b; }}
.metric span {{ color: #94a3b8; font-size: 13px; }}
.metric strong {{ display: block; font-size: 28px; margin-top: 4px; }}
pre {{ overflow: auto; border: 1px solid #263244; border-radius: 8px; padding: 16px; background: #101a2b; }}
</style>
</head>
<body>
<header><h1>Mito-Architect Report</h1></header>
<main>
<section class="grid">
<div class="metric"><span>Sample</span><strong>{sample}</strong></div>
<div class="metric"><span>Clusters</span><strong>{clusters}</strong></div>
<div class="metric"><span>SVs</span><strong>{svs}</strong></div>
<div class="metric"><span>Complex paths</span><strong>{complex_events}</strong></div>
<div class="metric"><span>Passed reads</span><strong>{passed}</strong></div>
</section>
<pre>{pretty}</pre>
</main>
</body>
</html>"#,
        sample = escape_html(sample),
        clusters = value
            .pointer("/clusters")
            .and_then(Value::as_array)
            .map_or(0, Vec::len),
        svs = value
            .pointer("/svs")
            .and_then(Value::as_array)
            .map_or(0, Vec::len),
        complex_events = value
            .pointer("/complex_events")
            .and_then(Value::as_array)
            .map_or(0, Vec::len),
        passed = value
            .pointer("/filter_stats/passed_reads")
            .and_then(Value::as_u64)
            .unwrap_or(0),
        pretty = escape_html(&pretty)
    )
}

pub(super) fn escape_html(value: &str) -> String {
    value
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&#39;")
}
