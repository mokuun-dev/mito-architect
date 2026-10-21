import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import test from 'node:test';

async function source(relativePath) {
  return readFile(new URL(relativePath, import.meta.url), 'utf8');
}

test('application shell preserves keyboard, theme, and reduced-motion contracts', async () => {
  const [app, html, css] = await Promise.all([
    source('../App.tsx'),
    source('../../index.html'),
    source('../index.css')
  ]);

  assert.match(app, /href="#main-content"/);
  assert.match(app, /aria-label=\{light \? 'Use dark theme' : 'Use light theme'\}/);
  assert.match(app, /localStorage\.setItem\('mito-theme'/);
  assert.match(html, /localStorage\.getItem\('mito-theme'\)/);
  assert.match(css, /@media \(prefers-reduced-motion: reduce\)/);
  assert.match(css, /:focus-visible/);
});

test('scientific visualizations retain accessible evidence fallbacks', async () => {
  const [resultPage, viewer, modules, upload] = await Promise.all([
    source('../routes/ResultPage.tsx'),
    source('../components/NglProteinViewer.tsx'),
    source('../components/AnalysisModulesPanel.tsx'),
    source('../components/UploadZone.tsx')
  ]);

  assert.match(resultPage, /id="main-content"/);
  assert.match(resultPage, /Exact values remain available in the linked textual projections/);
  assert.match(viewer, /role="img"/);
  assert.match(viewer, /not a pathogenicity verdict/);
  assert.match(modules, /EvidencePanel/);
  assert.match(modules, /VariantEvidenceTable/);
  assert.match(upload, /role="alert"/);
  assert.match(upload, /aria-describedby=/);
});
