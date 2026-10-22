import type { AnalysisModule } from '../../lib/store';

export const MODULES = [
  { id: 'overview', label: 'Overview' },
  { id: 'coverage', label: 'Coverage & QC' },
  { id: 'evidence', label: 'Molecules & phasing' },
  { id: 'architectures', label: 'Architectures' },
  { id: 'rearrangements', label: 'Rearrangements' },
  { id: 'haplogroups', label: 'Haplogroups' },
  { id: 'clinical', label: 'Annotations' },
  { id: 'protein', label: 'Protein' },
  { id: 'method', label: 'Method' }
] as const satisfies ReadonlyArray<{ id: AnalysisModule; label: string }>;

/** Maps URL values, including retired module IDs, onto current modules. */
export function parseAnalysisModule(value: string | null): AnalysisModule | undefined {
  if (value === 'variants') return 'overview';
  return MODULES.find((module) => module.id === value)?.id;
}

export function ModuleNavigation({ active, onSelect }: { active: AnalysisModule; onSelect: (module: AnalysisModule) => void }) {
  return (
    <nav aria-label="Analysis sections" className="tabbar">
      {MODULES.map(({ id, label }) => (
        <button
          key={id}
          type="button"
          onClick={() => onSelect(id)}
          aria-current={active === id ? 'page' : undefined}
          className="tab"
        >
          {label}
        </button>
      ))}
    </nav>
  );
}
