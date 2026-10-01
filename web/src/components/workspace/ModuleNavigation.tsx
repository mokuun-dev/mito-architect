import { Activity, Atom, BookOpen, ChartNoAxesCombined, Database, Dna, GitFork, LayoutDashboard, Route, ScanSearch } from 'lucide-react';
import type { AnalysisModule } from '../../lib/store';

export const MODULES = [
  { id: 'overview', label: 'Обзор', Icon: LayoutDashboard },
  { id: 'coverage', label: 'Качество и покрытие', Icon: Activity },
  { id: 'variants', label: 'Варианты', Icon: Dna },
  { id: 'evidence', label: 'Молекулы и связи', Icon: GitFork },
  { id: 'architectures', label: 'Архитектуры', Icon: ChartNoAxesCombined },
  { id: 'rearrangements', label: 'Перестройки', Icon: Route },
  { id: 'haplogroups', label: 'Гаплогруппы', Icon: ScanSearch },
  { id: 'clinical', label: 'Аннотации', Icon: Database },
  { id: 'protein', label: 'Белковый контекст', Icon: Atom },
  { id: 'method', label: 'Метод и источники', Icon: BookOpen }
] as const satisfies ReadonlyArray<{ id: AnalysisModule; label: string; Icon: typeof Dna }>;

export function isAnalysisModule(value: string | null): value is AnalysisModule {
  return MODULES.some((module) => module.id === value);
}

export function ModuleNavigation({ active, onSelect }: { active: AnalysisModule; onSelect: (module: AnalysisModule) => void }) {
  return <nav aria-label="Разделы анализа" className="workspace-nav">
    <div className="workspace-nav-title">Исследование</div>
    <div className="workspace-nav-list">
      {MODULES.map(({ id, label, Icon }) => <button
        key={id}
        type="button"
        onClick={() => onSelect(id)}
        aria-current={active === id ? 'page' : undefined}
        className={'workspace-nav-item' + (active === id ? ' workspace-nav-item-active' : '')}
      ><Icon className="h-[18px] w-[18px] shrink-0" aria-hidden /><span>{label}</span></button>)}
    </div>
  </nav>;
}
