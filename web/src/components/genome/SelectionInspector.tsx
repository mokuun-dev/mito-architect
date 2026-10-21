import type { MitoAnalysisData } from '@mito-architect/visualization-lib';
import { circularVariants } from '@mito-architect/visualization-lib';
import { ArrowRight, X } from 'lucide-react';
import { useMemo } from 'react';
import { useMitoStore } from '../../lib/store';

export function SelectionInspector({ data }: { data: MitoAnalysisData }) {
  const eventId = useMitoStore((state) => state.selectedEventId);
  const svId = useMitoStore((state) => state.selectedSvId);
  const geneName = useMitoStore((state) => state.selectedGene);
  const moleculeId = useMitoStore((state) => state.selectedMoleculeId);
  const clusterId = useMitoStore((state) => state.selectedCluster);
  const region = useMitoStore((state) => state.selectedRegion);
  const clear = useMitoStore((state) => state.clearInspection);
  const open = useMitoStore((state) => state.setActiveModule);
  const variants = useMemo(() => circularVariants(data), [data]);
  const variant = variants.find((item) => item.id === eventId);
  const sv = data.svs.find((item) => item.id === svId);
  const gene = data.genes.find((item) => item.name === geneName);
  const cluster = data.clusters.find((item) => item.id === clusterId);
  const molecule = data.molecules?.find((item) => item.id === moleculeId);
  const title = variant ? `${variant.position} ${variant.ref}→${variant.alt}` : sv ? `${sv.type} ${sv.start}–${sv.end}` : gene ? gene.name : moleculeId ? moleculeId : cluster ? cluster.label : region ? region.label ?? 'Выбранный участок' : 'Объект не выбран';
  const hasSelection = Boolean(eventId || svId || geneName || moleculeId || clusterId);

  return <section className="workspace-section" aria-label="Инспектор объекта" aria-live="polite">
    <div className="flex items-start justify-between gap-2 border-b border-line pb-3">
      <div className="min-w-0"><div className="text-xs font-medium text-muted">Выбранный объект</div><h2 className="mt-1 break-all text-base font-semibold">{title}</h2></div>
      <button type="button" onClick={clear} disabled={!hasSelection} aria-label="Снять выбор объекта" className="grid h-9 w-9 shrink-0 place-items-center rounded-md border border-line text-muted hover:border-aqua hover:text-aqua disabled:opacity-40"><X className="h-4 w-4" aria-hidden /></button>
    </div>
    {!hasSelection && !region && <p className="mt-4 text-sm leading-6 text-muted">Выберите ген или событие на карте. Здесь появятся координаты, подтверждение и переход к исходным данным.</p>}
    {variant && <div className="mt-4 space-y-3 text-sm">
      <Detail label="Тип" value={variant.type} /><Detail label="Координата" value={`${variant.position} bp`} />
      <Detail label="ALT" value={`${variant.support} молекул`} />
      <Detail label="Пригодный знаменатель" value={variant.callableDepth === undefined ? 'Не рассчитан' : `${variant.callableDepth} молекул`} />
      <Detail label="Доля ALT" value={variant.callableDepth === undefined || variant.callableDepth === 0 || variant.frequency === undefined ? 'Не рассчитана' : `${(variant.frequency * 100).toFixed(1)}%`} />
      {variant.callableDepth === undefined && <p className="border-l-2 border-amber pl-3 text-xs leading-5 text-muted">У результата нет пригодного знаменателя для этой позиции; поддержка события не равна его частоте в образце.</p>}
      {data.metadata.schema_version === '0.6' ? <Action onClick={() => open('evidence')}>Открыть молекулы</Action> : <p className="text-xs leading-5 text-muted">Связь события с молекулярным графом доступна в результате schema 0.6.</p>}<Action onClick={() => open('variants')}>Открыть таблицу вариантов</Action>
    </div>}
    {sv && <div className="mt-4 space-y-3 text-sm"><Detail label="Координаты" value={`${sv.start}–${sv.end} bp`} /><Detail label="Поддержка" value={`${sv.supporting_reads.length} записей`} /><p className="text-xs leading-5 text-muted">Поддержка структурного события не переводится в частоту без подходящего знаменателя отрицательных наблюдений.</p><Action onClick={() => open('rearrangements')}>Открыть перестройку</Action></div>}
    {gene && <div className="mt-4 space-y-3 text-sm"><Detail label="Диапазон" value={`${gene.start}–${gene.end} bp`} /><Detail label="Направление" value={gene.strand} /><Detail label="Тип" value={gene.biotype} /><Detail label="События в гене" value={variants.filter((item) => item.position >= gene.start && item.position <= gene.end).length.toString()} /><Action onClick={() => open('variants')}>Открыть варианты</Action></div>}
    {moleculeId && <div className="mt-4 space-y-3 text-sm"><Detail label="ID" value={moleculeId} /><Detail label="Найденная молекула" value={molecule ? 'Да' : 'В текущем результате не найдена'} /><Action onClick={() => open('evidence')}>Открыть молекулярные данные</Action></div>}
    {cluster && <div className="mt-4 space-y-3 text-sm"><Detail label="Записи в группе" value={cluster.size.toString()} /><p className="text-xs text-muted">Базовая группа по признакам; не является доказанным клеточным клоном.</p><Action onClick={() => open('architectures')}>Открыть сравнение групп</Action></div>}
    {!hasSelection && region && <div className="mt-4 text-sm"><Detail label="Диапазон" value={`${region.start}–${region.end} bp`} /><p className="mt-3 text-xs text-muted">Линейное окно выбранного участка показано под картой.</p></div>}
    {hasSelection && !variant && !sv && !gene && !moleculeId && !cluster && <p className="mt-4 text-sm text-amber">Этот объект отсутствует в загруженном результате.</p>}
  </section>;
}

function Detail({ label, value }: { label: string; value: string }) { return <div className="grid grid-cols-[minmax(0,1fr)_auto] items-start gap-2 border-b border-line/60 pb-2"><span className="text-muted">{label}</span><span className="max-w-44 break-all text-right font-mono tabular-nums">{value}</span></div>; }
function Action({ children, onClick }: { children: string; onClick: () => void }) { return <button type="button" onClick={onClick} className="flex min-h-10 w-full items-center justify-between border-b border-line text-left text-sm font-medium text-aqua hover:underline">{children}<ArrowRight className="h-4 w-4" aria-hidden /></button>; }
