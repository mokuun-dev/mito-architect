import { Download, ExternalLink, FileDown, Focus, Info, Maximize2, Pause, Play, RotateCcw } from 'lucide-react';
import { useCallback, useEffect, useMemo, useRef, useState, type KeyboardEvent as ReactKeyboardEvent, type ReactNode } from 'react';
import type { ProteinStructureMapping } from '@mito-architect/visualization-lib';
import type { Stage, StructureComponent } from 'ngl';
import {
  describeSubstitution,
  parseProteinChange,
  residueMatches,
  type ProteinChange
} from '../lib/aminoAcids';
import {
  complexStructure,
  normalizePdbId,
  rcsbCandidates,
  validChainId,
  type StructureSourceCandidate
} from '../lib/proteinStructureSources';
import {
  contactsToTsv,
  packingClass,
  percentileRank,
  secondaryStructureName,
  summarizeContacts,
  atomDistance,
  type AtomPosition,
  type NeighborAtom,
  type ResidueContact
} from '../lib/structureContext';

export interface ProteinViewerVariant {
  key: string;
  label: string;
  gene?: string;
  protein?: string;
  residue?: string;
  frequency: number;
  structure?: ProteinStructureMapping;
}

type RepresentationMode = 'cartoon' | 'tube' | 'backbone' | 'surface';
type ColorMode = 'sstruc' | 'bfactor' | 'hydrophobicity' | 'residueindex' | 'uniform';

interface NglProteinViewerProps {
  variant: ProteinViewerVariant;
  /** Other variants to mark on the same model; clicking one selects it. */
  companions?: ProteinViewerVariant[];
  onSelectVariant?: (key: string) => void;
}

interface LoadedModel {
  pdbId: string;
  label: string;
  entryUrl: string;
  coordinateDescription: string;
  component: StructureComponent;
  /** Residue name by "chain:number", so identity checks need no structure scan. */
  residues: Map<string, string>;
}

interface ResidueContext {
  contacts: ResidueContact[];
  packingNeighbors: number;
  secondaryStructure: string;
  bfactorPercentile?: number;
}

type ResidueCheck =
  | { state: 'none'; message: string }
  | { state: 'verified'; message: string; residueName: string }
  | { state: 'failed'; message: string };

const MODES: Array<{ value: RepresentationMode; label: string }> = [
  { value: 'cartoon', label: 'Cartoon' },
  { value: 'tube', label: 'Tube' },
  { value: 'backbone', label: 'Backbone' },
  { value: 'surface', label: 'Surface' }
];

const COLORS: Array<{ value: ColorMode; label: string }> = [
  { value: 'sstruc', label: 'Secondary structure' },
  { value: 'bfactor', label: 'Experimental B-factor' },
  { value: 'hydrophobicity', label: 'Hydrophobicity' },
  { value: 'residueindex', label: 'N- to C-terminus' },
  { value: 'uniform', label: 'Uniform' }
];

const THEME_BACKGROUND = { light: '#f5f6f8', dark: '#0d1219' };
const HIGHLIGHT_COLOR = '#f97316';
const HETEROPLASMIC_MARKER = '#a78bfa';
const HOMOPLASMIC_MARKER = '#94a3b8';
const CONTEXT_RADIUS = 5;
const ROTATE_STEP = 8;
const ZOOM_STEP = 0.12;

function themeBackground(): string {
  return document.documentElement.classList.contains('light') ? THEME_BACKGROUND.light : THEME_BACKGROUND.dark;
}

function prefersReducedMotion(): boolean {
  return window.matchMedia?.('(prefers-reduced-motion: reduce)').matches ?? false;
}

function hasWebGl(): boolean {
  try {
    const canvas = document.createElement('canvas');
    return Boolean(canvas.getContext('webgl2') ?? canvas.getContext('webgl'));
  } catch {
    return false;
  }
}

export default function NglProteinViewer({ variant, companions = [], onSelectVariant }: NglProteinViewerProps) {
  const frameRef = useRef<HTMLDivElement>(null);
  const hostRef = useRef<HTMLDivElement>(null);
  const stageRef = useRef<Stage>();
  const modelsRef = useRef(new Map<string, Promise<LoadedModel>>());
  const generationRef = useRef(0);
  const companionsRef = useRef<ProteinViewerVariant[]>(companions);
  const selectRef = useRef(onSelectVariant);
  const [mode, setMode] = useState<RepresentationMode>('cartoon');
  const [colorMode, setColorMode] = useState<ColorMode>('sstruc');
  const [showNeighborhood, setShowNeighborhood] = useState(true);
  const [showAssembly, setShowAssembly] = useState(false);
  const [showCompanions, setShowCompanions] = useState(true);
  const [spinning, setSpinning] = useState(false);
  const [clip, setClip] = useState(0);
  const [stageReady, setStageReady] = useState(false);
  const [stageFailure, setStageFailure] = useState<string>();
  const [loadStatus, setLoadStatus] = useState('Initializing WebGL');
  const [loadedModel, setLoadedModel] = useState<LoadedModel>();
  const [hover, setHover] = useState<{ text: string; x: number; y: number }>();
  const [picked, setPicked] = useState<string>();

  const structureId = variant.structure?.structure_id;
  const pdbId = structureId ? normalizePdbId(structureId) : undefined;
  const chain = validChainId(variant.structure?.chain);
  const change = useMemo(() => parseProteinChange(variant.residue), [variant.residue]);
  const residueNumber = variant.structure?.residue_index ?? change?.position;
  const mappedRoles = useMemo(() => complexStructure(pdbId)?.subunits ?? [], [pdbId]);
  companionsRef.current = companions;
  selectRef.current = onSelectVariant;

  const residueCheck = useMemo<ResidueCheck>(() => {
    if (!loadedModel) {
      return { state: 'none', message: loadStatus };
    }
    if (!chain) {
      return { state: 'failed', message: `The mapping names an invalid chain (${variant.structure?.chain ?? 'none'}); nothing is highlighted` };
    }
    if (residueNumber === undefined) {
      return { state: 'none', message: 'Loaded; this variant has no residue coordinate' };
    }
    return checkResidue(loadedModel, residueNumber, chain, change);
  }, [change, chain, loadStatus, loadedModel, residueNumber, variant.structure?.chain]);

  const residueSelection = residueCheck.state === 'verified' ? `${residueNumber}:${chain}` : undefined;

  const context = useMemo(
    () => (loadedModel && residueSelection && chain && residueNumber !== undefined
      ? residueContext(loadedModel, residueNumber, chain)
      : undefined),
    [chain, loadedModel, residueNumber, residueSelection]
  );
  // NGL selections have no distance operator, so the neighborhood is named residue
  // by residue from the distances computed above.
  const contextSelection = useMemo(() => {
    if (!residueSelection || !context) return undefined;
    return [residueSelection, ...context.contacts.map((contact) => `${contact.residueNumber}:${contact.chain}`)].join(' or ');
  }, [context, residueSelection]);

  const companionMarks = useMemo(() => {
    if (!loadedModel || !showCompanions) return { heteroplasmic: [] as string[], homoplasmic: [] as string[] };
    const heteroplasmic: string[] = [];
    const homoplasmic: string[] = [];
    for (const other of companions) {
      const otherChain = validChainId(other.structure?.chain);
      const number = other.structure?.residue_index;
      if (other.key === variant.key || !otherChain || number === undefined) continue;
      if (normalizePdbId(other.structure?.structure_id ?? '') !== loadedModel.pdbId) continue;
      const expected = parseProteinChange(other.residue)?.reference;
      if (expected && !residueMatches(expected, loadedModel.residues.get(`${otherChain}:${number}`))) continue;
      (other.frequency >= 0.95 ? homoplasmic : heteroplasmic).push(`(${number}:${otherChain} and .CA)`);
    }
    return { heteroplasmic, homoplasmic };
  }, [companions, loadedModel, showCompanions, variant.key]);

  // The stage lives as long as the component; models and representations come and go.
  useEffect(() => {
    if (!hostRef.current) return;
    if (!hasWebGl()) {
      setStageFailure('WebGL is not available in this browser. The residue facts below do not depend on it.');
      return;
    }
    let disposed = false;
    let resizeObserver: ResizeObserver | undefined;
    let themeObserver: MutationObserver | undefined;
    const models = modelsRef.current;

    import('ngl')
      .then(({ Stage: StageClass }) => {
        if (disposed || !hostRef.current) return;
        const stage = new StageClass(hostRef.current, {
          backgroundColor: themeBackground(),
          cameraType: 'perspective',
          quality: 'auto',
          sampleLevel: 0,
          tooltip: false
        });
        stageRef.current = stage;
        resizeObserver = new ResizeObserver(() => stage.handleResize());
        resizeObserver.observe(hostRef.current);
        themeObserver = new MutationObserver(() => stage.setParameters({ backgroundColor: themeBackground() }));
        themeObserver.observe(document.documentElement, { attributes: true, attributeFilter: ['class'] });
        stage.viewer.renderer.domElement.addEventListener('webglcontextlost', () => {
          setStageFailure('The graphics context was lost. Reload the page to restore the 3D view.');
        });

        stage.signals.hovered.add((proxy: PickedProxy | undefined) => {
          const atom = proxy?.atom;
          if (!atom) {
            setHover(undefined);
            return;
          }
          setHover({
            text: describeAtom(atom),
            x: proxy.canvasPosition?.x ?? 0,
            y: (hostRef.current?.clientHeight ?? 0) - (proxy.canvasPosition?.y ?? 0)
          });
        });
        stage.signals.clicked.add((proxy: PickedProxy | undefined) => {
          const atom = proxy?.atom;
          if (!atom) return;
          const target = companionsRef.current.find(
            (other) => other.structure?.chain === atom.chainname && other.structure?.residue_index === atom.resno
          );
          if (target && selectRef.current) {
            selectRef.current(target.key);
          } else {
            setPicked(describeAtom(atom));
          }
        });
        setStageReady(true);
      })
      .catch((error: unknown) => setStageFailure(errorMessage(error, 'The 3D library failed to load')));

    return () => {
      disposed = true;
      resizeObserver?.disconnect();
      themeObserver?.disconnect();
      stageRef.current?.dispose();
      stageRef.current = undefined;
      models.clear();
      setStageReady(false);
    };
  }, []);

  // A parsed model is kept for each structure, so moving between variants of one
  // complex is immediate and downloads and parses nothing again.
  useEffect(() => {
    const stage = stageRef.current;
    if (!stage || !stageReady) return;
    const generation = ++generationRef.current;
    setLoadedModel(undefined);
    setPicked(undefined);
    setSpinning(false);
    stage.setSpin(false);

    if (!structureId) {
      setLoadStatus('No curated structure mapping is available for this variant');
      for (const component of stage.compList) component.setVisibility(false);
      return;
    }
    if (!pdbId) {
      setLoadStatus(`Unsupported structure identifier: ${structureId}. Only four-character PDB IDs are enabled.`);
      return;
    }
    setLoadStatus(`Loading ${pdbId}`);
    let pending = modelsRef.current.get(pdbId);
    if (!pending) {
      pending = loadFirstAvailable(stage, rcsbCandidates(pdbId, chain ?? ''), pdbId);
      modelsRef.current.set(pdbId, pending);
      // A failed download is retried the next time this structure is opened.
      pending.catch(() => modelsRef.current.delete(pdbId));
    }
    pending
      .then((model) => {
        if (generation !== generationRef.current) return;
        for (const component of stage.compList) {
          component.setVisibility(component === model.component);
        }
        setLoadedModel(model);
      })
      .catch((error: unknown) => {
        if (generation === generationRef.current) {
          setLoadStatus(errorMessage(error, 'Could not load the curated structure'));
        }
      });
  }, [chain, pdbId, stageReady, structureId]);

  useEffect(() => {
    if (!loadedModel) return;
    const component = loadedModel.component;
    component.removeAllRepresentations();
    addProteinRepresentations(component, {
      mode,
      colorMode,
      chain,
      residueSelection,
      contextSelection,
      showNeighborhood,
      showAssembly,
      labelText: variant.residue ?? variant.label,
      residueNumber,
      marks: companionMarks
    });
  }, [chain, colorMode, companionMarks, contextSelection, loadedModel, mode, residueNumber, residueSelection, showAssembly, showNeighborhood, variant.label, variant.residue]);

  const focusTarget = useCallback(
    (duration: number) => {
      const component = loadedModel?.component;
      if (!component) return;
      if (residueSelection) {
        component.autoView(showNeighborhood && contextSelection ? contextSelection : residueSelection, duration);
      } else if (chain) {
        component.autoView(`:${chain}`, duration);
      } else {
        component.autoView('all', duration);
      }
    },
    [chain, contextSelection, loadedModel, residueSelection, showNeighborhood]
  );

  // Frames the chain, then the residue, so the first view shows where it sits.
  useEffect(() => {
    if (!loadedModel) return;
    const reduced = prefersReducedMotion();
    if (reduced || !residueSelection) {
      focusTarget(0);
      return;
    }
    loadedModel.component.autoView(chain ? `:${chain}` : 'all', 450);
    const timer = window.setTimeout(() => focusTarget(450), 520);
    return () => window.clearTimeout(timer);
    // The camera moves when the model or residue changes, not when context is toggled.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [loadedModel, residueSelection]);

  useEffect(() => {
    stageRef.current?.setSpin(spinning);
  }, [spinning]);

  useEffect(() => {
    stageRef.current?.setParameters({ clipNear: clip });
  }, [clip]);

  function resetView() {
    setClip(0);
    loadedModel?.component.autoView('all', prefersReducedMotion() ? 0 : 500);
  }

  function exportPng() {
    const stage = stageRef.current;
    if (!stage || !loadedModel) return;
    stage
      .makeImage({ factor: 2, antialias: true, trim: false, transparent: false })
      .then((blob) => saveBlob(blob, `${safeFileName(variant.protein ?? variant.gene ?? 'protein')}-${residueNumber ?? 'model'}.png`))
      .catch((error: unknown) => setLoadStatus(errorMessage(error, 'PNG export failed')));
  }

  function exportContacts() {
    if (!context) return;
    saveBlob(
      new Blob([contactsToTsv(context.contacts)], { type: 'text/tab-separated-values' }),
      `${safeFileName(variant.protein ?? variant.gene ?? 'protein')}-${residueNumber ?? 'residue'}-contacts.tsv`
    );
  }

  function toggleFullscreen() {
    if (frameRef.current) stageRef.current?.toggleFullscreen(frameRef.current);
  }

  function onKeyDown(event: ReactKeyboardEvent<HTMLDivElement>) {
    const stage = stageRef.current;
    if (!stage || event.target !== event.currentTarget) return;
    const controls = stage.trackballControls;
    const actions: Record<string, () => void> = {
      ArrowLeft: () => controls.rotate(-ROTATE_STEP, 0),
      ArrowRight: () => controls.rotate(ROTATE_STEP, 0),
      ArrowUp: () => controls.rotate(0, -ROTATE_STEP),
      ArrowDown: () => controls.rotate(0, ROTATE_STEP),
      '+': () => controls.zoom(ZOOM_STEP),
      '=': () => controls.zoom(ZOOM_STEP),
      '-': () => controls.zoom(-ZOOM_STEP),
      f: () => focusTarget(prefersReducedMotion() ? 0 : 400),
      r: resetView,
      s: () => setSpinning((current) => !current)
    };
    const action = actions[event.key];
    if (action) {
      event.preventDefault();
      action();
      stage.viewer.requestRender();
    }
  }

  const statusText = stageFailure ?? residueCheck.message;
  const heading = variant.protein ?? variant.gene ?? 'Protein model';
  const subunit = mappedRoles.find((role) => role.chain === chain);

  return (
    <div className="min-w-0 overflow-hidden rounded-xl border border-line bg-shell shadow-tool">
      <div className="flex flex-col gap-3 border-b border-line bg-panel2 px-4 py-3 lg:flex-row lg:items-center lg:justify-between">
        <div className="min-w-0 text-sm">
          <div className="flex flex-wrap items-center gap-2">
            <div className="truncate font-semibold">{heading}</div>
            {pdbId && <span className="status-pill font-mono">PDB {pdbId}{chain ? ` · chain ${chain}` : ''}</span>}
            <CheckBadge check={residueCheck} />
          </div>
          <div className="mt-0.5 text-xs text-muted" role="status" aria-live="polite">
            {statusText}
          </div>
        </div>
        <div className="data-scroll" role="group" aria-label="Structure representation">
          <div className="flex min-w-max gap-2 pb-1">
            {MODES.map((item) => (
              <TextButton key={item.value} active={mode === item.value} onClick={() => setMode(item.value)}>
                {item.label}
              </TextButton>
            ))}
          </div>
        </div>
      </div>

      <div
        ref={frameRef}
        className="relative bg-shell focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-[-2px] focus-visible:outline-accent"
        tabIndex={0}
        role="group"
        aria-label="3D view controls: arrow keys rotate, plus and minus zoom, F focuses the residue, R resets, S spins"
        onKeyDown={onKeyDown}
      >
        <div
          ref={hostRef}
          role="img"
          aria-label={`Interactive three-dimensional protein structure for ${heading}`}
          aria-describedby="protein-viewer-accessibility-note"
          className="h-[360px] w-full sm:h-[500px]"
        />
        {(!loadedModel || stageFailure) && (
          <div className="pointer-events-none absolute inset-0 grid place-items-center bg-shell/55 px-8 text-center">
            <div className="max-w-sm rounded-xl border border-line bg-panel/90 px-4 py-3 text-sm text-muted shadow-tool">
              {statusText}
            </div>
          </div>
        )}
        {hover && (
          <div
            className="pointer-events-none absolute z-10 max-w-[280px] rounded-md border border-line bg-panel px-2 py-1 text-xs shadow-pop"
            style={{ left: hover.x + 12, top: hover.y + 12 }}
          >
            {hover.text}
          </div>
        )}
      </div>

      <div className="grid gap-3 border-t border-line bg-panel px-4 py-3">
        <div className="flex flex-wrap items-center gap-2">
          <label className="inline-flex h-9 items-center gap-2 text-xs text-muted">
            Color
            <select
              value={colorMode}
              onChange={(event) => setColorMode(event.target.value as ColorMode)}
              className="h-9 rounded-lg border border-line bg-panel2 px-2 text-xs text-text"
            >
              {COLORS.map((item) => <option key={item.value} value={item.value}>{item.label}</option>)}
            </select>
          </label>
          <Toggle checked={showNeighborhood} onChange={setShowNeighborhood}>{`${CONTEXT_RADIUS} Å context`}</Toggle>
          <Toggle checked={showAssembly} onChange={setShowAssembly}>Full complex</Toggle>
          <Toggle checked={showCompanions} onChange={setShowCompanions}>Other variants</Toggle>
          <label className="inline-flex h-9 items-center gap-2 rounded-lg border border-line bg-panel2 px-3 text-xs text-muted">
            Clip front
            <input type="range" min={0} max={70} step={1} value={clip} onChange={(event) => setClip(Number(event.target.value))} aria-label="Clip the front of the model" className="w-24" />
          </label>
          <div className="ml-auto flex flex-wrap gap-2" role="group" aria-label="Protein viewer actions">
            <IconButton label="Focus residue (F)" onClick={() => focusTarget(prefersReducedMotion() ? 0 : 500)} icon={<Focus className="h-4 w-4" />} />
            <IconButton label="Reset view (R)" onClick={resetView} icon={<RotateCcw className="h-4 w-4" />} />
            <IconButton
              label={spinning ? 'Pause spin (S)' : 'Spin (S)'}
              onClick={() => setSpinning((current) => !current)}
              icon={spinning ? <Pause className="h-4 w-4" /> : <Play className="h-4 w-4" />}
            />
            <IconButton label="Fullscreen" onClick={toggleFullscreen} icon={<Maximize2 className="h-4 w-4" />} />
            <IconButton label="Export PNG" onClick={exportPng} icon={<Download className="h-4 w-4" />} />
          </div>
        </div>
        <dl className="flex flex-wrap gap-x-5 gap-y-1 text-xs text-muted">
          <LegendEntry color={HIGHLIGHT_COLOR} label="Selected residue" />
          <LegendEntry color={HETEROPLASMIC_MARKER} label="Other heteroplasmic variant" />
          <LegendEntry color={HOMOPLASMIC_MARKER} label="Other homoplasmic variant" />
          {loadedModel && (
            <div className="flex items-center gap-1">
              <dt className="sr-only">Source</dt>
              <dd>
                {loadedModel.coordinateDescription},{' '}
                <a href={loadedModel.entryUrl} target="_blank" rel="noreferrer" className="inline-flex items-center gap-1 text-accent hover:underline">
                  {loadedModel.label} <ExternalLink className="h-3 w-3" aria-hidden />
                </a>
              </dd>
            </div>
          )}
        </dl>
        {picked && <p className="text-xs text-muted"><span className="font-medium text-text">Selected in the model:</span> {picked}</p>}
      </div>

      <div className="grid gap-px border-t border-line bg-line lg:grid-cols-[minmax(0,1fr)_minmax(0,1fr)]">
        <div className="grid content-start gap-4 bg-panel px-4 py-3">
          <ResidueFacts context={context} check={residueCheck} subunit={subunit?.gene} residueNumber={residueNumber} chain={chain} />
          <SubstitutionFacts change={change} />
        </div>
        <ContactsTable context={context} onExport={exportContacts} />
      </div>

      <div id="protein-viewer-accessibility-note" className="flex items-start gap-2 border-t border-line bg-sky/10 px-4 py-3 text-xs leading-5 text-muted">
        <Info className="mt-0.5 h-4 w-4 shrink-0 text-sky" aria-hidden />
        <span>
          The 3D view is an exploratory structural projection, not a pathogenicity verdict. It shows the reference residue in an experimental model; the variant itself is not modelled.
          Every value shown here is also available as text.
        </span>
      </div>
    </div>
  );
}

function CheckBadge({ check }: { check: ResidueCheck }) {
  if (check.state === 'verified') return <span className="badge-neutral">Residue verified</span>;
  if (check.state === 'failed') return <span className="badge-warning">Residue not highlighted</span>;
  return null;
}

function ResidueFacts({ context, check, subunit, residueNumber, chain }: {
  context?: ResidueContext;
  check: ResidueCheck;
  subunit?: string;
  residueNumber?: number;
  chain?: string;
}) {
  return (
    <section aria-labelledby="residue-facts-title">
      <h4 id="residue-facts-title" className="section-kicker mb-2">Residue in the model</h4>
      {check.state === 'verified' && context ? (
        <dl className="grid grid-cols-2 gap-x-4 gap-y-2 text-[13px]">
          <Fact label="Residue" value={`${check.residueName} ${residueNumber}, chain ${chain}${subunit ? ` (${subunit})` : ''}`} />
          <Fact label="Secondary structure" value={context.secondaryStructure} />
          <Fact label="Packing" value={`${packingClass(context.packingNeighbors)}, ${context.packingNeighbors} Cα within 10 Å`} />
          <Fact label="B-factor in this chain" value={context.bfactorPercentile === undefined ? '—' : `${Math.round(context.bfactorPercentile * 100)}th percentile`} />
        </dl>
      ) : (
        <p className="text-[13px] text-muted">{check.message}</p>
      )}
      {check.state === 'verified' && (
        <p className="mt-2 text-xs leading-5 text-muted">
          Packing is a rough count of nearby residues, not solvent accessibility. In a membrane protein a loosely packed residue may still face lipid.
        </p>
      )}
    </section>
  );
}

function SubstitutionFacts({ change }: { change?: ProteinChange }) {
  if (!change) return null;
  const heading = <h4 id="substitution-title" className="section-kicker mb-2">Protein change</h4>;
  if (change.kind === 'missense' && change.reference && change.alternate) {
    const properties = describeSubstitution(change.reference, change.alternate);
    const notes = [
      properties.introducesProline && 'introduces proline, which rarely fits a helix',
      properties.losesGlycine && 'replaces glycine, which may be needed for backbone flexibility'
    ].filter(Boolean);
    return (
      <section aria-labelledby="substitution-title">
        {heading}
        <dl className="grid grid-cols-2 gap-x-4 gap-y-2 text-[13px]">
          <Fact label="Change" value={`${change.reference.name} to ${change.alternate.name}`} />
          <Fact label="Grantham distance" value={`${properties.grantham} (${properties.granthamClass})`} />
          <Fact label="Hydropathy change" value={signed(properties.hydropathyChange)} />
          <Fact label="Volume change" value={`${signed(properties.volumeChange)} Å³`} />
          <Fact label="Charge change" value={signed(properties.chargeChange)} />
        </dl>
        {notes.length > 0 && <p className="mt-2 text-xs leading-5 text-muted">This change {notes.join(' and ')}.</p>}
      </section>
    );
  }
  const description: Record<ProteinChange['kind'], string> = {
    missense: 'Amino acid substitution.',
    synonymous: 'The amino acid does not change.',
    nonsense: 'A premature stop codon: everything after this residue is lost.',
    frameshift: 'The reading frame shifts from here, so the protein downstream differs or is truncated.',
    inframe: 'Whole codons are gained or lost, and the frame is kept.',
    start_lost: 'The initiator codon is lost.',
    stop_lost: 'The stop codon is lost.',
    other: 'Protein change.'
  };
  return (
    <section aria-labelledby="substitution-title">
      {heading}
      <p className="text-[13px]">{description[change.kind]}</p>
    </section>
  );
}

function ContactsTable({ context, onExport }: { context?: ResidueContext; onExport: () => void }) {
  const rows = context?.contacts.filter((contact) => !contact.sequential) ?? [];
  const interChain = rows.filter((contact) => contact.interChain);
  return (
    <section className="min-w-0 bg-panel px-4 py-3" aria-labelledby="contacts-title">
      <div className="mb-2 flex items-center justify-between gap-2">
        <h4 id="contacts-title" className="section-kicker">{`Contacts within ${CONTEXT_RADIUS} Å`}</h4>
        {context && (
          <button type="button" className="inline-flex items-center gap-1 text-xs font-medium text-accent hover:underline" onClick={onExport}>
            <FileDown className="h-3.5 w-3.5" aria-hidden /> TSV
          </button>
        )}
      </div>
      {!context ? (
        <p className="text-[13px] text-muted">Available once the residue is verified.</p>
      ) : rows.length === 0 ? (
        <p className="text-[13px] text-muted">No contacts beyond the immediate sequence neighbors.</p>
      ) : (
        <>
          <p className="mb-2 text-[13px] text-muted">
            {rows.length} residues or ligands{interChain.length > 0 ? `, ${interChain.length} from other chains` : ', all in this chain'}.
          </p>
          <div className="max-h-52 overflow-auto rounded-lg border border-line scrollbar-thin">
            <table className="w-full text-left text-xs">
              <thead className="sticky top-0 border-b border-line bg-panel2 text-muted">
                <tr>
                  <th scope="col" className="px-2 py-1.5 font-semibold">Residue</th>
                  <th scope="col" className="px-2 py-1.5 font-semibold">Chain</th>
                  <th scope="col" className="px-2 py-1.5 text-right font-semibold">Å</th>
                </tr>
              </thead>
              <tbody className="tabular-nums">
                {rows.map((contact) => (
                  <tr key={`${contact.chain}:${contact.residueNumber}`} className="border-b border-line last:border-0">
                    <th scope="row" className="whitespace-nowrap px-2 py-1 font-mono font-medium">
                      {contact.residueName} {contact.residueNumber}
                    </th>
                    <td className="px-2 py-1">
                      {contact.chain}
                      {contact.subunit && <span className="text-muted"> · {contact.subunit}</span>}
                    </td>
                    <td className="px-2 py-1 text-right">{contact.distance.toFixed(1)}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </>
      )}
    </section>
  );
}

function Fact({ label, value }: { label: string; value: string }) {
  return (
    <div className="min-w-0">
      <dt className="text-xs text-muted">{label}</dt>
      <dd className="mt-0.5 break-words font-medium tabular-nums">{value}</dd>
    </div>
  );
}

function LegendEntry({ color, label }: { color: string; label: string }) {
  return (
    <div className="flex items-center gap-1.5">
      <dt className="sr-only">Marker</dt>
      <dd className="flex items-center gap-1.5">
        <span className="h-2.5 w-2.5 rounded-full" style={{ background: color }} aria-hidden />
        {label}
      </dd>
    </div>
  );
}

function Toggle({ checked, onChange, children }: { checked: boolean; onChange: (value: boolean) => void; children: string }) {
  return (
    <label className="inline-flex h-9 cursor-pointer items-center gap-2 rounded-lg border border-line bg-panel2 px-3 text-xs text-muted hover:border-accent/50">
      <input type="checkbox" checked={checked} onChange={(event) => onChange(event.target.checked)} />
      {children}
    </label>
  );
}

function TextButton({ active, children, onClick }: { active: boolean; children: string; onClick: () => void }) {
  return (
    <button
      type="button"
      onClick={onClick}
      aria-pressed={active}
      className={[
        'min-h-9 rounded-lg border px-3 py-1.5 text-xs font-semibold',
        active ? 'border-accent bg-accent/10 text-accent' : 'border-line bg-panel text-muted hover:border-accent/50 hover:text-text'
      ].join(' ')}
    >
      {children}
    </button>
  );
}

function IconButton({ icon, label, onClick }: { icon: ReactNode; label: string; onClick: () => void }) {
  return (
    <button
      type="button"
      title={label}
      aria-label={label}
      onClick={onClick}
      className="inline-flex h-9 w-9 items-center justify-center rounded-lg border border-line bg-panel text-muted hover:border-accent hover:text-accent"
    >
      {icon}
    </button>
  );
}

interface RepresentationOptions {
  mode: RepresentationMode;
  colorMode: ColorMode;
  chain?: string;
  residueSelection?: string;
  contextSelection?: string;
  showNeighborhood: boolean;
  showAssembly: boolean;
  labelText: string;
  residueNumber?: number;
  marks: { heteroplasmic: string[]; homoplasmic: string[] };
}

function addProteinRepresentations(component: StructureComponent, options: RepresentationOptions) {
  const { mode, colorMode, chain, residueSelection, contextSelection, showNeighborhood, showAssembly, marks } = options;
  const coloring = colorMode === 'uniform'
    ? { colorScheme: 'uniform', colorValue: '#38bdf8' }
    : { colorScheme: colorMode };
  const modelSelection = chain ? `:${chain}` : 'protein';

  if (showAssembly && chain) {
    // Surrounding subunits stay recognizable but recede behind the mapped chain.
    component.addRepresentation('cartoon', {
      sele: `protein and not :${chain}`,
      colorScheme: 'uniform',
      colorValue: '#94a3b8',
      opacity: 0.22,
      side: 'front',
      depthWrite: false
    } as never);
  }
  if (mode === 'surface') {
    component.addRepresentation('cartoon', { sele: modelSelection, ...coloring, opacity: 0.9 } as never);
    component.addRepresentation('surface', {
      sele: modelSelection,
      ...coloring,
      opacity: 0.3,
      side: 'front',
      depthWrite: false,
      probeRadius: 1.4,
      surfaceType: 'av',
      useWorker: true
    } as never);
  } else if (mode === 'cartoon') {
    component.addRepresentation('cartoon', { sele: modelSelection, ...coloring, aspectRatio: 4.5, radiusScale: 0.9 } as never);
  } else if (mode === 'tube') {
    component.addRepresentation('tube', { sele: modelSelection, ...coloring, radiusScale: 0.8 } as never);
  } else {
    component.addRepresentation('backbone', { sele: modelSelection, ...coloring, radiusScale: 0.8 } as never);
  }

  const addMarks = (selections: string[], color: string) => {
    if (selections.length === 0) return;
    component.addRepresentation('spacefill', {
      sele: selections.join(' or '),
      colorScheme: 'uniform',
      colorValue: color,
      radiusType: 'size',
      radiusSize: 1.5
    } as never);
  };
  addMarks(marks.heteroplasmic, HETEROPLASMIC_MARKER);
  addMarks(marks.homoplasmic, HOMOPLASMIC_MARKER);

  if (!residueSelection) return;
  if (showNeighborhood && contextSelection) {
    component.addRepresentation('licorice', {
      sele: `(${contextSelection}) and not (${residueSelection})`,
      colorScheme: 'element',
      radiusScale: 0.6,
      opacity: 0.8
    } as never);
  }
  component.addRepresentation('ball+stick', {
    sele: residueSelection,
    colorScheme: 'uniform',
    colorValue: HIGHLIGHT_COLOR,
    radiusScale: 1.8
  } as never);
  const light = document.documentElement.classList.contains('light');
  component.addRepresentation('label', {
    sele: `${residueSelection} and .CA`,
    labelType: 'text',
    // NGL looks label text up by atom index, not by position in a list.
    labelText: alphaCarbonLabel(component, options.residueNumber, chain, options.labelText),
    labelGrouping: 'atom',
    colorScheme: 'uniform',
    colorValue: light ? '#0f172a' : '#f8fafc',
    fixedSize: true,
    labelSize: 2.2,
    showBackground: true,
    backgroundColor: light ? '#ffffff' : '#0d1219',
    backgroundOpacity: 0.9,
    backgroundMargin: 2,
    yOffset: 1.2,
    zOffset: 12
  } as never);
}

function alphaCarbonLabel(
  component: StructureComponent,
  residueNumber: number | undefined,
  chain: string | undefined,
  text: string
): Record<number, string> {
  const labels: Record<number, string> = {};
  component.structure.eachAtom((atom) => {
    if (atom.resno === residueNumber && atom.chainname === chain && atom.atomname === 'CA') {
      labels[atom.index] = text;
    }
  });
  return labels;
}

/** Fails closed: the position must exist and hold the variant's reference residue. */
function checkResidue(model: LoadedModel, residueNumber: number, chain: string, change: ProteinChange | undefined): ResidueCheck {
  const place = `${residueNumber}:${chain}`;
  const found = model.residues.get(`${chain}:${residueNumber}`);
  if (!found) {
    return { state: 'failed', message: `Residue ${place} has no coordinates in this model; nothing is highlighted` };
  }
  const expected = change?.reference;
  if (expected && !residueMatches(expected, found)) {
    return { state: 'failed', message: `The model has ${found} at ${place}, but the variant expects ${expected.three.toUpperCase()}; nothing is highlighted` };
  }
  return {
    state: 'verified',
    residueName: found,
    message: expected ? `${found} ${place} confirmed against the variant` : `Residue ${place} mapped; no expected residue to compare`
  };
}

function residueContext(model: LoadedModel, residueNumber: number, chain: string): ResidueContext {
  const structure = model.component.structure;
  const target: AtomPosition[] = [];
  let center: AtomPosition | undefined;
  let secondaryStructure = '';
  let ownBfactor: number | undefined;
  structure.eachAtom((atom) => {
    if (atom.resno !== residueNumber || atom.chainname !== chain || atom.element === 'H') return;
    target.push({ x: atom.x, y: atom.y, z: atom.z });
    if (atom.atomname === 'CA') {
      center = { x: atom.x, y: atom.y, z: atom.z };
      secondaryStructure = atom.sstruc;
      ownBfactor = atom.bfactor;
    }
  });
  const origin = center ?? target[0];
  const neighbors: NeighborAtom[] = [];
  const chainBfactors: number[] = [];
  let packingNeighbors = 0;
  const mappedGene = new Map(complexStructure(model.pdbId)?.subunits.map((role) => [role.chain, role.gene]) ?? []);
  if (origin) {
    structure.eachAtom((atom) => {
      if (atom.element === 'H' || atom.isWater()) return;
      const position = { x: atom.x, y: atom.y, z: atom.z };
      const fromCenter = atomDistance(position, origin);
      if (atom.atomname === 'CA' && atom.isProtein()) {
        if (atom.chainname === chain) chainBfactors.push(atom.bfactor);
        if (fromCenter <= 10 && !(atom.resno === residueNumber && atom.chainname === chain)) packingNeighbors += 1;
      }
      // Any contact lies within the cutoff plus the residue's own reach from Cα.
      if (fromCenter <= CONTEXT_RADIUS + 12) {
        neighbors.push({
          ...position,
          chain: atom.chainname,
          residueNumber: atom.resno,
          residueName: atom.resname,
          subunit: mappedGene.get(atom.chainname) ?? entityName(atom)
        });
      }
    });
  }
  return {
    contacts: summarizeContacts({ chain, residueNumber, atoms: target }, neighbors, CONTEXT_RADIUS),
    packingNeighbors,
    secondaryStructure: secondaryStructureName(secondaryStructure),
    bfactorPercentile: ownBfactor === undefined || chainBfactors.length === 0 ? undefined : percentileRank(ownBfactor, chainBfactors)
  };
}

interface PickedAtom {
  resname: string;
  resno: number;
  chainname: string;
  entity?: { description?: string };
}

interface PickedProxy {
  atom?: PickedAtom;
  canvasPosition?: { x: number; y: number };
}

function entityName(atom: { entity?: { description?: string } }): string | undefined {
  const description = atom.entity?.description;
  return description && description.length > 0 ? description : undefined;
}

function describeAtom(atom: PickedAtom): string {
  const name = entityName(atom);
  return `${atom.resname} ${atom.resno}, chain ${atom.chainname}${name ? `, ${name}` : ''}`;
}

async function loadFirstAvailable(stage: Stage, candidates: StructureSourceCandidate[], pdbId: string): Promise<LoadedModel> {
  const failures: string[] = [];
  for (const candidate of candidates) {
    try {
      const component = await stage.loadFile(candidate.url, { ext: candidate.ext, defaultRepresentation: false });
      if (!component || !('structure' in component)) {
        failures.push(`${candidate.label}: no structure returned`);
        continue;
      }
      const structureComponent = component as StructureComponent;
      const residues = new Map<string, string>();
      structureComponent.structure.eachResidue((residue) => {
        residues.set(`${residue.chainname}:${residue.resno}`, residue.resname);
      });
      return {
        pdbId,
        component: structureComponent,
        residues,
        label: candidate.label,
        entryUrl: candidate.entryUrl,
        coordinateDescription: candidate.coordinateDescription
      };
    } catch (error: unknown) {
      failures.push(`${candidate.label}: ${errorMessage(error, 'request failed')}`);
    }
  }
  throw new Error(`Curated RCSB structure unavailable. ${failures.join('; ')}. Restore the bundled cache with scripts/fetch_protein_structures.sh.`);
}

function signed(value: number): string {
  return value > 0 ? `+${value}` : String(value);
}

function saveBlob(blob: Blob, name: string) {
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = name;
  link.click();
  window.setTimeout(() => URL.revokeObjectURL(url), 0);
}

function safeFileName(value: string): string {
  return value.replace(/[^A-Za-z0-9._-]+/g, '_').replace(/^_+|_+$/g, '') || 'protein';
}

function errorMessage(error: unknown, fallback: string): string {
  return error instanceof Error && error.message ? error.message : fallback;
}
