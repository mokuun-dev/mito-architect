export type StructureFileExtension = 'pdb' | 'cif' | 'bcif';

export interface StructureSourceCandidate {
  url: string;
  ext: StructureFileExtension;
  label: string;
  entryUrl: string;
  chain: string;
  coordinateDescription: string;
  confidenceSemantics: 'bfactor';
}

export interface MitochondrialSubunit {
  gene: string;
  chain: string;
}

export interface ComplexStructure {
  pdbId: string;
  complex: string;
  /** Resolution in Å, as reported by the RCSB entry. */
  resolution: number;
  subunits: MitochondrialSubunit[];
}

/**
 * Human structures covering all 13 mtDNA-encoded proteins. Chains and numbering
 * are verified against the rCRS translation by scripts/build_structure_catalog.py,
 * and a test checks this table against core/data/mt_structure_catalog.tsv.
 */
export const COMPLEX_STRUCTURES: ComplexStructure[] = [
  {
    pdbId: '9I4I',
    complex: 'Complex I NADH dehydrogenase',
    resolution: 2.63,
    subunits: [
      { gene: 'MT-ND1', chain: 's' },
      { gene: 'MT-ND2', chain: 'i' },
      { gene: 'MT-ND3', chain: 'j' },
      { gene: 'MT-ND4L', chain: 'k' },
      { gene: 'MT-ND4', chain: 'r' },
      { gene: 'MT-ND5', chain: 'l' },
      { gene: 'MT-ND6', chain: 'm' }
    ]
  },
  {
    pdbId: '9HZL',
    complex: 'Complex III cytochrome bc1',
    resolution: 2.52,
    subunits: [{ gene: 'MT-CYB', chain: 'J' }]
  },
  {
    pdbId: '9I7U',
    complex: 'Complex IV cytochrome c oxidase',
    resolution: 3.15,
    subunits: [
      { gene: 'MT-CO1', chain: 'A' },
      { gene: 'MT-CO2', chain: 'B' },
      { gene: 'MT-CO3', chain: 'C' }
    ]
  },
  {
    pdbId: '8H9S',
    complex: 'Complex V ATP synthase',
    resolution: 2.53,
    subunits: [
      { gene: 'MT-ATP8', chain: 'Q' },
      { gene: 'MT-ATP6', chain: 'N' }
    ]
  }
];

export function complexStructure(pdbId: string | undefined): ComplexStructure | undefined {
  const normalized = pdbId ? normalizePdbId(pdbId) : undefined;
  return normalized ? COMPLEX_STRUCTURES.find((item) => item.pdbId === normalized) : undefined;
}

export function normalizePdbId(structureId: string): string | undefined {
  const trimmed = structureId.trim();
  return /^[0-9][A-Z0-9]{3}$/i.test(trimmed) ? trimmed.toUpperCase() : undefined;
}

/**
 * Chain IDs reach NGL inside a selection string. Author chains can be several
 * characters long in large cryo-EM entries, so an invalid ID is refused here
 * rather than replaced: falling back to another chain would highlight the wrong
 * protein without any sign of the error.
 */
export function validChainId(chain: string | undefined): string | undefined {
  return chain && /^[A-Za-z0-9]{1,4}$/.test(chain) ? chain : undefined;
}

export function rcsbCandidates(pdbId: string, chain: string): StructureSourceCandidate[] {
  const normalized = normalizePdbId(pdbId);
  if (!normalized) {
    return [];
  }
  const lowerId = normalized.toLowerCase();
  const common = {
    ext: 'bcif' as const,
    entryUrl: `https://www.rcsb.org/structure/${normalized}`,
    chain,
    coordinateDescription: 'experimental cryo-EM coordinates',
    confidenceSemantics: 'bfactor' as const
  };
  return [
    {
      ...common,
      url: `/structures/${lowerId}.bcif`,
      label: `RCSB PDB ${normalized} (local cache)`
    },
    {
      ...common,
      url: `https://models.rcsb.org/${lowerId}.bcif`,
      label: `RCSB PDB ${normalized}`
    },
    {
      ...common,
      url: `https://files.rcsb.org/download/${normalized}.cif`,
      ext: 'cif',
      label: `RCSB PDB ${normalized} (mmCIF fallback)`
    }
  ];
}
