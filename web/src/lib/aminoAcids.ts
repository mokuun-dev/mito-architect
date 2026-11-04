/**
 * Amino-acid properties for describing a protein change next to its 3D context.
 * These describe the chemistry of the substitution, not its pathogenicity.
 */

export interface AminoAcid {
  one: string;
  three: string;
  name: string;
  /** Kyte and Doolittle (1982) hydropathy index. */
  hydropathy: number;
  /** Side-chain charge at pH 7. */
  charge: -1 | 0 | 1;
  /** Grantham (1974) composition, polarity and volume. */
  composition: number;
  polarity: number;
  volume: number;
}

const TABLE: AminoAcid[] = [
  { one: 'A', three: 'Ala', name: 'Alanine', hydropathy: 1.8, charge: 0, composition: 0, polarity: 8.1, volume: 31 },
  { one: 'R', three: 'Arg', name: 'Arginine', hydropathy: -4.5, charge: 1, composition: 0.65, polarity: 10.5, volume: 124 },
  { one: 'N', three: 'Asn', name: 'Asparagine', hydropathy: -3.5, charge: 0, composition: 1.33, polarity: 11.6, volume: 56 },
  { one: 'D', three: 'Asp', name: 'Aspartate', hydropathy: -3.5, charge: -1, composition: 1.38, polarity: 13, volume: 54 },
  { one: 'C', three: 'Cys', name: 'Cysteine', hydropathy: 2.5, charge: 0, composition: 2.75, polarity: 5.5, volume: 55 },
  { one: 'Q', three: 'Gln', name: 'Glutamine', hydropathy: -3.5, charge: 0, composition: 0.89, polarity: 10.5, volume: 85 },
  { one: 'E', three: 'Glu', name: 'Glutamate', hydropathy: -3.5, charge: -1, composition: 0.92, polarity: 12.3, volume: 83 },
  { one: 'G', three: 'Gly', name: 'Glycine', hydropathy: -0.4, charge: 0, composition: 0.74, polarity: 9, volume: 3 },
  { one: 'H', three: 'His', name: 'Histidine', hydropathy: -3.2, charge: 0, composition: 0.58, polarity: 10.4, volume: 96 },
  { one: 'I', three: 'Ile', name: 'Isoleucine', hydropathy: 4.5, charge: 0, composition: 0, polarity: 5.2, volume: 111 },
  { one: 'L', three: 'Leu', name: 'Leucine', hydropathy: 3.8, charge: 0, composition: 0, polarity: 4.9, volume: 111 },
  { one: 'K', three: 'Lys', name: 'Lysine', hydropathy: -3.9, charge: 1, composition: 0.33, polarity: 11.3, volume: 119 },
  { one: 'M', three: 'Met', name: 'Methionine', hydropathy: 1.9, charge: 0, composition: 0, polarity: 5.7, volume: 105 },
  { one: 'F', three: 'Phe', name: 'Phenylalanine', hydropathy: 2.8, charge: 0, composition: 0, polarity: 5.2, volume: 132 },
  { one: 'P', three: 'Pro', name: 'Proline', hydropathy: -1.6, charge: 0, composition: 0.39, polarity: 8, volume: 32.5 },
  { one: 'S', three: 'Ser', name: 'Serine', hydropathy: -0.8, charge: 0, composition: 1.42, polarity: 9.2, volume: 32 },
  { one: 'T', three: 'Thr', name: 'Threonine', hydropathy: -0.7, charge: 0, composition: 0.71, polarity: 8.6, volume: 61 },
  { one: 'W', three: 'Trp', name: 'Tryptophan', hydropathy: -0.9, charge: 0, composition: 0.13, polarity: 5.4, volume: 170 },
  { one: 'Y', three: 'Tyr', name: 'Tyrosine', hydropathy: -1.3, charge: 0, composition: 0.2, polarity: 6.2, volume: 136 },
  { one: 'V', three: 'Val', name: 'Valine', hydropathy: 4.2, charge: 0, composition: 0, polarity: 5.9, volume: 84 }
];

const BY_THREE = new Map(TABLE.map((item) => [item.three.toUpperCase(), item]));
const BY_ONE = new Map(TABLE.map((item) => [item.one, item]));

export function aminoAcid(code: string | undefined): AminoAcid | undefined {
  if (!code) return undefined;
  return code.length === 1 ? BY_ONE.get(code.toUpperCase()) : BY_THREE.get(code.toUpperCase());
}

/** Grantham (1974) distance, using the published formula and weights. */
export function granthamDistance(first: AminoAcid, second: AminoAcid): number {
  const composition = first.composition - second.composition;
  const polarity = first.polarity - second.polarity;
  const volume = first.volume - second.volume;
  return Math.round(
    50.723 * Math.sqrt(1.833 * composition * composition + 0.1018 * polarity * polarity + 0.000399 * volume * volume)
  );
}

export type ProteinChangeKind = 'missense' | 'synonymous' | 'nonsense' | 'frameshift' | 'inframe' | 'start_lost' | 'stop_lost' | 'other';

export interface ProteinChange {
  kind: ProteinChangeKind;
  position: number;
  reference?: AminoAcid;
  alternate?: AminoAcid;
}

/** Parses the HGVS protein forms the engine emits, such as p.Leu156Arg. */
export function parseProteinChange(residue: string | undefined): ProteinChange | undefined {
  if (!residue) return undefined;
  const match = /^p\.(?:\(?)([A-Z][a-z]{2})(\d+)([A-Za-z]{3}|=|\?|fs|del|ins|extTer\?)?/.exec(residue.trim());
  if (!match) return undefined;
  const reference = aminoAcid(match[1]);
  const position = Number(match[2]);
  if (!Number.isSafeInteger(position) || position <= 0) return undefined;
  const suffix = match[3] ?? '';
  const tail = residue.trim().slice(match[0].length);
  if (match[1] === 'Ter') {
    return { kind: 'stop_lost', position, alternate: aminoAcid(suffix) };
  }
  if (suffix === '=') return { kind: 'synonymous', position, reference, alternate: reference };
  if (suffix === '?') return { kind: 'start_lost', position, reference };
  if (suffix === 'Ter') return { kind: 'nonsense', position, reference };
  if (suffix === 'fs') return { kind: 'frameshift', position, reference };
  if (suffix === 'del' || suffix === 'ins') return { kind: 'inframe', position, reference };
  if (tail.length > 0) return { kind: 'other', position, reference };
  const alternate = aminoAcid(suffix);
  return alternate ? { kind: 'missense', position, reference, alternate } : { kind: 'other', position, reference };
}

export interface SubstitutionProperties {
  grantham: number;
  /** Conventional bands from Li, Wu and Luo (1984). */
  granthamClass: 'conservative' | 'moderately conservative' | 'moderately radical' | 'radical';
  hydropathyChange: number;
  volumeChange: number;
  chargeChange: number;
  introducesProline: boolean;
  losesGlycine: boolean;
}

export function describeSubstitution(reference: AminoAcid, alternate: AminoAcid): SubstitutionProperties {
  const grantham = granthamDistance(reference, alternate);
  return {
    grantham,
    granthamClass:
      grantham <= 50
        ? 'conservative'
        : grantham <= 100
          ? 'moderately conservative'
          : grantham <= 150
            ? 'moderately radical'
            : 'radical',
    hydropathyChange: Number((alternate.hydropathy - reference.hydropathy).toFixed(1)),
    volumeChange: Number((alternate.volume - reference.volume).toFixed(1)),
    chargeChange: alternate.charge - reference.charge,
    introducesProline: alternate.one === 'P' && reference.one !== 'P',
    losesGlycine: reference.one === 'G' && alternate.one !== 'G'
  };
}

/** Three-letter residue name as it appears in structure files, such as LEU. */
export function structureResidueName(acid: AminoAcid | undefined): string | undefined {
  return acid?.three.toUpperCase();
}

const INITIATOR_NAMES = new Set(['MET', 'FME', 'MSE']);

/**
 * Whether a residue found in coordinates is the residue a variant expects. The
 * mitochondrial initiator can be modelled as N-formylmethionine, so those names
 * are the same residue here.
 */
export function residueMatches(expected: AminoAcid | undefined, found: string | undefined): boolean {
  if (!expected || !found) return false;
  const name = found.toUpperCase();
  return expected.three.toUpperCase() === 'MET' ? INITIATOR_NAMES.has(name) : expected.three.toUpperCase() === name;
}
