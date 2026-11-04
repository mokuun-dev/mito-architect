/** Summaries of a residue's surroundings, computed from coordinates at load. */

export interface AtomPosition {
  x: number;
  y: number;
  z: number;
}

export interface NeighborAtom extends AtomPosition {
  chain: string;
  residueNumber: number;
  residueName: string;
  subunit?: string;
}

export interface ResidueContact {
  chain: string;
  residueNumber: number;
  residueName: string;
  subunit?: string;
  distance: number;
  /** A neighbor in the same chain within two residues is bonded context. */
  sequential: boolean;
  interChain: boolean;
}

export function atomDistance(a: AtomPosition, b: AtomPosition): number {
  const dx = a.x - b.x;
  const dy = a.y - b.y;
  const dz = a.z - b.z;
  return Math.sqrt(dx * dx + dy * dy + dz * dz);
}

/**
 * One entry per neighboring residue with its closest heavy-atom distance. Sorted
 * by distance, then chain and number, so the order does not depend on atom order.
 */
export function summarizeContacts(
  target: { chain: string; residueNumber: number; atoms: AtomPosition[] },
  neighbors: NeighborAtom[],
  cutoff = 5
): ResidueContact[] {
  const closest = new Map<string, ResidueContact>();
  for (const atom of neighbors) {
    if (atom.chain === target.chain && atom.residueNumber === target.residueNumber) continue;
    let distance = Number.POSITIVE_INFINITY;
    for (const targetAtom of target.atoms) {
      distance = Math.min(distance, atomDistance(atom, targetAtom));
    }
    if (distance > cutoff) continue;
    const key = `${atom.chain}:${atom.residueNumber}`;
    const existing = closest.get(key);
    if (!existing || distance < existing.distance) {
      const interChain = atom.chain !== target.chain;
      closest.set(key, {
        chain: atom.chain,
        residueNumber: atom.residueNumber,
        residueName: atom.residueName,
        subunit: atom.subunit,
        distance,
        interChain,
        sequential: !interChain && Math.abs(atom.residueNumber - target.residueNumber) <= 2
      });
    }
  }
  return [...closest.values()].sort(
    (a, b) => a.distance - b.distance || a.chain.localeCompare(b.chain) || a.residueNumber - b.residueNumber
  );
}

export type PackingClass = 'densely packed' | 'intermediate' | 'loosely packed';

/**
 * C-alpha atoms within 10 Å measure how crowded a residue is. It is a rough
 * packing estimate with fixed bands, not a solvent-accessibility calculation,
 * and in a membrane protein a low value may still face lipid.
 */
export function packingClass(neighborCount: number): PackingClass {
  if (neighborCount >= 22) return 'densely packed';
  if (neighborCount >= 14) return 'intermediate';
  return 'loosely packed';
}

export function secondaryStructureName(code: string | undefined): string {
  switch (code) {
    case 'h':
    case 'g':
    case 'i':
      return 'helix';
    case 'e':
    case 'b':
      return 'strand';
    case 't':
    case 's':
    case 'l':
    case '':
    case undefined:
      return 'loop or turn';
    default:
      return 'unassigned';
  }
}

export function percentileRank(value: number, values: number[]): number {
  if (values.length === 0) return 0;
  let below = 0;
  for (const item of values) {
    if (item < value) below += 1;
  }
  return below / values.length;
}

export function contactsToTsv(contacts: ResidueContact[]): string {
  const header = ['chain', 'subunit', 'residue', 'number', 'distance_angstrom', 'relation'];
  const rows = contacts.map((contact) => [
    contact.chain,
    contact.subunit ?? '',
    contact.residueName,
    String(contact.residueNumber),
    contact.distance.toFixed(2),
    contact.interChain ? 'other_chain' : contact.sequential ? 'sequence_neighbor' : 'same_chain'
  ]);
  return [header, ...rows].map((row) => row.join('\t')).join('\n') + '\n';
}
