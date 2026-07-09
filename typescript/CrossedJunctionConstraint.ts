/**
 * CrossedJunctionConstraint  --  PHASE 2 DRAFT for EternaJS.
 *
 * Soft constraint: rewards structures in which every multiway junction (interior
 * multiloop, or the exterior loop when it joins >= 3 stems) is "crossed" -- i.e.
 * a residue bordering the junction participates in a pseudoknot (crossed) pair.
 *
 * Modeled on MinimumCrossedPercentConstraint.ts. This is a faithful port of the
 * reference Python prototype (../python/crossed_junctions.py); the algorithm has
 * been validated there against the OpenKnotAIDesignData targets. Before merging:
 *   - confirm the evaluate() signature against the target EternaJS branch
 *     (dev uses `ConstraintContext`; older branches pass `undoBlocks[]`);
 *   - add a Bitmaps icon and register in the constraint factory / ConstraintClasses;
 *   - port the python/ unit + example checks into the EternaJS test suite.
 *
 * NOTE: EternaJS SecStruct already provides getCrossedPairs() and stems(); only
 * the junction parser below is genuinely new. If preferred, this whole file can be
 * reduced by moving separateLayers()/findJunctions() onto SecStruct as methods.
 */

import EPars from 'eterna/EPars';
import BitmapManager from 'eterna/resources/BitmapManager';
import Bitmaps from 'eterna/resources/Bitmaps';
import SecStruct from 'eterna/rnatypes/SecStruct';
import ConstraintBox, {ConstraintBoxConfig} from '../ConstraintBox';
import Constraint, {BaseConstraintStatus, ConstraintContext} from '../Constraint';

interface CrossedJunctionStatus extends BaseConstraintStatus {
    totalJunctions: number;
    crossedJunctions: number;
}

type BP = [number, number];

interface Junction {
    kind: 'multiloop' | 'exterior';
    nway: number;
    members: Set<number>;
}

// Residues that participate in any crossed pair, from the full structure.
function crossedResidues(pairs: SecStruct): Set<number> {
    const crossed = pairs.getCrossedPairs();
    const res = new Set<number>();
    for (let i = 0; i < crossed.length; i++) {
        if (crossed.pairingPartner(i) >= 0) res.add(i);
    }
    return res;
}

// Split pairs into non-crossing layers EXACTLY as Eterna does. Rather than
// re-deriving the layering (which is easy to get subtly wrong), we reuse
// SecStruct.getParenthesis() -- the same routine filterForPseudoknots() /
// onlyPseudoknots() rely on -- and split its output by bracket character:
// layer 0 = '()' (the backbone), 1 = '[]', 2 = '{}', 3 = '<>', 4+ = letter pairs.
// This guarantees our layers match the game and never diverge.
function separateLayers(pairs: SecStruct): number[][] {
    const dbn = pairs.getParenthesis({pseudoknots: true});
    const charsL = ['(', '[', '{', '<'];
    const charsR = [')', ']', '}', '>'];
    for (let i = 0; i < 26; i++) {
        charsL.push(String.fromCharCode(i + 97)); // a-z open
        charsR.push(String.fromCharCode(i + 65)); // A-Z close
    }

    const layers: number[][] = [];
    for (let d = 0; d < charsL.length; d++) {
        if (!dbn.includes(charsL[d])) {
            if (d < 4) continue;   // a fixed bracket may be unused; keep scanning
            break;                 // no more letter layers in use
        }
        // Keep only this layer's characters, as a plain nested () string.
        const filtered = Array.from(dbn, (c) => (
            c === charsL[d] ? '(' : (c === charsR[d] ? ')' : '.')
        )).join('');
        layers.push(SecStruct.fromParens(filtered, false).pairs);
    }
    return layers;
}

// Immediate contents of the loop spanning the open interval (lo, hi).
function immediateChildren(layer: number[], lo: number, hi: number): [BP[], number[]] {
    const children: BP[] = [];
    const unpaired: number[] = [];
    let k = lo + 1;
    while (k < hi) {
        const p = layer[k];
        if (p === -1) { unpaired.push(k); k += 1; } else if (p > k) { children.push([k, p]); k = p + 1; } else { k += 1; }
    }
    return [children, unpaired];
}

// Qualifying junctions (>=3-way multiloops + >=3-stem exterior loop) in one layer.
function findJunctions(layer: number[]): Junction[] {
    const n = layer.length;
    const out: Junction[] = [];

    for (let i = 0; i < n; i++) {
        const j = layer[i];
        if (j <= i) continue;
        const [children, unpaired] = immediateChildren(layer, i, j);
        if (children.length >= 2) {
            const members = new Set<number>([i, j, ...unpaired]);
            for (const [a, b] of children) { members.add(a); members.add(b); }
            out.push({kind: 'multiloop', nway: children.length + 1, members});
        }
    }

    const [extChildren, extUnpaired] = immediateChildren(layer, -1, n);
    if (extChildren.length >= 3) {
        const members = new Set<number>(extUnpaired);
        for (const [a, b] of extChildren) { members.add(a); members.add(b); }
        out.push({kind: 'exterior', nway: extChildren.length, members});
    }
    return out;
}

export default class CrossedJunctionConstraint extends Constraint<CrossedJunctionStatus> {
    public static readonly NAME = 'CROSSED_JUNCTION';
    // Minimum fraction of qualifying junctions that must be crossed (1 => all).
    public readonly minFraction: number;

    constructor(minFraction = 1) {
        super();
        this.minFraction = minFraction;
    }

    public evaluate(context: ConstraintContext): CrossedJunctionStatus {
        const undoBlock = context.undoBlocks[0];
        const pseudoknots = (undoBlock.targetConditions !== undefined
            && undoBlock.targetConditions['type'] === 'pseudoknot');
        const pairs = undoBlock.getPairs(EPars.DEFAULT_TEMPERATURE, pseudoknots);

        const crossed = crossedResidues(pairs);
        const layers = separateLayers(pairs);

        let total = 0;
        let ok = 0;
        for (const layer of layers) {
            for (const jn of findJunctions(layer)) {
                total += 1;
                for (const r of jn.members) {
                    if (crossed.has(r)) { ok += 1; break; }
                }
            }
        }

        // No qualifying junction => trivially satisfied.
        const fraction = total === 0 ? 1 : ok / total;
        return {
            satisfied: fraction >= this.minFraction,
            totalJunctions: total,
            crossedJunctions: ok
        };
    }

    public getConstraintBoxConfig(status: CrossedJunctionStatus): ConstraintBoxConfig {
        const statText = ConstraintBox.createTextStyle()
            .append(`${status.crossedJunctions}`, {fill: status.satisfied ? 0x00aa00 : 0xaa0000})
            .append(`/${status.totalJunctions}`);

        return {
            satisfied: status.satisfied,
            tooltip: 'Every multiway junction must be crossed by a pseudoknot.',
            statText,
            icon: BitmapManager.getBitmap(Bitmaps.PseudoknotReqIcon), // TODO: bespoke icon
            drawBG: true,
            showOutline: true
        };
    }

    public serialize(): [string, string] {
        return [
            CrossedJunctionConstraint.NAME,
            this.minFraction.toString()
        ];
    }
}
