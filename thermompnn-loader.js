import * as ort from 'onnxruntime-web';

const BASE = import.meta.env.BASE_URL;
const MODEL_URL = BASE + 'thermompnn.onnx';

// Point ORT to the wasm files we placed in public/
ort.env.wasm.wasmPaths = BASE;
ort.env.wasm.numThreads = 1; // single-threaded — no SharedArrayBuffer needed

const ALPHABET = 'ACDEFGHIKLMNPQRSTVWYX';

const AA3TO1 = {
    ALA:'A', ARG:'R', ASN:'N', ASP:'D', CYS:'C', GLN:'Q', GLU:'E',
    GLY:'G', HIS:'H', ILE:'I', LEU:'L', LYS:'K', MET:'M', PHE:'F',
    PRO:'P', SER:'S', THR:'T', TRP:'W', TYR:'Y', VAL:'V',
    MSE:'M', SEC:'C', PYL:'K', UNK:'X',
};

let _session = null;

async function getSession() {
    if (!_session) {
        _session = await ort.InferenceSession.create(MODEL_URL, {
            executionProviders: ['wasm'],
        });
    }
    return _session;
}

/**
 * Parse backbone atom coordinates (N, CA, C, O) from PDB text.
 * Returns residues in chain order with all four backbone atoms present.
 */
function parseBackbone(pdbText) {
    const map = {};
    for (const line of pdbText.split('\n')) {
        if (!line.startsWith('ATOM  ') && !line.startsWith('ATOM ')) continue;
        const name  = line.substring(12, 16).trim();
        if (!['N', 'CA', 'C', 'O'].includes(name)) continue;
        const resn  = line.substring(17, 20).trim();
        const chain = (line.substring(21, 22).trim()) || 'A';
        const resi  = parseInt(line.substring(22, 26).trim());
        const icode = line.substring(26, 27).trim();
        const x     = parseFloat(line.substring(30, 38));
        const y     = parseFloat(line.substring(38, 46));
        const z     = parseFloat(line.substring(46, 54));
        if (isNaN(resi) || isNaN(x)) continue;
        const aa    = AA3TO1[resn] ?? 'X';
        const key   = `${chain}:${resi}:${icode}`;
        if (!map[key]) map[key] = { chain, resi, aa, N: null, CA: null, C: null, O: null };
        map[key][name] = [x, y, z];
    }
    return Object.values(map)
        .filter(r => r.N && r.CA && r.C && r.O)
        .sort((a, b) => a.resi - b.resi);
}

/**
 * Run ThermoMPNN encoder on a PDB string.
 * Returns { logits: Float32Array[L*21], residues, resNumToIdx }.
 */
export async function computeLogits(pdbText) {
    const residues = parseBackbone(pdbText);
    if (residues.length < 4) throw new Error('ThermoMPNN: fewer than 4 backbone residues found');

    const L = residues.length;
    const X           = new Float32Array(L * 4 * 3);
    const S           = new BigInt64Array(L);
    const mask        = new Float32Array(L).fill(1.0);
    const residueIdx  = new BigInt64Array(L);
    const chainEnc    = new BigInt64Array(L).fill(1n);

    // Build a map from PDB residue number → array index
    const resNumToIdx = {};
    for (let i = 0; i < L; i++) {
        const r = residues[i];
        const base = i * 12;
        X[base + 0] = r.N[0];  X[base + 1] = r.N[1];  X[base + 2] = r.N[2];
        X[base + 3] = r.CA[0]; X[base + 4] = r.CA[1]; X[base + 5] = r.CA[2];
        X[base + 6] = r.C[0];  X[base + 7] = r.C[1];  X[base + 8] = r.C[2];
        X[base + 9] = r.O[0];  X[base + 10] = r.O[1]; X[base + 11] = r.O[2];
        const aaIdx = ALPHABET.indexOf(r.aa);
        S[i] = BigInt(aaIdx >= 0 ? aaIdx : ALPHABET.indexOf('X'));
        residueIdx[i] = BigInt(r.resi);
        resNumToIdx[r.resi] = i;
    }

    const sess = await getSession();
    const result = await sess.run({
        X:              new ort.Tensor('float32', X,          [1, L, 4, 3]),
        S:              new ort.Tensor('int64',   S,          [1, L]),
        mask:           new ort.Tensor('float32', mask,       [1, L]),
        residue_idx:    new ort.Tensor('int64',   residueIdx, [1, L]),
        chain_encoding: new ort.Tensor('int64',   chainEnc,   [1, L]),
    });

    return {
        logits: result.logits.data,  // Float32Array length L*21
        residues,
        resNumToIdx,
    };
}

/**
 * Compute ΔΔG for one mutation given pre-computed logits.
 * position: 0-based index in the residues array
 * wtAA / mutAA: single-letter amino acid codes
 * Convention: positive ΔΔG = destabilising, negative = stabilising
 */
export function ddG(logits, position, wtAA, mutAA) {
    const wi = ALPHABET.indexOf(wtAA);
    const mi = ALPHABET.indexOf(mutAA);
    if (wi < 0 || mi < 0) return null;
    return logits[position * 21 + mi] - logits[position * 21 + wi];
}
