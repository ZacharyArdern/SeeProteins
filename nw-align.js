// Needleman-Wunsch sequence alignment + Kabsch superposition + TM-score
import { kabsch } from './kabsch.js';

const AA3TO1 = {
    ALA:'A', ARG:'R', ASN:'N', ASP:'D', CYS:'C', GLN:'Q', GLU:'E',
    GLY:'G', HIS:'H', ILE:'I', LEU:'L', LYS:'K', MET:'M', PHE:'F',
    PRO:'P', SER:'S', THR:'T', TRP:'W', TYR:'Y', VAL:'V',
    MSE:'M', SEC:'U', PYL:'O',
};

function extractResidues(pdbText) {
    const residues = [];
    const seen = new Set();
    for (const line of pdbText.split('\n')) {
        if (!line.startsWith('ATOM') && !line.startsWith('HETATM')) continue;
        if (line.substring(12, 16).trim() !== 'CA') continue;
        const resNum = parseInt(line.substring(22, 26));
        if (seen.has(resNum)) continue;
        seen.add(resNum);
        residues.push({ resNum, aa: AA3TO1[line.substring(17, 20).trim()] ?? 'X' });
    }
    return residues;
}

function extractCACoords(pdbText) {
    const coords = new Map();
    for (const line of pdbText.split('\n')) {
        if (!line.startsWith('ATOM') && !line.startsWith('HETATM')) continue;
        if (line.substring(12, 16).trim() !== 'CA') continue;
        const resNum = parseInt(line.substring(22, 26));
        if (!coords.has(resNum))
            coords.set(resNum, [
                parseFloat(line.substring(30, 38)),
                parseFloat(line.substring(38, 46)),
                parseFloat(line.substring(46, 54)),
            ]);
    }
    return coords;
}

function nw(seqA, seqB, match = 2, mismatch = -1, gap = -2) {
    const m = seqA.length, n = seqB.length;
    const dp = Array.from({ length: m + 1 }, () => new Float32Array(n + 1));
    for (let i = 0; i <= m; i++) dp[i][0] = i * gap;
    for (let j = 0; j <= n; j++) dp[0][j] = j * gap;
    for (let i = 1; i <= m; i++)
        for (let j = 1; j <= n; j++)
            dp[i][j] = Math.max(
                dp[i-1][j-1] + (seqA[i-1] === seqB[j-1] ? match : mismatch),
                dp[i-1][j] + gap,
                dp[i][j-1] + gap,
            );
    let qAln = '', tAln = '', i = m, j = n;
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 &&
            dp[i][j] === dp[i-1][j-1] + (seqA[i-1] === seqB[j-1] ? match : mismatch)) {
            qAln = seqA[i-1] + qAln; tAln = seqB[j-1] + tAln; i--; j--;
        } else if (i > 0 && dp[i][j] === dp[i-1][j] + gap) {
            qAln = seqA[i-1] + qAln; tAln = '-' + tAln; i--;
        } else {
            qAln = '-' + qAln; tAln = seqB[j-1] + tAln; j--;
        }
    }
    return { qAln, tAln };
}

export function nwRigidAlign(pdbA, pdbB) {
    const resA = extractResidues(pdbA);
    const resB = extractResidues(pdbB);
    const caA  = extractCACoords(pdbA);
    const caB  = extractCACoords(pdbB);

    const { qAln, tAln } = nw(
        resA.map(r => r.aa).join(''),
        resB.map(r => r.aa).join(''),
    );

    // Collect matched CA pairs: P = source (B), Q = target (A)
    const P = [], Q = [];
    let ai = 0, bi = 0;
    for (let col = 0; col < qAln.length; col++) {
        const qa = qAln[col], ta = tAln[col];
        if (qa !== '-' && ta !== '-') {
            const cA = caA.get(resA[ai].resNum);
            const cB = caB.get(resB[bi].resNum);
            if (cA && cB) { P.push(cB); Q.push(cA); }
        }
        if (qa !== '-') ai++;
        if (ta !== '-') bi++;
    }

    const { R, t } = kabsch(P, Q);

    // RMSD and TM-score (normalised by query length)
    const Lref = resA.length;
    const d0   = Math.max(0.5, 1.24 * Math.pow(Math.max(Lref - 15, 1), 1/3) - 1.8);
    let sumSq = 0, tmSum = 0;
    for (let k = 0; k < P.length; k++) {
        const p = P[k];
        const nx = R[0][0]*p[0]+R[0][1]*p[1]+R[0][2]*p[2]+t[0];
        const ny = R[1][0]*p[0]+R[1][1]*p[1]+R[1][2]*p[2]+t[1];
        const nz = R[2][0]*p[0]+R[2][1]*p[1]+R[2][2]*p[2]+t[2];
        const q = Q[k];
        const dsq = (nx-q[0])**2 + (ny-q[1])**2 + (nz-q[2])**2;
        sumSq += dsq;
        tmSum += 1 / (1 + dsq / (d0 * d0));
    }

    return {
        tmScore: P.length > 0 ? tmSum / Lref : 0,
        rmsd:    P.length > 0 ? Math.sqrt(sumSq / P.length) : 0,
        alnLen:  P.length,
        qAln, tAln,
        R, t,
    };
}
