// Kabsch alignment + flexible PDB transform — pure JS, no dependencies

// ── Matrix math (3×3) ──────────────────────────────────────────────────────────

function zeros33() { return [[0,0,0],[0,0,0],[0,0,0]]; }
function eye3()    { return [[1,0,0],[0,1,0],[0,0,1]]; }

function mmul(A, B) {
    const C = zeros33();
    for (let i = 0; i < 3; i++)
        for (let k = 0; k < 3; k++)
            for (let j = 0; j < 3; j++)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

function mtr(A) {
    return [[A[0][0],A[1][0],A[2][0]],
            [A[0][1],A[1][1],A[2][1]],
            [A[0][2],A[1][2],A[2][2]]];
}

function det3(A) {
    return A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])
          -A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
          +A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
}

function mv3(M, v) {
    return [
        M[0][0]*v[0]+M[0][1]*v[1]+M[0][2]*v[2],
        M[1][0]*v[0]+M[1][1]*v[1]+M[1][2]*v[2],
        M[2][0]*v[0]+M[2][1]*v[1]+M[2][2]*v[2]
    ];
}

// ── Jacobi eigenvalue for 3×3 symmetric matrix ─────────────────────────────────
// Returns {vals, vecs} where vecs[k] is the k-th eigenvector.

function jacobi3(S) {
    const a = S.map(r => [...r]);
    const V = eye3();

    for (let sweep = 0; sweep < 50; sweep++) {
        let off = 0;
        for (let i = 0; i < 3; i++)
            for (let j = i+1; j < 3; j++) off += a[i][j]*a[i][j];
        if (off < 1e-30) break;

        for (let p = 0; p < 2; p++) {
            for (let q = p+1; q < 3; q++) {
                if (Math.abs(a[p][q]) < 1e-15) continue;
                const tau = (a[q][q]-a[p][p]) / (2*a[p][q]);
                const t   = (tau >= 0 ? 1 : -1) / (Math.abs(tau)+Math.sqrt(1+tau*tau));
                const c   = 1 / Math.sqrt(1+t*t);
                const s   = t * c;
                const app=a[p][p], aqq=a[q][q], apq=a[p][q];
                a[p][p] = c*c*app - 2*s*c*apq + s*s*aqq;
                a[q][q] = s*s*app + 2*s*c*apq + c*c*aqq;
                a[p][q] = a[q][p] = 0;
                for (let r = 0; r < 3; r++) {
                    if (r===p||r===q) continue;
                    const apr=a[p][r], aqr=a[q][r];
                    a[p][r]=a[r][p]= c*apr - s*aqr;
                    a[q][r]=a[r][q]= s*apr + c*aqr;
                }
                for (let r = 0; r < 3; r++) {
                    const vp=V[r][p], vq=V[r][q];
                    V[r][p] = c*vp - s*vq;
                    V[r][q] = s*vp + c*vq;
                }
            }
        }
    }

    return {
        vals: [a[0][0], a[1][1], a[2][2]],
        vecs: [0,1,2].map(k => [V[0][k], V[1][k], V[2][k]])
    };
}

// ── Kabsch algorithm ───────────────────────────────────────────────────────────
// P: source [[x,y,z],...], Q: target [[x,y,z],...]
// Returns {R, t} such that R @ p + t ≈ q for each pair.

export function kabsch(P, Q) {
    const n = P.length;
    if (n < 3) return { R: eye3(), t: [0,0,0] };

    const cP=[0,0,0], cQ=[0,0,0];
    for (let i=0;i<n;i++) for (let k=0;k<3;k++) { cP[k]+=P[i][k]; cQ[k]+=Q[i][k]; }
    for (let k=0;k<3;k++) { cP[k]/=n; cQ[k]/=n; }

    const H = zeros33();
    for (let i=0;i<n;i++)
        for (let j=0;j<3;j++)
            for (let k=0;k<3;k++)
                H[j][k] += (P[i][j]-cP[j]) * (Q[i][k]-cQ[k]);

    const { vals, vecs } = jacobi3(mmul(mtr(H), H));
    const order = [0,1,2].sort((a,b) => vals[b]-vals[a]);
    const s     = order.map(i => Math.sqrt(Math.max(0, vals[i])));

    // Build V matrix (sorted eigenvectors as columns)
    const Vm = zeros33();
    for (let k=0;k<3;k++) for (let r=0;r<3;r++) Vm[r][k] = vecs[order[k]][r];

    // U columns: H V[:,k] / s[k]
    const Um = zeros33();
    for (let k=0;k<3;k++) {
        if (s[k] < 1e-10) continue;
        for (let i=0;i<3;i++) {
            let sum=0;
            for (let j=0;j<3;j++) sum += H[i][j]*Vm[j][k];
            Um[i][k] = sum/s[k];
        }
    }
    // Complete degenerate axis via cross product
    if (s[2] < 1e-10) {
        Um[0][2] = Um[1][0]*Um[2][1]-Um[2][0]*Um[1][1];
        Um[1][2] = Um[2][0]*Um[0][1]-Um[0][0]*Um[2][1];
        Um[2][2] = Um[0][0]*Um[1][1]-Um[1][0]*Um[0][1];
    }

    let R = mmul(Vm, mtr(Um));

    // Fix reflection
    if (det3(R) < 0) {
        for (let r=0;r<3;r++) Vm[r][2] *= -1;
        R = mmul(Vm, mtr(Um));
    }

    const RcP = mv3(R, cP);
    const t   = cQ.map((q,k) => q - RcP[k]);

    return { R, t };
}

// ── PDB utilities ──────────────────────────────────────────────────────────────

export function parseCACoords(pdbText) {
    const out = [];
    for (const line of pdbText.split('\n')) {
        if (!line.startsWith('ATOM') && !line.startsWith('HETATM')) continue;
        if (line.substring(12, 16).trim() !== 'CA') continue;
        const resNum = parseInt(line.substring(22, 26));
        const x = parseFloat(line.substring(30, 38));
        const y = parseFloat(line.substring(38, 46));
        const z = parseFloat(line.substring(46, 54));
        if (!isNaN(resNum) && !isNaN(x)) out.push({ resNum, xyz: [x, y, z] });
    }
    return out;
}

// Split residue list into segments at hinge positions.
// Hinge at residue h means segment boundary before h.
export function defineSegments(commonRes, hinges) {
    const sorted = [...hinges].sort((a,b) => a-b);
    const breaks = [commonRes[0], ...sorted, commonRes[commonRes.length-1]+1];
    return breaks.slice(0,-1).map((lo, i) => {
        const hi  = breaks[i+1];
        const res = commonRes.filter(r => r >= lo && r < hi);
        return { idx: i, residues: res, minRes: res[0], maxRes: res[res.length-1] };
    }).filter(s => s.residues.length >= 3);
}

// Apply per-segment Kabsch transforms to pdbBt atoms.
export function flexibleAlignPDB(pdbA, pdbBt, hingeResidues) {
    const caA  = parseCACoords(pdbA);
    const caBt = parseCACoords(pdbBt);

    const mapA  = new Map(caA.map(c  => [c.resNum, c.xyz]));
    const mapBt = new Map(caBt.map(c => [c.resNum, c.xyz]));

    const common = [...mapA.keys()].filter(r => mapBt.has(r)).sort((a,b) => a-b);
    const segments = defineSegments(common, hingeResidues);

    // Per-segment Kabsch: rotate pdbBt segment onto pdbA
    const transforms = segments.map(seg => {
        const P = seg.residues.map(r => mapBt.get(r)).filter(Boolean);
        const Q = seg.residues.map(r => mapA.get(r)).filter(Boolean);
        return kabsch(P, Q);
    });

    // Build residue → segment index map
    const resToSeg = new Map();
    segments.forEach((seg, si) => seg.residues.forEach(r => resToSeg.set(r, si)));

    // Transform all atoms in pdbBt
    const lines = pdbBt.split('\n').map(line => {
        if (!line.startsWith('ATOM') && !line.startsWith('HETATM')) return line;
        const resNum = parseInt(line.substring(22, 26));
        const si = resToSeg.get(resNum);
        if (si === undefined) return line;
        const { R, t } = transforms[si];
        const x = parseFloat(line.substring(30, 38));
        const y = parseFloat(line.substring(38, 46));
        const z = parseFloat(line.substring(46, 54));
        if (isNaN(x)) return line;
        const nx = R[0][0]*x+R[0][1]*y+R[0][2]*z+t[0];
        const ny = R[1][0]*x+R[1][1]*y+R[1][2]*z+t[1];
        const nz = R[2][0]*x+R[2][1]*y+R[2][2]*z+t[2];
        return line.substring(0,30) +
            nx.toFixed(3).padStart(8) +
            ny.toFixed(3).padStart(8) +
            nz.toFixed(3).padStart(8) +
            line.substring(54);
    });

    // Compute RMSDh after flexible alignment
    let sumSq = 0, cnt = 0;
    segments.forEach((seg, si) => {
        const { R, t } = transforms[si];
        seg.residues.forEach(r => {
            const p = mapBt.get(r), q = mapA.get(r);
            if (!p || !q) return;
            const nx = R[0][0]*p[0]+R[0][1]*p[1]+R[0][2]*p[2]+t[0];
            const ny = R[1][0]*p[0]+R[1][1]*p[1]+R[1][2]*p[2]+t[1];
            const nz = R[2][0]*p[0]+R[2][1]*p[1]+R[2][2]*p[2]+t[2];
            sumSq += (nx-q[0])**2 + (ny-q[1])**2 + (nz-q[2])**2;
            cnt++;
        });
    });
    const rmsdh = cnt > 0 ? Math.sqrt(sumSq/cnt) : null;

    return { pdb: lines.join('\n'), segments, transforms, rmsdh };
}
