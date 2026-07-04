import createTmalign from './node_modules/tmalign-wasm/tmalign-wasm.js'

const WASM_URL = '/tmalign-wasm.wasm';

export function tmalign(pdb1, pdb2, alignment = null) {
    return new Promise((resolve, reject) => {
        let buffer = "";
        createTmalign({
            locateFile: () => WASM_URL,
            print: (msg) => buffer += msg + "\n"
        }).then((instance) => {
            const cmd = ['/pdb1.pdb', '/pdb2.pdb', '-m', '/matrix.txt'];
            instance.FS.writeFile('/pdb1.pdb', pdb1);
            instance.FS.writeFile('/pdb2.pdb', pdb2);
            if (alignment) {
                cmd.push('-I', '/aln.fa');
                instance.FS.writeFile('/aln.fa', alignment);
            }
            const err = instance.callMain(cmd);
            if (err == 0) {
                const matrix = instance.FS.readFile('/matrix.txt', { encoding: 'utf8' });
                resolve({ output: buffer, matrix });
            } else {
                reject(new Error(`TM-align exited with code ${err}`));
            }
        });
    });
}

export function parse(output) {
    const lines = output.split('\n');
    let chain1, chain2, tmScore, rmsd;
    for (let i = 0; i < lines.length; i++) {
        const line = lines[i];
        if (line.startsWith('Aligned length='))
            rmsd = parseFloat(line.match(/RMSD=\s*(\d*\.\d*),/)[1]);
        if (line.startsWith('Name of Chain_1:'))
            chain1 = line.split(' ')[3].trim();
        if (line.startsWith('Name of Chain_2:'))
            chain2 = line.split(' ')[3].trim();
        if (line.startsWith('TM-score=') && line.includes('Chain_1'))
            tmScore = parseFloat(line.split(' ')[1]);
        if (line.startsWith('(":" denotes')) {
            const q = lines[i + 1].trim();
            const t = lines[i + 3].trim();
            let qAln = '', tAln = '', op = '', cigar = '', len = 0;
            let qPos = 0, tPos = 0, qStart = 0, tStart = 0, firstM = true, lastMatch = 0;
            for (let j = 0; j < q.length; j++) {
                const qc = q[j], tc = t[j];
                const gap = qc === '-' || tc === '-';
                const match = !gap;
                if (match) {
                    if (op !== 'M' && len) { cigar += len + op; len = 0; }
                    if (firstM) { qStart = qPos; tStart = tPos; cigar = ''; firstM = false; }
                    op = 'M'; len++;
                    qAln += qc; tAln += tc;
                    qPos++; tPos++;
                    lastMatch = j;
                } else if (qc === '-') {
                    if (op !== 'D' && len) { cigar += len + op; len = 0; }
                    op = 'D'; len++;
                    qAln += '-'; tAln += tc;
                    tPos++;
                } else {
                    if (op !== 'I' && len) { cigar += len + op; len = 0; }
                    op = 'I'; len++;
                    qAln += qc; tAln += '-';
                    qPos++;
                }
            }
            if (len) cigar += len + op;
            const lastM = cigar.lastIndexOf('M');
            return {
                query: chain1, target: chain2,
                tmScore, rmsd,
                cigar: cigar.substring(0, lastM + 1),
                qAln: qAln.substring(0, lastMatch + 1),
                tAln: tAln.substring(0, lastMatch + 1),
            };
        }
    }
    return null;
}

export function parseMatrix(matrix) {
    const rows = matrix.split('\n').slice(2, 5).map(line =>
        line.split(/\s+/).slice(1).map(Number)
    );
    return {
        t: [rows[0][0], rows[1][0], rows[2][0]],
        u: [rows[0].slice(1), rows[1].slice(1), rows[2].slice(1)],
    };
}
