import createBICExact from './bicexact.js';

const WASM_URL = '/bicexact.wasm';

function detectChain(pdbText) {
    for (const line of pdbText.split('\n'))
        if (line.startsWith('ATOM')) return line.substring(21, 22).trim() || 'A';
    return 'A';
}

function parse(output) {
    let rmsd = null, rmsdh = null, hingeCount = 0, hingeResidues = [];
    for (const line of output.split('\n')) {
        if (line.startsWith('RMSD:'))
            rmsd = parseFloat(line.split(':')[1]);
        else if (line.startsWith('RMSDhk is:'))
            rmsdh = parseFloat(line.split(':')[1]);
        else if (line.startsWith('The number of hinges is:'))
            hingeCount = parseInt(line.split(':')[1]);
        else if (line.startsWith('Hinge index is:')) {
            const raw = line.replace('Hinge index is:', '').trim();
            hingeResidues = raw
                ? raw.split(':').map(s => parseInt(s.trim())).filter(n => !isNaN(n))
                : [];
        }
    }
    return { rmsd, rmsdh, hingeCount, hingeResidues };
}

export async function runBICExact(pdb1, pdb2, criterion = 'bic', mode = 'lh') {
    const chain1 = detectChain(pdb1);
    const chain2 = detectChain(pdb2);
    let output = '';
    const instance = await createBICExact({
        locateFile: () => WASM_URL,
        print:    s => { output += s + '\n'; },
        printErr: () => {}
    });
    instance.FS.writeFile('/pdb1.pdb', pdb1);
    instance.FS.writeFile('/pdb2.pdb', pdb2);
    instance.callMain(['/pdb1.pdb', '/pdb2.pdb', chain1, chain2, criterion, mode]);
    return parse(output);
}
