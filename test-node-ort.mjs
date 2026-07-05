// Test onnxruntime-web WASM backend from Node.js
import * as ort from 'onnxruntime-web';
import { readFile } from 'fs/promises';
import { pathToFileURL } from 'url';

ort.env.wasm.numThreads = 1;

const modelData = await readFile('public/thermompnn.onnx');

console.log('Creating session...');
try {
    const sess = await ort.InferenceSession.create(modelData, {
        executionProviders: ['wasm'],
    });
    console.log('Session OK, inputs:', sess.inputNames.join(', '));
} catch(e) {
    console.error('Session create ERROR:', e.message);
    console.error(e.stack);
}
