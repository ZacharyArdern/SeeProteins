const puppeteer = require('puppeteer');

(async () => {
    const browser = await puppeteer.launch({
        headless: 'new',
        args: ['--no-sandbox', '--disable-web-security'],
    });
    const page = await browser.newPage();

    const logs = [];
    page.on('console', m => logs.push(`[${m.type()}] ${m.text()}`));
    page.on('pageerror', e => logs.push(`[pageerror] ${e.message}\n${e.stack}`));
    page.on('requestfailed', r => logs.push(`[fail] ${r.url()} — ${r.failure()?.errorText}`));
    page.on('response', r => { const s = r.status(); if (s >= 400 || r.url().includes('.wasm') || r.url().includes('.onnx')) logs.push(`[resp ${s}] ${r.url()}`); });

    // Navigate to a small dedicated ORT test page served by Vite (same origin = no CORS)
    await page.goto('http://localhost:5174/ort-test.html', { waitUntil: 'domcontentloaded' });

    // Wait up to 30s for the test to complete
    await page.waitForFunction(
        () => window._ortTestDone === true,
        { timeout: 30000 }
    ).catch(() => logs.push('[timeout] test did not complete within 30s'));

    await browser.close();

    for (const l of logs) console.log(l);
})();
