import { defineConfig } from 'vite'
import { readFileSync } from 'fs'
import { resolve, join } from 'path'

export default defineConfig({
  plugins: [
    {
      name: 'serve-wasm-binary',
      configureServer(server) {
        server.middlewares.use((req, res, next) => {
          if (req.url?.includes('.wasm') && !req.url.includes('?')) {
            const filePath = resolve(join(process.cwd(), req.url.split('?')[0]))
            try {
              const data = readFileSync(filePath)
              res.setHeader('Content-Type', 'application/wasm')
              res.setHeader('Cache-Control', 'no-cache')
              res.end(data)
              return
            } catch {}
          }
          next()
        })
      }
    }
  ]
})
