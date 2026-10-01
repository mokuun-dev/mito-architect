import type { Config } from 'tailwindcss';

export default {
  content: ['./index.html', './src/**/*.{ts,tsx}'],
  darkMode: 'class',
  theme: {
    extend: {
      fontFamily: {
        sans: ['Inter', 'IBM Plex Sans', 'ui-sans-serif', 'system-ui', 'sans-serif'],
        mono: ['JetBrains Mono', 'Fira Code', 'ui-monospace', 'SFMono-Regular', 'monospace']
      },
      colors: {
        shell: 'rgb(var(--color-shell) / <alpha-value>)',
        panel: 'rgb(var(--color-panel) / <alpha-value>)',
        panel2: 'rgb(var(--color-panel-2) / <alpha-value>)',
        line: 'rgb(var(--color-line) / <alpha-value>)',
        text: 'rgb(var(--color-text) / <alpha-value>)',
        muted: 'rgb(var(--color-muted) / <alpha-value>)',
        aqua: 'rgb(var(--color-aqua) / <alpha-value>)',
        sky: 'rgb(var(--color-sky) / <alpha-value>)',
        magenta: 'rgb(var(--color-magenta) / <alpha-value>)',
        amber: 'rgb(var(--color-amber) / <alpha-value>)',
        coral: 'rgb(var(--color-coral) / <alpha-value>)',
        leaf: 'rgb(var(--color-leaf) / <alpha-value>)'
      },
      boxShadow: {
        tool: 'var(--shadow-tool)',
        focus: 'var(--shadow-focus)'
      }
    }
  },
  plugins: []
} satisfies Config;
