import type { Config } from 'tailwindcss';

const token = (name: string) => `rgb(var(--color-${name}) / <alpha-value>)`;

export default {
  content: ['./index.html', './src/**/*.{ts,tsx}'],
  darkMode: 'class',
  theme: {
    extend: {
      fontFamily: {
        sans: ['Inter Variable', 'ui-sans-serif', 'system-ui', 'sans-serif'],
        mono: ['JetBrains Mono Variable', 'ui-monospace', 'SFMono-Regular', 'monospace']
      },
      colors: {
        shell: token('shell'),
        panel: token('panel'),
        panel2: token('panel-2'),
        line: { DEFAULT: token('line'), strong: token('line-strong') },
        text: token('text'),
        muted: token('muted'),
        subtle: token('subtle'),
        accent: { DEFAULT: token('accent'), strong: token('accent-strong') },
        onaccent: token('on-accent'),
        gene: { coding: token('gene-coding'), trna: token('gene-trna'), rrna: token('gene-rrna') },
        variant: token('variant'),
        sv: token('sv'),
        coverage: token('coverage'),
        sky: token('sky'),
        magenta: token('magenta'),
        amber: token('amber'),
        coral: token('coral'),
        leaf: token('leaf')
      },
      boxShadow: {
        card: 'var(--shadow-card)',
        pop: 'var(--shadow-pop)',
        tool: 'var(--shadow-tool)',
        focus: 'var(--shadow-focus)'
      }
    }
  },
  plugins: []
} satisfies Config;
