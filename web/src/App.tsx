import { Dna, FlaskConical, Moon, Sun } from 'lucide-react';
import { useEffect, useState } from 'react';
import { Link, NavLink, Route, Routes } from 'react-router-dom';
import Home from './routes/Home';
import ResultPage from './routes/ResultPage';
import UploadPage from './routes/UploadPage';

export default function App() {
  const [light, setLight] = useState(() => document.documentElement.classList.contains('light'));

  useEffect(() => {
    document.documentElement.classList.toggle('light', light);
    try {
      window.localStorage.setItem('mito-theme', light ? 'light' : 'dark');
    } catch {
      // Theme selection still applies for the current session when storage is unavailable.
    }
    document.querySelector('meta[name="theme-color"]')?.setAttribute('content', light ? '#ffffff' : '#141b24');
  }, [light]);

  return (
    <div className="min-h-screen bg-shell text-text">
      <a
        href="#main-content"
        className="fixed left-4 top-3 z-50 -translate-y-24 rounded-lg bg-accent px-4 py-2 text-sm font-semibold text-onaccent shadow-focus focus:translate-y-0"
      >
        Skip to analysis content
      </a>
      <header className="sticky top-0 z-30 border-b border-line bg-panel/90 backdrop-blur supports-[backdrop-filter]:bg-panel/75">
        <div className="mx-auto flex h-14 max-w-[1600px] items-center gap-6 px-4 sm:px-6">
          <Link to="/" className="flex min-w-0 shrink-0 items-center gap-2.5 rounded-lg">
            <span className="grid h-8 w-8 place-items-center rounded-lg bg-accent text-onaccent">
              <Dna className="h-[18px] w-[18px]" aria-hidden />
            </span>
            <span className="text-[15px] font-semibold tracking-tight">Mito-Architect</span>
          </Link>
          <nav aria-label="Primary" className="hidden items-center gap-1 sm:flex">
            <HeaderLink to="/" end>New analysis</HeaderLink>
            <HeaderLink to="/result/demo">Demo result</HeaderLink>
          </nav>
          <div className="ml-auto flex items-center gap-3">
            <span className="hidden items-center gap-1.5 text-xs font-medium text-muted md:inline-flex" title="Results are for research use and are not clinical diagnoses.">
              <FlaskConical className="h-3.5 w-3.5" aria-hidden />
              Research use only
            </span>
            <button
              type="button"
              onClick={() => setLight((value) => !value)}
              className="btn-icon"
              title={light ? 'Use dark theme' : 'Use light theme'}
              aria-label={light ? 'Use dark theme' : 'Use light theme'}
              aria-pressed={!light}
            >
              {light ? <Moon className="h-4 w-4" aria-hidden /> : <Sun className="h-4 w-4" aria-hidden />}
            </button>
          </div>
        </div>
      </header>
      <Routes>
        <Route path="/" element={<Home />} />
        <Route path="/upload" element={<UploadPage />} />
        <Route path="/result/:jobId" element={<ResultPage />} />
      </Routes>
    </div>
  );
}

function HeaderLink({ to, end, children }: { to: string; end?: boolean; children: string }) {
  return (
    <NavLink
      to={to}
      end={end}
      className={({ isActive }) =>
        'rounded-lg px-3 py-1.5 text-sm font-medium ' + (isActive ? 'bg-panel2 text-text' : 'text-muted hover:text-text')
      }
    >
      {children}
    </NavLink>
  );
}
