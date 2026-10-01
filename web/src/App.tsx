import { CircleDot, Dna, Moon, Sun } from 'lucide-react';
import { useEffect, useState } from 'react';
import { Link, Route, Routes } from 'react-router-dom';
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
    document.querySelector('meta[name="theme-color"]')?.setAttribute('content', light ? '#f4f6f8' : '#111820');
  }, [light]);

  return (
    <div className="min-h-screen bg-shell text-text">
      <a
        href="#main-content"
        className="fixed left-4 top-3 z-50 -translate-y-24 rounded-md bg-aqua px-4 py-2 text-sm font-semibold text-shell shadow-focus focus:translate-y-0"
      >
        Skip to analysis content
      </a>
      <header className="sticky top-0 z-20 border-b border-line bg-panel">
        <div className="mx-auto flex h-14 max-w-[1800px] items-center justify-between gap-3 px-4 sm:px-5">
          <Link to="/" className="flex min-w-0 items-center gap-3 rounded-md">
            <span className="grid h-8 w-8 shrink-0 place-items-center rounded-md bg-aqua/10">
              <Dna className="h-5 w-5 text-aqua" aria-hidden />
            </span>
            <div className="min-w-0">
              <div className="truncate text-sm font-semibold tracking-tight">Mito Architect</div>
            </div>
          </Link>
          <div className="flex items-center gap-2">
            <span className="hidden items-center gap-2 text-xs text-muted md:inline-flex">
              <CircleDot className="h-3 w-3 text-leaf" aria-hidden />
              Для исследований
            </span>
            <button
              type="button"
              onClick={() => setLight((value) => !value)}
              className="inline-flex h-10 w-10 items-center justify-center rounded-md border border-line bg-panel text-muted hover:border-aqua hover:text-text"
              title={light ? 'Use dark theme' : 'Use light theme'}
              aria-label={light ? 'Use dark theme' : 'Use light theme'}
              aria-pressed={light}
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
