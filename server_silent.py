"""Silent background static HTTP server for YeastWGS Studio (Windows pythonw compatible)."""

from __future__ import annotations

import functools
import http.server
import os
import shutil
import signal
import socket
import subprocess
import sys
import threading
import time
import webbrowser

HOST = "127.0.0.1"
PORT = 4173
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DIST_DIR = os.path.join(BASE_DIR, "dist")
PID_FILE = os.path.join(BASE_DIR, ".viewer_server.pid")


def find_chrome_executable() -> str | None:
    """Locate Google Chrome executable across Windows Registry and standard OS paths."""
    if sys.platform == "win32":
        try:
            import winreg

            for root in (winreg.HKEY_CURRENT_USER, winreg.HKEY_LOCAL_MACHINE):
                try:
                    with winreg.OpenKey(
                        root,
                        r"SOFTWARE\Microsoft\Windows\CurrentVersion\App Paths\chrome.exe",
                    ) as key:
                        val, _ = winreg.QueryValueEx(key, None)
                        if isinstance(val, str) and os.path.isfile(val):
                            return val
                except OSError:
                    pass
        except ImportError:
            pass

        win_candidates = [
            os.path.join(
                os.environ.get("PROGRAMFILES", r"C:\Program Files"),
                "Google",
                "Chrome",
                "Application",
                "chrome.exe",
            ),
            os.path.join(
                os.environ.get("PROGRAMFILES(X86)", r"C:\Program Files (x86)"),
                "Google",
                "Chrome",
                "Application",
                "chrome.exe",
            ),
            os.path.join(
                os.environ.get("LOCALAPPDATA", ""),
                "Google",
                "Chrome",
                "Application",
                "chrome.exe",
            ),
        ]
        for cand in win_candidates:
            if cand and os.path.isfile(cand):
                return cand
    elif sys.platform == "darwin":
        mac_chrome = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
        if os.path.isfile(mac_chrome):
            return mac_chrome

    for name in ("chrome", "google-chrome", "chromium", "chromium-browser"):
        found = shutil.which(name)
        if found:
            return found
    return None


def open_in_chrome(url: str) -> None:
    """Open URL preferentially in Google Chrome without showing a console window."""
    chrome_exe = find_chrome_executable()
    if chrome_exe:
        try:
            kwargs: dict[str, object] = {
                "stdin": subprocess.DEVNULL,
                "stdout": subprocess.DEVNULL,
                "stderr": subprocess.DEVNULL,
            }
            if sys.platform == "win32":
                kwargs["creationflags"] = (
                    getattr(subprocess, "CREATE_NO_WINDOW", 0x08000000)
                    | getattr(subprocess, "DETACHED_PROCESS", 0x00000008)
                )
            subprocess.Popen([chrome_exe, url], **kwargs)
            return
        except OSError:
            pass
    webbrowser.open(url)


class SilentStaticHandler(http.server.SimpleHTTPRequestHandler):
    """Static file handler that never writes to None stdout/stderr under pythonw.exe."""

    extensions_map = {
        **http.server.SimpleHTTPRequestHandler.extensions_map,
        ".js": "application/javascript; charset=utf-8",
        ".mjs": "application/javascript; charset=utf-8",
        ".css": "text/css; charset=utf-8",
        ".json": "application/json; charset=utf-8",
        ".svg": "image/svg+xml",
        ".wasm": "application/wasm",
        ".gz": "application/gzip",
    }

    def end_headers(self) -> None:
        self.send_header("Cache-Control", "no-store, no-cache, must-revalidate")
        self.send_header("Pragma", "no-cache")
        super().end_headers()

    def log_message(self, format: str, *args: object) -> None:
        # pythonw.exe sets sys.stderr to None; suppress request logs silently.
        return

    def handle_error(self, request: object, client_address: object) -> None:
        # Suppress BrokenPipeError / ConnectionResetError tracebacks cleanly.
        return


class QuietThreadingHTTPServer(http.server.ThreadingHTTPServer):
    daemon_threads = True
    allow_reuse_address = True

    def handle_error(self, request: object, client_address: object) -> None:
        return


def is_server_running(host: str = HOST, port: int = PORT) -> bool:
    try:
        with socket.create_connection((host, port), timeout=0.5):
            return True
    except OSError:
        return False


def stop_existing_server() -> bool:
    if not os.path.exists(PID_FILE):
        return False
    try:
        with open(PID_FILE, "r", encoding="utf-8") as fh:
            pid = int(fh.read().strip())
        os.kill(pid, signal.SIGTERM)
        for _ in range(20):
            if not is_server_running():
                break
            time.sleep(0.1)
    except (OSError, ValueError):
        pass
    try:
        if os.path.exists(PID_FILE):
            os.remove(PID_FILE)
    except OSError:
        pass
    return not is_server_running()


def main() -> int:
    if "--stop" in sys.argv:
        stop_existing_server()
        return 0

    if not os.path.isfile(os.path.join(DIST_DIR, "index.html")):
        return 1

    url = f"http://{HOST}:{PORT}"
    if is_server_running():
        if "--no-browser" not in sys.argv:
            open_in_chrome(url)
        return 0

    handler_cls = functools.partial(SilentStaticHandler, directory=DIST_DIR)
    httpd = QuietThreadingHTTPServer((HOST, PORT), handler_cls)

    try:
        with open(PID_FILE, "w", encoding="utf-8") as fh:
            fh.write(str(os.getpid()))
    except OSError:
        pass

    if "--no-browser" not in sys.argv:
        t = threading.Timer(0.1, open_in_chrome, args=(url,))
        t.daemon = True
        t.start()

    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        httpd.server_close()
        try:
            if os.path.exists(PID_FILE):
                with open(PID_FILE, "r", encoding="utf-8") as fh:
                    if fh.read().strip() == str(os.getpid()):
                        os.remove(PID_FILE)
        except OSError:
            pass
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
