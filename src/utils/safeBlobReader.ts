export const SAFE_IO_SLICE_BYTES = 512 * 1024; // 512 KB
export const SAFE_IO_SLICE_TIMEOUT_MS = 8_000; // 8s per 512 KB slice

export type MainThreadBlobRelay = (start: number, end: number) => Promise<Uint8Array>;

export interface SafeBlobReadOptions {
  signal?: AbortSignal;
  timeoutMs?: number;
  onProgress?: (loadedBytes: number, totalBytes: number) => void;
  mainThreadRelay?: MainThreadBlobRelay;
}

// Worker/앱 스레드 내 모든 로컬 File 디스크 I/O(BAM·BAI·FASTA)를 동시성 1로 직렬화하여
// Windows 파일 락 경합과 DLP/보안 필터 병목(중간 멈춤)을 차단합니다.
let globalIoQueue: Promise<void> = Promise.resolve();

export function runGlobalSerializedIo<T>(fn: () => Promise<T>): Promise<T> {
  const next = globalIoQueue.then(fn, fn);
  globalIoQueue = next.then(
    () => undefined,
    () => undefined,
  );
  return next;
}

function isAbortError(err: unknown): boolean {
  if (!err || typeof err !== 'object') return false;
  return (err as { name?: string }).name === 'AbortError';
}

function throwIfAborted(signal?: AbortSignal): void {
  if (!signal) return;
  if (typeof signal.throwIfAborted === 'function') {
    signal.throwIfAborted();
  } else if (signal.aborted) {
    const e = new Error('작업이 취소되었습니다.');
    e.name = 'AbortError';
    throw e;
  }
}

export async function raceWithTimeout<T>(
  promise: Promise<T>,
  timeoutMs: number,
  timeoutMessage: string,
  signal?: AbortSignal,
): Promise<T> {
  throwIfAborted(signal);
  let timerId: ReturnType<typeof setTimeout> | undefined;
  let onAbort: (() => void) | undefined;

  const guardPromise = new Promise<never>((_, reject) => {
    timerId = setTimeout(() => {
      reject(new Error(timeoutMessage));
    }, timeoutMs);

    if (signal) {
      onAbort = () => {
        const err = new Error('작업이 취소되었습니다.');
        err.name = 'AbortError';
        reject(err);
      };
      signal.addEventListener('abort', onAbort, { once: true });
    }
  });

  try {
    return await Promise.race([promise, guardPromise]);
  } finally {
    if (timerId !== undefined) clearTimeout(timerId);
    if (signal && onAbort) {
      signal.removeEventListener('abort', onAbort);
    }
  }
}

function readViaFileReader(slice: Blob): Promise<Uint8Array> {
  return new Promise<Uint8Array>((resolve, reject) => {
    if (typeof FileReader === 'undefined') {
      reject(new Error('FileReader unavailable'));
      return;
    }
    const reader = new FileReader();
    reader.onload = () => {
      if (reader.result instanceof ArrayBuffer) {
        resolve(new Uint8Array(reader.result));
      } else {
        reject(new Error('FileReader returned non-ArrayBuffer result'));
      }
    };
    reader.onerror = () => {
      reject(reader.error ?? new Error('FileReader error'));
    };
    reader.onabort = () => {
      const err = new Error('작업이 취소되었습니다.');
      err.name = 'AbortError';
      reject(err);
    };
    reader.readAsArrayBuffer(slice);
  });
}

function formatReadableIoError(fileName: string, err: unknown): Error {
  if (isAbortError(err)) {
    return err instanceof Error ? err : new Error('작업이 취소되었습니다.');
  }
  const raw = err instanceof Error ? err.message : String(err);
  const name = err && typeof err === 'object' ? (err as { name?: string }).name : '';
  if (name === 'NotReadableError' || /could not be read|permission/i.test(raw)) {
    return new Error(
      `'${fileName}' 파일을 읽을 수 없습니다(OS 잠금·권한 변경 또는 보안 DLP 간섭). 다른 프로그램에서 파일을 사용 중인지 확인하거나 로컬 SSD로 복사한 뒤 다시 선택해주세요.`,
    );
  }
  return new Error(
    `'${fileName}' 파일 읽기 응답이 지연되고 있습니다 (${raw}). 네트워크/OneDrive 폴더이거나 보안(DLP) 검사가 걸린 경우 로컬 SSD로 복사해 다시 시도해주세요.`,
  );
}

/**
 * 단일 File 슬라이스를 전역 직렬화 큐 안에서 안전하게 읽습니다.
 * 1차: slice.arrayBuffer() -> 2차: FileReader 폴백 -> 3차: 메인 스레드 Blob 중계 폴백
 */
export async function readBlobSliceSafe(
  file: File,
  start: number,
  end: number,
  opts: SafeBlobReadOptions = {},
): Promise<Uint8Array> {
  const signal = opts.signal;
  throwIfAborted(signal);

  const clampedStart = Math.max(0, Math.min(file.size, Math.floor(start)));
  const clampedEnd = Math.max(clampedStart, Math.min(file.size, Math.floor(end)));
  if (clampedEnd <= clampedStart) return new Uint8Array(0);

  const timeoutMs = opts.timeoutMs ?? SAFE_IO_SLICE_TIMEOUT_MS;
  const secLabel = `${Math.round(timeoutMs / 1000)}초 초과`;

  return runGlobalSerializedIo(async () => {
    throwIfAborted(signal);
    const slice = file.slice(clampedStart, clampedEnd);

    const attempts: Array<{ label: string; run: () => Promise<Uint8Array> }> = [
      {
        label: 'slice.arrayBuffer()',
        run: () => slice.arrayBuffer().then((ab) => new Uint8Array(ab)),
      },
    ];
    if (typeof FileReader !== 'undefined') {
      attempts.push({
        label: 'FileReader',
        run: () => readViaFileReader(slice),
      });
    }
    if (opts.mainThreadRelay) {
      const relay = opts.mainThreadRelay;
      attempts.push({
        label: 'mainThreadRelay',
        run: () => relay(clampedStart, clampedEnd),
      });
    }

    let lastErr: unknown;
    for (const attempt of attempts) {
      throwIfAborted(signal);
      try {
        const bytes = await raceWithTimeout(
          attempt.run(),
          timeoutMs,
          `${secLabel}: ${attempt.label}`,
          signal,
        );
        throwIfAborted(signal);
        return bytes;
      } catch (err) {
        if (isAbortError(err)) throw err;
        lastErr = err;
      }
    }

    throw formatReadableIoError(file.name, lastErr);
  });
}

/**
 * 큰 구간(.bai 전체 또는 대형 BAM 청크)을 512 KB 단위로 분할·스트리밍 읽기하여
 * 진행률 보고, 타임아웃, 취소 응답성을 보장합니다.
 */
export async function readBlobRangeInChunks(
  file: File,
  start: number,
  end: number,
  opts: SafeBlobReadOptions = {},
): Promise<Uint8Array> {
  const clampedStart = Math.max(0, Math.min(file.size, Math.floor(start)));
  const clampedEnd = Math.max(clampedStart, Math.min(file.size, Math.floor(end)));
  const totalLen = clampedEnd - clampedStart;
  if (totalLen <= 0) return new Uint8Array(0);

  if (totalLen <= SAFE_IO_SLICE_BYTES) {
    const single = await readBlobSliceSafe(file, clampedStart, clampedEnd, opts);
    opts.onProgress?.(single.byteLength, totalLen);
    return single;
  }

  const out = new Uint8Array(totalLen);
  let offset = 0;
  opts.onProgress?.(0, totalLen);

  while (offset < totalLen) {
    throwIfAborted(opts.signal);
    const s = clampedStart + offset;
    const e = Math.min(clampedEnd, s + SAFE_IO_SLICE_BYTES);
    const chunk = await readBlobSliceSafe(file, s, e, opts);
    if (chunk.byteLength === 0) break;
    const copyLen = Math.min(chunk.byteLength, totalLen - offset);
    out.set(copyLen === chunk.byteLength ? chunk : chunk.subarray(0, copyLen), offset);
    offset += copyLen;
    opts.onProgress?.(offset, totalLen);
  }

  return offset === totalLen ? out : out.subarray(0, offset);
}
