export interface PairedStrainFiles {
  id: string;
  label: string;
  bamFile: File;
  indexFile: File;
}

export interface MatchedWorkspaceFiles {
  fastaFile?: File;
  faiFile?: File;
  gffFile?: File;
  pairedStrains: PairedStrainFiles[];
  unpairedBams: File[];
  unpairedIndexes: File[];
  unrecognized: File[];
}

let strainSeqCounter = 0;

function stripBamSuffix(name: string): string {
  return name.replace(/\.bam$/i, '');
}

function stripIndexSuffix(name: string): string {
  return name.replace(/\.bam\.(bai|csi)$/i, '').replace(/\.(bai|csi)$/i, '');
}

export function matchDroppedFiles(files: File[]): MatchedWorkspaceFiles {
  let fastaFile: File | undefined;
  let faiFile: File | undefined;
  let gffFile: File | undefined;

  const bamFiles: File[] = [];
  const indexFiles: File[] = [];
  const unrecognized: File[] = [];

  for (const f of files) {
    const lower = f.name.toLowerCase();
    if (/\.(fa|fasta|fna|fsa)(\.gz)?$/.test(lower)) {
      fastaFile = f;
    } else if (/\.fai$/.test(lower)) {
      faiFile = f;
    } else if (/\.(gff3?|gtf)(\.gz)?$/.test(lower)) {
      gffFile = f;
    } else if (/\.bam$/.test(lower)) {
      bamFiles.push(f);
    } else if (/\.(bai|csi)$/.test(lower)) {
      indexFiles.push(f);
    } else {
      unrecognized.push(f);
    }
  }

  const usedIndexes = new Set<File>();
  const pairedStrains: PairedStrainFiles[] = [];
  const unpairedBams: File[] = [];

  for (const bam of bamFiles) {
    const bamLower = bam.name.toLowerCase();
    const bamStem = stripBamSuffix(bamLower);

    // 1순위: foo.bam.bai / foo.bam.csi, 2순위: foo.bai / foo.csi
    const exactMatch =
      indexFiles.find(
        (idx) =>
          !usedIndexes.has(idx) &&
          (idx.name.toLowerCase() === `${bamLower}.bai` ||
            idx.name.toLowerCase() === `${bamLower}.csi`),
      ) ??
      indexFiles.find(
        (idx) => !usedIndexes.has(idx) && stripIndexSuffix(idx.name.toLowerCase()) === bamStem,
      );

    if (exactMatch) {
      usedIndexes.add(exactMatch);
      pairedStrains.push({
        id: `strain_${++strainSeqCounter}_${bam.name}`,
        label: stripBamSuffix(bam.name),
        bamFile: bam,
        indexFile: exactMatch,
      });
    } else {
      unpairedBams.push(bam);
    }
  }

  const unpairedIndexes = indexFiles.filter((idx) => !usedIndexes.has(idx));

  return {
    fastaFile,
    faiFile,
    gffFile,
    pairedStrains,
    unpairedBams,
    unpairedIndexes,
    unrecognized,
  };
}

export function formatBytes(bytes: number): string {
  if (bytes < 1024) return `${bytes} B`;
  if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(1)} KB`;
  if (bytes < 1024 * 1024 * 1024) return `${(bytes / (1024 * 1024)).toFixed(1)} MB`;
  return `${(bytes / (1024 * 1024 * 1024)).toFixed(2)} GB`;
}

