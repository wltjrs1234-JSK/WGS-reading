import { describe, expect, it } from 'vitest';
import {
  buildAnalysisFromCalls,
  buildThreeWayAnalysis,
  computeConsensusForRegion,
  computeRegionMappingStats,
  sliceSplicedAnalysis,
} from '../consensus';
import { createYeastDemoDataset } from '../demoData';
import { IndexedFastaSource } from '../fastaIndex';
import {
  compareExtractedProteins,
  evaluateStrainProtein,
  reverseComplement,
  translateDna,
} from '../geneticCode';
import { AnnotationIndex, parseGff3Text } from '../gff3Parser';
import { parseCigar } from '../cigar';
import { packGeneFeatures } from '../geneLayout';
import { deriveAlleleTicksFromCalls, packReadsWithHeap } from '../readPacking';
import {
  isMitochondrialChrom,
  normalizeYeastChromName,
  resolveMatchingRefName,
  sortYeastContigs,
  YEAST_CHROMOSOMES,
} from '../yeastAliases';
import { createBuiltInS288CAnnotationIndex } from '../s288cBuiltIn';
import { matchDroppedFiles } from '../../utils/fileMatcher';
import { DEFAULT_CONSENSUS_OPTIONS, type AlignedRead } from '../types';

describe('yeastAliases', () => {
  it('normalizes roman, arabic, and RefSeq chromosome names', () => {
    expect(normalizeYeastChromName('V')).toBe('chrV');
    expect(normalizeYeastChromName('chr5')).toBe('chrV');
    expect(normalizeYeastChromName('NC_001137.3')).toBe('chrV');
    expect(normalizeYeastChromName('MT')).toBe('chrM');
    expect(normalizeYeastChromName('chrmt')).toBe('chrM');
    expect(normalizeYeastChromName('NC_001224.1')).toBe('chrM');
  });

  it('detects mitochondrial chromosomes accurately', () => {
    expect(isMitochondrialChrom('chrM')).toBe(true);
    expect(isMitochondrialChrom('Mito')).toBe(true);
    expect(isMitochondrialChrom('NC_001224.1')).toBe(true);
    expect(isMitochondrialChrom('chrV')).toBe(false);
  });

  it('matches Reference names against BAM contigs using aliases and overrides', () => {
    const bamContigs = ['I', 'V', 'MT', 'custom_contig'];
    expect(resolveMatchingRefName('chrV', bamContigs)).toBe('V');
    expect(resolveMatchingRefName('chrM', bamContigs)).toBe('MT');
    expect(
      resolveMatchingRefName('chrV', bamContigs, { chrV: 'custom_contig' }),
    ).toBe('custom_contig');
  });
});

describe('geneticCode (Nuclear Table 1 vs Yeast Mitochondrial Table 3)', () => {
  it('computes reverse complement preserving IUPAC and gap characters', () => {
    expect(reverseComplement('ATGRY-')).toBe('-RYCAT');
    expect(reverseComplement('ATGAAACCC')).toBe('GGGTTTCAT');
  });

  it('translates chrM using Yeast Mitochondrial Table 3 (TGA->W, ATA->M, CTT->T)', () => {
    const seq = 'ATGATACTTTGATAA'; // M, I/M, L/T, */W, *
    expect(translateDna(seq, 1, 0)).toBe('MIL**');
    expect(translateDna(seq, 3, 0)).toBe('MMTW*');
  });

  it('classifies synonymous, missense, nonsense, and frameshift protein changes', () => {
    const cmp = compareExtractedProteins(
      'ATGGCTAAATAA', // M A K *
      'ATGGTTAAATAA', // M V K * (Missense)
      'ATGGCTAAAATAA', // +1bp (Frameshift)
      1,
      0,
    );
    expect(cmp.s1Status).toBe('Missense');
    expect(cmp.s1Changes).toEqual(['A2V']);
    expect(cmp.s2Status).toBe('Frameshift');
  });
});

describe('CIGAR-aware Consensus & 3-Way Analysis', () => {
  const ref = 'ATGCGATCGATCGATC'; // 16 bp

  function makeRead(id: string, start: number, end: number, seq: string, cigar: string): AlignedRead {
    return {
      id,
      name: id,
      start,
      end,
      strand: '+',
      seq,
      cigar,
      mapq: 60,
      flags: 0,
    };
  }

  it('handles SNP, Insertion, Deletion, and aligns all 3 rows identically', () => {
    // Strain 1: SNP at pos 2 (G->A) + Insertion 'TT' after pos 5 -> 6M2I10M
    const s1Seq = 'ATACGATTTCGATCGATC'; // len 18
    const s1Reads = Array.from({ length: 8 }, (_, i) =>
      makeRead(`s1_${i}`, 0, 16, s1Seq, '6M2I10M'),
    );

    // Strain 2: 3bp Deletion at pos 8..11 (GAT) -> 8M3D5M
    const s2Seq = 'ATGCGATCCGATC'; // len 13
    const s2Reads = Array.from({ length: 8 }, (_, i) =>
      makeRead(`s2_${i}`, 0, 16, s2Seq, '8M3D5M'),
    );

    const analysis = buildThreeWayAnalysis(
      'chrI',
      0,
      16,
      ref,
      s1Reads,
      s2Reads,
      DEFAULT_CONSENSUS_OPTIONS,
    );

    // Ungapped outputs reflect each strain's true physical DNA sequence
    expect(analysis.ungappedRef).toBe('ATGCGATCGATCGATC');
    expect(analysis.ungappedS1).toBe(s1Seq);
    expect(analysis.ungappedS2).toBe(s2Seq);

    // Aligned MSA rows must all share exact equal length (16 + 2bp insertion = 18)
    expect(analysis.alignedRef.length).toBe(18);
    expect(analysis.alignedS1.length).toBe(18);
    expect(analysis.alignedS2.length).toBe(18);
    expect(analysis.alignedRef).toBe('ATGCGA--TCGATCGATC');
    expect(analysis.alignedS1).toBe(s1Seq);
    expect(analysis.alignedS2).toBe('ATGCGA--TC---CGATC');

    expect(analysis.variants.some((v) => v.category === 'strain1_only' && v.type === 'SNP')).toBe(true);
    expect(analysis.variants.some((v) => v.category === 'strain1_only' && v.type === 'INS')).toBe(true);
    expect(analysis.variants.some((v) => v.category === 'strain2_only' && v.type === 'DEL')).toBe(true);
  });

  it('calls heterozygous IUPAC bases in diploid mode', () => {
    const reads: AlignedRead[] = [
      ...Array.from({ length: 5 }, (_, i) => makeRead(`a_${i}`, 0, 16, 'ATGCGATCGATCGATC', '16M')),
      ...Array.from({ length: 5 }, (_, i) => makeRead(`b_${i}`, 0, 16, 'ATACGATCGATCGATC', '16M')),
    ];
    const calls = computeConsensusForRegion(ref, reads, 0, 16, {
      ...DEFAULT_CONSENSUS_OPTIONS,
      ploidyMode: 'diploid',
      hetThreshold: 0.25,
    });
    expect(calls[2].call).toBe('R'); // A + G -> R
    expect(calls[2].isHeterozygous).toBe(true);
  });
});

describe('IndexedFastaSource (without .fai and with .fai)', () => {
  it('builds in-memory index and slices multi-line FASTA accurately', async () => {
    const fastaText = '>chrI synthetic\nATGCCCGG\nGTTTAAAC\n>chrM mito\nATATATAT\nGGGG\n';
    const fastaFile = new File([fastaText], 'test.fasta', { type: 'text/plain' });

    const source = await IndexedFastaSource.fromFiles(fastaFile);
    expect(source.getContigs()).toEqual([
      { name: 'chrI', length: 16 },
      { name: 'chrM', length: 12 },
    ]);

    expect(await source.getSequence('chrI', 4, 12)).toBe('CCGGGTTT');
    expect(await source.getSequence('chrM', 6, 12)).toBe('ATGGGG');
  });

  it('uses provided .fai file cleanly', async () => {
    const fastaText = '>chrI\nAAAA\nCCCC\n';
    const faiText = 'chrI\t8\t6\t4\t5\n';
    const fastaFile = new File([fastaText], 'test.fa');
    const faiFile = new File([faiText], 'test.fa.fai');

    const source = await IndexedFastaSource.fromFiles(fastaFile, faiFile);
    expect(await source.getSequence('chrI', 2, 6)).toBe('AACC');
  });
});

describe('GFF3 Parser & Demo Dataset', () => {
  it('parses SGD GFF3 genes and supports gene search', () => {
    const gff = `##gff-version 3
chrV\tSGD\tgene\t1001\t1804\t.\t+\t.\tID=YEL021W;Name=YEL021W;gene=URA3;Note=Test%20note
chrV\tSGD\tCDS\t1001\t1804\t.\t+\t0\tParent=YEL021W
`;
    const idx = new AnnotationIndex(parseGff3Text(gff));
    const hits = idx.searchGenes('ura3');
    expect(hits).toHaveLength(1);
    expect(hits[0].start).toBe(1000);
    expect(hits[0].end).toBe(1804);
    expect(hits[0].subFeatures).toHaveLength(1);
  });

  it('creates S. cerevisiae demo dataset with all 17 chromosomes and URA3/chrM variants', async () => {
    const demo = createYeastDemoDataset();
    expect(demo.reference.getContigs()).toHaveLength(17);
    const refSeq = await demo.reference.getSequence('chrV', 1000, 1300);
    const s1Chrom = await demo.strain1.resolveContigName('chrV');
    const s2Chrom = await demo.strain2.resolveContigName('chrV');
    expect(s1Chrom).toBe('chrV');
    expect(s2Chrom).toBe('V'); // alias resolved!

    const s1Reads = await demo.strain1.getRecords(s1Chrom!, 1000, 1300);
    const s2Reads = await demo.strain2.getRecords(s2Chrom!, 1000, 1300);
    const analysis = buildThreeWayAnalysis(
      'chrV',
      1000,
      1300,
      refSeq,
      s1Reads,
      s2Reads,
      DEFAULT_CONSENSUS_OPTIONS,
    );

    expect(analysis.variants.length).toBeGreaterThanOrEqual(4);
    expect(analysis.ungappedS1.length).toBe(303); // +3bp insertion
    expect(analysis.ungappedS2.length).toBe(294); // -6bp deletion

    // chrM COX1 mitochondrial Table 3 synonymous check (TGA -> TGG = W -> W)
    const mRef = await demo.reference.getSequence('chrM', 500, 740);
    const mS1Chrom = await demo.strain1.resolveContigName('chrM');
    const mS2Chrom = await demo.strain2.resolveContigName('chrM');
    expect(mS2Chrom).toBe('MT');
    const mAnalysis = buildThreeWayAnalysis(
      'chrM',
      500,
      740,
      mRef,
      await demo.strain1.getRecords(mS1Chrom!, 500, 740),
      await demo.strain2.getRecords(mS2Chrom!, 500, 740),
      DEFAULT_CONSENSUS_OPTIONS,
    );
    const mitoProt = compareExtractedProteins(
      mAnalysis.ungappedRef,
      mAnalysis.ungappedS1,
      mAnalysis.ungappedS2,
      3,
      0,
    );
    expect(mitoProt.s1Status).toBe('Synonymous');
  });
});

describe('fileMatcher', () => {
  it('pairs .bam with .bam.bai and .bai automatically and detects missing index', () => {
    const files = [
      new File([''], 'S288C.fasta'),
      new File([''], 'S288C.fasta.fai'),
      new File([''], 'sgd.gff3'),
      new File([''], 'strain1.bam'),
      new File([''], 'strain1.bam.bai'),
      new File([''], 'strain2.bam'),
      new File([''], 'strain2.bai'),
      new File([''], 'strain3.bam'),
    ];
    const matched = matchDroppedFiles(files);
    expect(matched.fastaFile?.name).toBe('S288C.fasta');
    expect(matched.faiFile?.name).toBe('S288C.fasta.fai');
    expect(matched.gffFile?.name).toBe('sgd.gff3');
    expect(matched.pairedStrains).toHaveLength(2);
    expect(matched.unpairedBams.map((f) => f.name)).toEqual(['strain3.bam']);
  });
});

describe('Single-strain 2-way comparison, full S. cerevisiae aliases, irregular FASTA, and built-in S288C genes', () => {
  it('supports single-strain 2-way comparison when s2Reads is null', () => {
    const ref = 'ATGCGATCGATCGATC';
    const s1Reads: AlignedRead[] = Array.from({ length: 6 }, (_, i) => ({
      id: `r_${i}`,
      name: `r_${i}`,
      start: 0,
      end: 16,
      strand: '+',
      seq: 'ATACGATCGATCGATC', // SNP at pos 2 (G->A)
      cigar: '16M',
      mapq: 60,
      flags: 0,
    }));

    const analysis = buildThreeWayAnalysis(
      'chrI',
      0,
      16,
      ref,
      s1Reads,
      null,
      DEFAULT_CONSENSUS_OPTIONS,
    );

    expect(analysis.hasStrain2).toBe(false);
    expect(analysis.s2Calls).toEqual([]);
    expect(analysis.ungappedS2).toBe('');
    expect(analysis.alignedS2).toBe('');
    expect(analysis.variants).toHaveLength(1);
    expect(analysis.variants[0].category).toBe('strain1_only');
    expect(analysis.variants[0].type).toBe('SNP');
    expect(analysis.variants[0].start).toBe(2);
  });

  it('normalizes all S. cerevisiae NCBI/GenBank/SGD/zero-padded names and length-based fallback', () => {
    expect(normalizeYeastChromName('ref|NC_001133|')).toBe('chrI');
    expect(normalizeYeastChromName('gi|330443520|ref|NC_001140.6|')).toBe('chrVIII');
    expect(normalizeYeastChromName('BK006934.2')).toBe('chrVIII');
    expect(normalizeYeastChromName('BK006949.2')).toBe('chrXVI');
    expect(normalizeYeastChromName('chr01')).toBe('chrI');
    expect(normalizeYeastChromName('AJ011856.1')).toBe('chrM');

    // Every canonical S. cerevisiae chromosome resolves against NCBI RefSeq/Roman names
    const romanList = YEAST_CHROMOSOMES.map((c) => c.roman);
    for (const def of YEAST_CHROMOSOMES) {
      expect(resolveMatchingRefName(def.refseq, romanList)).toBe(def.roman);
    }

    // Length-based fallback for custom contig names
    const customContigs = [
      { name: 'scaf_1', length: 230218 },
      { name: 'scaf_8', length: 562643 },
    ];
    expect(
      resolveMatchingRefName('chrVIII', ['scaf_1', 'scaf_8'], undefined, 562643, customContigs),
    ).toBe('scaf_8');
  });

  it('reads FASTA accurately even with UTF-8 BOM and irregular line widths', async () => {
    const bomFasta = '\uFEFF>chrI\nATGCG\nATCGATCGATC\n\n>chrII\nAAAACCC\nCGGGG\n';
    const fastaFile = new File([bomFasta], 'irregular.fasta');
    const source = await IndexedFastaSource.fromFiles(fastaFile);

    expect(source.getContigs()).toEqual([
      { name: 'chrI', length: 16 },
      { name: 'chrII', length: 12 },
    ]);
    expect(await source.getSequence('chrI', 3, 12)).toBe('CGATCGATC');
    expect(await source.getSequence('chrII', 2, 10)).toBe('AACCCCGG');
  });

  it('searches built-in S288C R64 genes including URA3, ACT1, and COX1', () => {
    const idx = createBuiltInS288CAnnotationIndex();
    expect(idx.totalCount).toBeGreaterThan(7000);

    const ura3 = idx.searchGenes('URA3');
    expect(ura3[0].name).toBe('URA3');
    expect(ura3[0].chrom).toBe('chrV');

    const cox1 = idx.searchGenes('COX1');
    expect(cox1[0].name).toBe('COX1');
    expect(cox1[0].chrom).toBe('chrM');

    const act1 = idx.searchGenes('ACT1');
    expect(act1[0].name).toBe('ACT1');
    expect(act1[0].subFeatures.length).toBeGreaterThan(1); // spliced yeast gene!
  });

  it('sorts S. cerevisiae chromosomes in canonical biological order (I..XVI, chrM)', () => {
    const unsorted = [
      { name: 'chrM', length: 85779 },
      { name: 'chrIX', length: 439888 },
      { name: 'chrV', length: 576874 },
      { name: 'chrI', length: 230218 },
      { name: 'chrII', length: 813184 },
    ];
    expect(sortYeastContigs(unsorted).map((c) => c.name)).toEqual([
      'chrI',
      'chrII',
      'chrV',
      'chrIX',
      'chrM',
    ]);
  });

  it('packs overlapping and tightly adjacent genes onto separate non-colliding lanes', () => {
    const idx = createBuiltInS288CAnnotationIndex();
    const feats = idx.getFeaturesInRange('chrV', 435_890, 438_397);
    expect(feats.length).toBeGreaterThanOrEqual(2);

    const { packed, totalRows } = packGeneFeatures(feats, 0.4);
    expect(totalRows).toBeGreaterThanOrEqual(2);

    // Every pair on the same row must have non-overlapping coordinate ranges
    for (let i = 0; i < packed.length; i++) {
      for (let j = i + 1; j < packed.length; j++) {
        if (packed[i].row !== packed[j].row) continue;
        const a = packed[i].feature;
        const b = packed[j].feature;
        expect(a.end <= b.start || b.end <= a.start).toBe(true);
      }
    }

    // Also test explicit overlapping synthetic annotations on opposite strands
    const overlapping = [
      {
        id: 'g1',
        name: 'YER133W-A',
        chrom: 'chrV',
        start: 437_200,
        end: 438_350,
        strand: '+' as const,
        type: 'gene',
        subFeatures: [],
      },
      {
        id: 'g2',
        name: 'YER135C',
        chrom: 'chrV',
        start: 437_600,
        end: 438_500,
        strand: '-' as const,
        type: 'gene',
        subFeatures: [],
      },
      {
        id: 'g3',
        name: 'MDP-1',
        chrom: 'chrV',
        start: 436_300,
        end: 437_195,
        strand: '-' as const,
        type: 'gene',
        subFeatures: [],
      },
    ];
    const syn = packGeneFeatures(overlapping, 0.45);
    const rowsByName = new Map(syn.packed.map((p) => [p.feature.name, p.row]));
    expect(rowsByName.get('YER133W-A')).not.toBe(rowsByName.get('YER135C'));
    expect(syn.totalRows).toBeGreaterThanOrEqual(2);
  });

  it('matches IGV RefSeq CDS bounds for YER134C and YER135C and places them on one non-overlapping row', () => {
    const idx = createBuiltInS288CAnnotationIndex();
    const y134 = idx.searchGenes('YER134C').find((g) => g.name === 'YER134C');
    const byAlias = idx.searchGenes('MDP-1');
    const y135 = idx.searchGenes('YER135C').find((g) => g.name === 'YER135C');

    expect(y134).toBeDefined();
    expect(y134?.start).toBe(437_266);
    expect(y134?.end).toBe(437_803);
    expect(y134?.strand).toBe('-');
    expect(byAlias.some((g) => g.name === 'YER134C')).toBe(true);

    expect(y135).toBeDefined();
    expect(y135?.start).toBe(437_951);
    expect(y135?.end).toBe(438_344);
    expect(y135?.strand).toBe('-');

    const refseqVisible = idx
      .getFeaturesInRange('chrV', 436_267, 439_105)
      .filter((f) => f.type !== 'dubious_gene');
    const packed = packGeneFeatures(refseqVisible, 0.45);
    const p134 = packed.packed.find((p) => p.feature.name === 'YER134C');
    const p135 = packed.packed.find((p) => p.feature.name === 'YER135C');
    expect(p134?.row).toBe(0);
    expect(p135?.row).toBe(0);
  });

  it('trims SGD GFF3 UTR gene bounds to exact CDS intervals and flags Dubious ORFs', () => {
    const gff = `##gff-version 3
chrV\tSGD\tgene\t436308\t437979\t.\t-\t.\tID=YER134C;Name=YER134C;Alias=MDP-1,Mg-dependent%20acid%20phosphatase
chrV\tSGD\tCDS\t437267\t437803\t.\t-\t0\tParent=YER134C_id001;orf_classification=Verified
chrV\tSGD\tmRNA\t436308\t437979\t.\t-\t.\tID=YER134C_id001;Parent=YER134C
chrV\tSGD\tgene\t437191\t437532\t.\t+\t.\tID=YER133W-A;Name=YER133W-A
chrV\tSGD\tCDS\t437191\t437532\t.\t+\t0\tParent=YER133W-A_mRNA;orf_classification=Dubious
chrV\tSGD\tmRNA\t437191\t437532\t.\t+\t.\tID=YER133W-A_mRNA;Parent=YER133W-A
`;
    const feats = parseGff3Text(gff);
    const yer134c = feats.find((f) => f.name === 'YER134C');
    const dubious = feats.find((f) => f.name === 'YER133W-A');
    expect(yer134c?.start).toBe(437_266);
    expect(yer134c?.end).toBe(437_803);
    expect(yer134c?.note).toContain('MDP-1');
    expect(dubious?.type).toBe('dubious_gene');
  });

  it('filters secondary (0x100), QC-fail (0x200), and duplicate (0x400) alignments in consensus', () => {
    const ref = 'ATGCGATCGATCGATC';
    const goodReads: AlignedRead[] = Array.from({ length: 5 }, (_, i) => ({
      id: `good_${i}`,
      name: `good_${i}`,
      start: 0,
      end: 16,
      strand: '+',
      seq: ref,
      cigar: '16M',
      mapq: 60,
      flags: 0x43,
    }));
    const excludedReads: AlignedRead[] = [0x100, 0x200, 0x400].flatMap((flag, idx) =>
      Array.from({ length: 4 }, (_, j) => ({
        id: `bad_${idx}_${j}`,
        name: `bad_${idx}_${j}`,
        start: 0,
        end: 16,
        strand: '+',
        seq: 'TTGCGATCGATCGATC', // A->T SNP only in excluded flags
        cigar: '16M',
        mapq: 60,
        flags: flag,
      })),
    );

    const calls = computeConsensusForRegion(
      ref,
      [...goodReads, ...excludedReads],
      0,
      16,
      DEFAULT_CONSENSUS_OPTIONS,
    );
    expect(calls[0].depth).toBe(5);
    expect(calls[0].call).toBe('A');
  });

  it('packs 400x high-depth reads with min-heap keeping <=28 visible rows and exact totalRows', () => {
    const denseReads: AlignedRead[] = Array.from({ length: 400 }, (_, i) => ({
      id: `d_${i}`,
      name: `d_${i}`,
      start: 100 + (i % 4) * 5,
      end: 250 + (i % 4) * 5,
      strand: '+',
      seq: 'A'.repeat(150),
      cigar: '150M',
      mapq: 60,
      flags: 0x43,
    }));

    const { visiblePacked, totalRows } = packReadsWithHeap(denseReads, 28);
    expect(totalRows).toBe(400);
    expect(visiblePacked.every((r) => r.row < 28)).toBe(true);
    expect(visiblePacked).toHaveLength(28);
  });

  it('parses CIGAR without regex and derives allele ticks cleanly at 20% threshold', () => {
    expect(parseCigar('5S95M1I3M2D2H')).toEqual([
      { op: 'S', len: 5 },
      { op: 'M', len: 95 },
      { op: 'I', len: 1 },
      { op: 'M', len: 3 },
      { op: 'D', len: 2 },
      { op: 'H', len: 2 },
    ]);

    const ref = 'ATGCGATCGATCGATC';
    const reads: AlignedRead[] = [
      ...Array.from({ length: 8 }, (_, i) => ({
        id: `r_${i}`,
        name: `r_${i}`,
        start: 0,
        end: 16,
        strand: '+' as const,
        seq: ref,
        cigar: '16M',
        mapq: 60,
        flags: 0,
      })),
      ...Array.from({ length: 2 }, (_, i) => ({
        id: `v_${i}`,
        name: `v_${i}`,
        start: 0,
        end: 16,
        strand: '+' as const,
        seq: 'TTGCGATCGATCGATC', // 20% T at pos 0
        cigar: '16M',
        mapq: 60,
        flags: 0,
      })),
    ];

    const compact = buildThreeWayAnalysis('chrI', 0, 16, ref, reads, null, DEFAULT_CONSENSUS_OPTIONS, true);
    expect(compact.columns).toEqual([]);
    expect(compact.s1Calls).toHaveLength(16);
    const ticks = deriveAlleleTicksFromCalls(compact.s1Calls);
    expect(ticks.some((t) => t.pos === 0 && t.counts.T === 2)).toBe(true);
  });

  it('compares full >3kb coding DNA and full >1000aa protein sequences without 12-change truncation (e.g. T511I)', () => {
    // 1036 codons = 3108 bp (like MET10): starts ATG, T (ACT) at codon 511, ends TAA
    const codons = Array.from({ length: 1036 }, (_, idx) => {
      if (idx === 0) return 'ATG';
      if (idx === 510) return 'ACT'; // T -> will mutate to ATT (I)
      if (idx === 1035) return 'TAA';
      return idx % 2 === 0 ? 'GCT' : 'GTT'; // A / V
    });
    const refDna = codons.join('');
    expect(refDna.length).toBe(3108);

    // Strain 1: 15 missense mutations across the gene including T511I (>12 changes should NOT truncate)
    const s1Codons = [...codons];
    const missenseIndices = [5, 20, 40, 80, 120, 180, 240, 300, 360, 420, 480, 510, 700, 900, 1020];
    for (const ci of missenseIndices) {
      s1Codons[ci] = ci === 510 ? 'ATT' : 'TTT';
    }
    const s1Dna = s1Codons.join('');

    // Strain 2: only synonymous codon at 100 and T511I
    const s2Codons = [...codons];
    s2Codons[100] = 'GCC'; // A -> A synonymous
    s2Codons[510] = 'ATT'; // T511I
    const s2Dna = s2Codons.join('');

    const cmp = compareExtractedProteins(refDna, s1Dna, s2Dna, 1, 0);
    expect(cmp.refAA.length).toBe(1036);
    expect(cmp.s1AA.length).toBe(1036);
    expect(cmp.s2AA.length).toBe(1036);
    expect(cmp.s1Status).toBe('Missense');
    expect(cmp.s1Changes).toHaveLength(15);
    expect(cmp.s1Changes).toContain('T511I');
    expect(cmp.s2Changes).toEqual(['T511I']);
    expect(cmp.s2SynonymousPositions).toEqual([101]);

    const t511iDetail = cmp.s2Details.find((d) => d.label === 'T511I');
    expect(t511iDetail).toMatchObject({
      aaPos1: 511,
      refAA: 'T',
      strainAA: 'I',
      refCodon: 'ACT',
      strainCodon: 'ATT',
      isAmbiguousGap: false,
    });
  });

  it('splices multi-exon minus-strand CDS cleanly, reverses orientation, and reports full mapping coverage stats', () => {
    // Minus-strand 2-exon gene: exon1 [0, 6) + intron [6, 12) + exon2 [12, 18)
    // Coding 5'->3' on minus strand should be ATG (rc CAT at [15,18)), GCT (rc AGC at [12,15)), GAA (rc TTC at [3,6)), TAA (rc TTA at [0,3))
    // Forward genomic: Exon1 = TTA TTC ('TTATTC', rc='GAATAA'), Intron = GGGGGG, Exon2 = AGC CAT ('AGCCAT', rc='ATGGCT')
    // Reverse complement of Exon1 + Exon2 = rc('TTATTCAGCCAT') = 'ATGGCTGAATAA' => M A E *
    const refSeq = 'TTATTCGGGGGGAGCCAT';
    const mutSeq = 'TTATTCGGGGGGAGTCAT'; // Exon2 pos 14 (G->T), so rc becomes ATG GCT ... Wait: 'AGTCAT' rc='ATGACT' => T (M T E *)
    const reads1: AlignedRead[] = Array.from({ length: 6 }, (_, i) => ({
      id: `ok_${i}`,
      name: `ok_${i}`,
      start: 0,
      end: 18,
      strand: '+',
      seq: mutSeq,
      cigar: '18M',
      mapq: 60,
      flags: 0,
    }));
    // Strain 2 covers only exon1 ([0,6)) -> exon2 is unmapped low-depth gap
    const reads2: AlignedRead[] = Array.from({ length: 6 }, (_, i) => ({
      id: `partial_${i}`,
      name: `partial_${i}`,
      start: 0,
      end: 6,
      strand: '+',
      seq: refSeq.slice(0, 6),
      cigar: '6M',
      mapq: 60,
      flags: 0,
    }));

    const fullAnalysis = buildThreeWayAnalysis(
      'chrI',
      0,
      18,
      refSeq,
      reads1,
      reads2,
      DEFAULT_CONSENSUS_OPTIONS,
    );
    const spliced = sliceSplicedAnalysis(fullAnalysis, [
      { start: 0, end: 6 },
      { start: 12, end: 18 },
    ]);

    expect(spliced.ungappedRef).toBe('TTATTCAGCCAT');
    expect(spliced.ungappedS1).toBe('TTATTCAGTCAT');

    const rcRef = reverseComplement(spliced.ungappedRef);
    const rcS1 = reverseComplement(spliced.ungappedS1);
    const rcS2 = reverseComplement(spliced.ungappedS2);
    expect(rcRef).toBe('ATGGCTGAATAA');
    expect(rcS1).toBe('ATGACTGAATAA');

    const ev1 = evaluateStrainProtein(rcRef, rcS1, 1, 0);
    expect(ev1.aa).toBe('MTE*');
    expect(ev1.changes).toEqual(['A2T']);

    const s1Stats = computeRegionMappingStats(spliced.s1Calls, spliced.refSeq);
    const s2Stats = computeRegionMappingStats(spliced.s2Calls, spliced.refSeq);
    expect(s1Stats.coveragePct).toBe(100);
    expect(s1Stats.meanDepth).toBe(6);
    expect(s2Stats.coveragePct).toBe(50);
    expect(s2Stats.lowDepthBp).toBe(6);

    const ev2 = evaluateStrainProtein(rcRef, rcS2, 1, 0);
    expect(ev2.ambiguousXCount).toBe(2); // 6 unmapped N bases -> 2 X codons
  });

  it('produces identical chunked consensus calls and boundary insertions via buildAnalysisFromCalls', () => {
    const ref = 'ATGCGATCGATCGATC';
    const s1Seq = 'ATGCGATTGCGATCGATC'; // 2bp insertion 'TT' anchored at pos 5 (chunk boundary [0,6) vs [6,16))
    const reads: AlignedRead[] = Array.from({ length: 6 }, (_, i) => ({
      id: `c_${i}`,
      name: `c_${i}`,
      start: 0,
      end: 16,
      strand: '+',
      seq: s1Seq,
      cigar: '6M2I10M',
      mapq: 60,
      flags: 0,
    }));

    const single = buildThreeWayAnalysis('chrI', 0, 16, ref, reads, null, DEFAULT_CONSENSUS_OPTIONS);
    const cA = computeConsensusForRegion(ref.slice(0, 6), reads, 0, 6, DEFAULT_CONSENSUS_OPTIONS);
    const cB = computeConsensusForRegion(ref.slice(6, 16), reads, 6, 16, DEFAULT_CONSENSUS_OPTIONS);
    const chunked = buildAnalysisFromCalls('chrI', 0, 16, ref, [...cA, ...cB], null);

    expect(chunked.alignedRef).toBe(single.alignedRef);
    expect(chunked.alignedS1).toBe(single.alignedS1);
    expect(chunked.ungappedS1).toBe(single.ungappedS1);
    expect(chunked.variants).toEqual(single.variants);
  });

  it('auto-recovers when external .fai has mismatched CRLF/offset and supports calmd (=) BAM bases', async () => {
    // Windows CRLF FASTA given a stale LF-based .fai
    const crlfFasta = '>chrI\r\nATGC\r\nGATC\r\n';
    const staleLfFai = 'chrI\t8\t6\t4\t5\n'; // wrong offset (6 vs 7) and lineWidth (5 vs 6)
    const src = await IndexedFastaSource.fromFiles(
      new File([crlfFasta], 'win.fa'),
      new File([staleLfFai], 'win.fa.fai'),
    );
    expect(await src.getSequence('chrI', 0, 8)).toBe('ATGCGATC');

    // calmd -e uses '=' for bases matching Reference
    const calmdReads: AlignedRead[] = Array.from({ length: 6 }, (_, i) => ({
      id: `calmd_${i}`,
      name: `calmd_${i}`,
      start: 0,
      end: 8,
      strand: '+',
      seq: '==A=====', // SNP at pos 2 (G->A), rest match Ref via '='
      cigar: '8M',
      mapq: 60,
      flags: 0,
    }));
    const calls = computeConsensusForRegion('ATGCGATC', calmdReads, 0, 8, DEFAULT_CONSENSUS_OPTIONS);
    expect(calls.map((c) => c.call).join('')).toBe('ATACGATC');

    // Sequential batches with same filename yield distinct strain IDs
    const batch1 = matchDroppedFiles([new File([''], 's.bam'), new File([''], 's.bam.bai')]);
    const batch2 = matchDroppedFiles([new File([''], 's.bam'), new File([''], 's.bam.bai')]);
    expect(batch1.pairedStrains[0].id).not.toBe(batch2.pairedStrains[0].id);
  });
});
