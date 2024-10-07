```bash
mkdir -p $outdir/genomescope/tmp
cd $outdir/genomescope
FastK \
	-k31 \								# kmer size
	-T8 \								# Threads
	-t1 \								# Produce table of sorted k-mers & counts
	-Ptmp \								# Directory for temporary files
	-N$sample.k31 \							# Output prefix
	$reads
```

```bash
Histex \
	-h1:32767 \ 							# Set range of counts for histogram
	-G \ 								# Output an ASCII format histogram especially for GeneScope.FK
	$specimen.k31 \ 						# Prefix of input hist file
	| tee k31.hist.txt \ 						# Save output of previous command while passing it on to the next
	| docker run dovetailg/genomescope \
		genomescope.R \
			-o . \ 						# Save output in current directory
			-n $sample.k31 \ 				# Output prefix
			-k31 # kmer size
```

```bash
KatGC \
	-T8 \								# Threads
	$sample.k31 \							# Input prefix
	$sample.k31_gc							# Output prefix
```

```bash
PloidyPlot \
	-T8 \								# Threads
	-Ptmp \								# Directory for temporary files
	-kv \								# Keep het-mer table for re-use / Verbose mode
	-o$sample.k31_ploidy \						# Output prefix
	$sample.k31 \
	2>k31_ploidy.log
```

```bash
cd /$outdir
hifiasm \
	-t 8 \								# Threads
	-o $sample \							# Output prefix					
	--primary \							# Output primary and alternate assemblies
	$reads
```

```bash
mkdir -p $outdir/$sample.p_ctg.ccs.merquryk
cd $outdir/$sample.p_ctg.ccs.merquryk
MerquryFK \
	-T6 \								# Threads
	$outdir/genomescope/$sample.k31 \				# Prefix for kmer counts
	$outdir/$sample.p_ctg.fa.gz \					# Primary assembly
	$outdir/$sample.a_ctg.fa.gz \					# Alternate assembly
	$sample.ccs							# Output prefix
```

```bash
docker run staphb/busco busco \
	busco \
		--metaeuk \						# We'll use metaEuk for gene discovery as it is much less memory intensive compared to miniProt							
		--tar \							# Compress subdirectories to save space
		--in $outdir/$sample.p_ctg.fa \				# Input assembly
		--cpu 8 \						# Threads
		--out $specimen.p_ctg.basidiomycota_odb10.busco \	# Output directory
		--mode genome \						# Type of assembly
		--lineage_dataset basidiomycota_odb10			# Lineage database to query
```

```bash
samtools view \							# **Decompress file**
	-u $hic_cram \
	| samtools fastq \					# **Convert to fastq**
		-F0xB00 \					# Filter to remove supplementary alignments, reads not passing filters, and secondary alignments
		-nt \						# Don't append read pair info to read names / Copy header lines
		- \
	| bwa-mem2 mem \					# **Align reads to assembly**
		-t8 \						# Threads
		-5 \						# For split alignment, take the alignment with the smallest coordinate as primary
		-S \						# Skip mate rescue
		-P \						# Skip pairing
		-C \						# Append FASTA/FASTQ comment to SAM output
		-p \						# Smart pairing
		$rg_lines \					# Insert retained header lines
		$outdir/purging/$sample.purged.fa - \
	| samtools fixmate \				# **Fix mate informatio**
		-m \						# Add mate score tag
		-p \						# Disable FR proper pair check
		-u \						# Uncompressed output
		- - \
	| samtools sort \				# **Sort alignment**
		--write-index \
		-l1 \						# Set compression level
		-@8 \						# Threads
		-T $outdir/scaffolding/$sample.sort.tmp \	# Temporary output
		-o $outdir/scaffolding/$sample.bam \		# Output file
		- 	
```

```bash
samtools markdup \
	--write-index \
	-c \							# Clear previous duplicate settings and tags.
	-@8 \							# Threads
	-T $outfile.mkdup.tmp \					# Temporary output
	-f $outfile.metrics.txt \				# Output stats
	$outdir/scaffolding/$sample.bam \
	$outdir/scaffolding/$sample.mkdup.bam
```

```bash
samtools stats \
	-@8 \							# Threads
	-F0xB00 \						# Filter to remove supplementary alignments, reads not passing filters, and secondary alignments
	$outdir/scaffolding/$sample.mkdup.bam \
	> $outdir/scaffolding/$sample.mkdup.bam.stats
plot_bamstats \
	-p $outdir/scaffolding/ 				# Output prefix
	$outdir/scaffolding/$sample.mkdup.bam.stats
```
