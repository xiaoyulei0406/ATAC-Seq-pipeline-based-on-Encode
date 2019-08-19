import "atacseq.wdl" as atac

workflow ATAC {
	
	Array[Array[Array[File]]] sample_info
	Array[String] sample_name
	String? adapter
	Array[String] adapters_R1
	Array[String] adapters_R2

	# temporary 2-dim adapters array [rep_id][merge_id]
	#Array[Array[String]] adapters_R1 = 
	#Array[Array[String]] adapters_R2 = 


	# endedness for input data
	Boolean? paired_end				# to define endedness for all replciates
									#	if defined, this will override individual endedness below
	Array[Boolean] paired_ends	# to define endedness for individual replicate

	# individual genome parameters
	File bowtie2_idx_tar
	File Ref_fa
	File chrsz
	String gensz	
	File blacklist

	# individual genome parameters for ATAqC
	File? tss 						# TSS BED file
	File? dnase 					# open chromatin region BED file
	File? prom 						# promoter region BED file
	File? enh 						# enhancer region BED file
	File? reg2map 					# file with cell type signals
	File? reg2map_bed 				# file of regions used to generate reg2map signals
	File? roadmap_meta 				# roadmap metedata

	### pipeline type
	# parameters for pipeline
	String pipeline_type = 'atac'	# atac (default), dnase
									# tn5 shiting will be enabled for atac only
	Boolean need_trimandmerge = false 
									#	start from multi lanes
	Boolean only_need_trim = true	#	start from replications	
	Boolean no_trim = false			
	Boolean align_only = false		# disable all post-align analyses (peak-calling, overlap, idr, ...)
	Boolean true_rep_only = false 	# disable all analyses for pseudo replicates
									# 	if activated, overlap and idr will be disabled

	# parameters for trim_adapter

	Boolean auto_detect_adapter = false
									# automatically detect/trim adapters
									# 	can detect three adapters only
									# 	see /src/detect_adapter.py for details
	String cutadapt_param = '-e 0.1 -m 5' 
									# cutadapt parameters (err_rate=0.1, min_trim_len=5)
	## parameters for align (align FASTQs and create raw BAM)
	#String aligner = 'bowtie2' 		# bowtie2, custom
	Int multimapping = 4			# for samples with multimapping reads
	String bowtie2_param_se = '--local'
									# params for bowtie2 (single-ended samples)
	String bowtie2_param_pe = '-X2000 --mm --local' 
									# params for bowtie2 (paired-ended samples)

	# parameters for filter (filter/dedup raw BAM)
	String dup_marker = 'picard'	# picard.jar MarkDuplicates (picard) or 
									# sambamba markdup (sambamba)
	Int mapq_thresh = 30			# threshold for low MAPQ reads removal
	Boolean no_dup_removal = false	# no dupe reads removal when filtering BAM
									# dup.qc and pbc.qc will be empty files
									# and nodup_bam in the output is 
									# filtered bam with dupes

	# parameters for bam2ta (convert filtered/deduped BAM to TAG-ALIGN)
	String mito_chr_name = 'chrM' 	# name of mito chromosome. THIS IS NOT A REG-EX! you can define only one chromosome name for mito.
	String regex_filter_reads = 'chrM' 	# Perl-style regular expression pattern for chr name to filter out reads
                        			# those reads will be excluded from peak calling
	Int subsample_reads = 0			# number of reads to subsample TAGALIGN
									# 0 for no subsampling. this affects all downstream analysis

	# parameters for cross-correlation analysis
	Boolean enable_xcor = false 	# enable cross-corr analysis
	Int xcor_subsample_reads = 25000000	
									# number of reads to subsample TAG-ALIGN
									# 	this will be used for cross-corr only
									# 	will not affect any downstream analyses

	# parameters for blacklist filtering peaks
	Boolean keep_irregular_chr_in_bfilt_peak = false 
									# peaks with irregular chr name will not be filtered out
									# 	in bfilt_peak (blacklist filtered peak) file
									# 	(e.g. chr1_AABBCC, AABR07024382.1, ...)
									# 	reg-ex pattern for "regular" chr name is chr[\dXY]+\b

	# parameters for peak calling
	String peak_type = 'narrowPeak'
	Int cap_num_peak = 300000		# cap number of raw peaks called
	Float pval_thresh = 0.01		# p.value threshold for peak caller
	Int smooth_win = 73				# size of smoothing window for peak caller
	
	# parameters for signal tracks
	Boolean enable_count_signal_track = false # generate count signal track

	# parameters for IDR
	Boolean enable_idr = true 		# enable IDR analysis on raw peaks
	Float idr_thresh = 0.05			# IDR threshold
	String idr_rank = 'p.value' 	# IDR ranking method (p.value, q.value, score)


	# parameters for ATAqC
	Boolean disable_ataqc = false 	# disable ATAqC (extra annotation-based analysis)

		# resources 	
	#	these variables will be automatically ignored if they are not supported by platform
	# 	"disks" is for cloud platforms (Google Cloud Platform, DNAnexus) only
	Int trim_adapter_cpu = 2
	Int trim_adapter_mem_mb = 12000
	Int trim_adapter_time_hr = 24
	String trim_adapter_disks = "local-disk 100 HDD"

	Int bowtie2_cpu = 4
	Int bowtie2_mem_mb = 20000
	Int bowtie2_time_hr = 48
	String bowtie2_disks = "local-disk 200 HDD"

	Int filter_cpu = 2
	Int filter_mem_mb = 20000
	Int filter_time_hr = 24
	String filter_disks = "local-disk 400 HDD"

	Int bam2ta_cpu = 2
	Int bam2ta_mem_mb = 10000
	Int bam2ta_time_hr = 6
	String bam2ta_disks = "local-disk 100 HDD"

	Int spr_mem_mb = 16000

	Int xcor_cpu = 2
	Int xcor_mem_mb = 16000
	Int xcor_time_hr = 6
	String xcor_disks = "local-disk 100 HDD"

	Int macs2_mem_mb = 16000
	Int macs2_time_hr = 24
	String macs2_disks = "local-disk 200 HDD"

	Int ataqc_mem_mb = 16000
	Int ataqc_mem_java_mb = 15000
	Int ataqc_time_hr = 24
	String ataqc_disks = "local-disk 200 HDD"

	####################### pipeline starts here #######################
	# DO NOT DEFINE ANY VARIABLES DECLARED BELOW IN AN INPUT JSON FILE #
	# THEY ARE TEMPORARY/INTERMEDIATE SYSTEM VARIABLES                 #
	####################### pipeline starts here #######################


	# align each replicate
	scatter (i in range(length(sample_info))) {
		###1.trim_adaper
		if (only_need_trim){
			call atac.only_trim { input :
				fq1 = sample_info[i][0][0],
				fq2 = sample_info[i][0][1],
				adapter = adapter,
				prefix = sample_name[i],
				adapters_R1 = adapters_R1[i],
				adapters_R2 = adapters_R2[i],
				paired_end = paired_end,
				auto_detect_adapter = auto_detect_adapter,
				cutadapt_param = cutadapt_param,
				# resource
				cpu = trim_adapter_cpu,
				mem_mb = trim_adapter_mem_mb,
				time_hr = trim_adapter_time_hr,
				disks = trim_adapter_disks,
			}
			###2.Alignment Bowtie2
			call atac.bowtie2 { input :
				fastq_R1 = only_trim.trimmed_fq1,
				fastq_R2 = only_trim.trimmed_fq2,
				paired_end = paired_end,
				#aligner = aligner,
				multimapping = multimapping,
				idx_tar = bowtie2_idx_tar,
				bowtie2_param_se = bowtie2_param_se,
				bowtie2_param_pe = bowtie2_param_pe,
				# resource
				cpu = bowtie2_cpu,
				mem_mb = bowtie2_mem_mb,
				time_hr = bowtie2_time_hr,
				disks = bowtie2_disks,
			}
		}
		if (no_trim){
			call atac.bowtie2 as bowtie_only { input :
				fastq_R1 = sample_info[i][0][0],
				fastq_R2 = sample_info[i][0][1],
				paired_end = paired_end,
				#aligner = aligner,
				multimapping = multimapping,
				idx_tar = bowtie2_idx_tar,
				bowtie2_param_se = bowtie2_param_se,
				bowtie2_param_pe = bowtie2_param_pe,
				# resource
				cpu = bowtie2_cpu,
				mem_mb = bowtie2_mem_mb,
				time_hr = bowtie2_time_hr,
				disks = bowtie2_disks,
			}
		}
		File bowtie_bam = if no_trim then bowtie_only.bam else bowtie2.bam
		###3. filtering	
		call atac.filter {input:
			bam = bowtie_bam,
			paired_end = paired_end,
			prefix = sample_name[i],
			dup_marker = dup_marker,
			mapq_thresh = mapq_thresh,
			multimapping = multimapping,
			no_dup_removal = no_dup_removal,
			mito_chr_name = mito_chr_name,

			cpu = bam2ta_cpu,
			mem_mb = bam2ta_mem_mb,
			time_hr = bam2ta_time_hr,
			disks = bam2ta_disks,
		}

		###4. bamtotagalign: convert bam to tagalign and subsample it if necessary
		call atac.bam2ta {input :
			bam = filter.nodup_bam,
			disable_tn5_shift = if pipeline_type=='atac' then false else true,
			paired_end = paired_end,
			subsample = subsample_reads,
			regex_grep_v_ta = regex_filter_reads,
			mito_chr_name = mito_chr_name,

			cpu = bam2ta_cpu,
			mem_mb = bam2ta_mem_mb,
			time_hr = bam2ta_time_hr,
			disks = bam2ta_disks,
		}
		# let's call peaks
		if ( enable_xcor ) {
			#5. subsample tagalign (non-mito) and cross-correlation analysis
			call atac.xcor { input :
				ta = bam2ta.ta,
				subsample = xcor_subsample_reads,
				paired_end = paired_end,
				mito_chr_name = mito_chr_name,

				cpu = xcor_cpu,
				mem_mb = xcor_mem_mb,
				time_hr = xcor_time_hr,
				disks = xcor_disks,
			}
		}
		#6. call peaks on tagalign
		call atac.macs2 { input :
			ta = bam2ta.ta,
			gensz = gensz,
			chrsz = chrsz,
			cap_num_peak = cap_num_peak,
			pval_thresh = pval_thresh,
			smooth_win = smooth_win,
			blacklist = blacklist,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,

			mem_mb = macs2_mem_mb,
			disks = macs2_disks,
			time_hr = macs2_time_hr,
		}
		#7. call MACS2 peaks for true replicates to get signal tracks 
		call atac.macs2_signal_track { input :
			ta = bam2ta.ta,
			gensz = gensz,
			chrsz = chrsz,
			pval_thresh = pval_thresh,
			smooth_win = smooth_win,

			mem_mb = macs2_mem_mb,
			disks = macs2_disks,
			time_hr = macs2_time_hr,
		}

		if ( enable_count_signal_track ) {
			#11. generate count signal track
			call atac.count_signal_track { input :
				ta = bam2ta.ta,
				chrsz = chrsz,
			}
		}
	}
	# if there are TAs for ALL replicates then pool them

	Array[File] tas_ = bam2ta.ta

	#8. make two self pseudo replicates per true replicate
	if ( true_rep_only == false) {
		scatter (ta in tas_) {
			call atac.spr { input :
				ta = ta,
				paired_end = paired_end,
				mem_mb = spr_mem_mb
			}
		}
	}
	if ( true_rep_only == false ) {
		scatter(i in range(length(tas_))) {
		#9. call peaks on 1st pseudo replicated tagalign 
			call atac.macs2 as macs2_pr1 { input :
				ta = spr.ta_pr1[i],
				gensz = gensz,
				chrsz = chrsz,
				cap_num_peak = cap_num_peak,
				pval_thresh = pval_thresh,
				smooth_win = smooth_win,
				blacklist = blacklist,
				keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,

				mem_mb = macs2_mem_mb,
				disks = macs2_disks,
				time_hr = macs2_time_hr,
			}
			#10. call peaks on 2nd pseudo replicated tagalign 
			call atac.macs2 as macs2_pr2 { input :
				ta = spr.ta_pr2[i],
				gensz = gensz,
				chrsz = chrsz,
				cap_num_peak = cap_num_peak,
				pval_thresh = pval_thresh,
				smooth_win = smooth_win,
				blacklist = blacklist,
				keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,

				mem_mb = macs2_mem_mb,
				disks = macs2_disks,
				time_hr = macs2_time_hr,
			}
		}
			
	}

	if ( length(tas_)>1 ) {
		# pool tagaligns from true replicates
		call atac.pool_ta { input :
			tas = tas_,
		}
	}

	if ( (true_rep_only==false) && length(tas_)>1 ) {
		# pool tagaligns from pseudo replicates
		call atac.pool_ta as pool_ta_pr1 { input :
			tas = spr.ta_pr1,
		}
		call atac.pool_ta as pool_ta_pr2 { input :
			tas = spr.ta_pr2,
		}
	}

	if ( length(tas_)>1 ) {
		# call peaks on pooled replicate
		# always call MACS2 peaks for pooled replicate to get signal tracks
		call atac.macs2 as macs2_pooled { input :
			ta = pool_ta.ta_pooled,
			gensz = gensz,
			chrsz = chrsz,
			cap_num_peak = cap_num_peak,
			pval_thresh = pval_thresh,
			smooth_win = smooth_win,
			blacklist = blacklist,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,

			mem_mb = macs2_mem_mb,
			disks = macs2_disks,
			time_hr = macs2_time_hr,
		}
	}
	if (enable_count_signal_track && length(tas_)>1 ) {
		call atac.count_signal_track as count_signal_track_pooled { input :
			ta = pool_ta.ta_pooled,
			gensz = gensz,
			chrsz = chrsz,
			pval_thresh = pval_thresh,
			smooth_win = smooth_win,

			mem_mb = macs2_mem_mb,
			disks = macs2_disks,
			time_hr = macs2_time_hr,
		}
	}

	if (length(tas_)>1){
		call atac.macs2_signal_track as macs2_signal_track_pooled { input :
			ta = pool_ta.ta_pooled,
			chrsz = chrsz,
			gensz = gensz,
			pval_thresh = pval_thresh,
			smooth_win = smooth_win,

			mem_mb = macs2_mem_mb,
			disks = macs2_disks,
			time_hr = macs2_time_hr,
		}
	}


	if ( (true_rep_only== false) && length(tas_)>1 ) {
		# call peaks on 1st pooled pseudo replicates
		call atac.macs2 as macs2_ppr1 { input :
			ta = pool_ta_pr1.ta_pooled,
			gensz = gensz,
			chrsz = chrsz,
			cap_num_peak = cap_num_peak,
			pval_thresh = pval_thresh,
			smooth_win = smooth_win,
			blacklist = blacklist,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,

			mem_mb = macs2_mem_mb,
			disks = macs2_disks,
			time_hr = macs2_time_hr,
		}
	}
	if ( (true_rep_only== false) && length(tas_)>1 ) {
		# call peaks on 2nd pooled pseudo replicates
		call atac.macs2 as macs2_ppr2 { input :
			ta = pool_ta_pr2.ta_pooled,
			gensz = gensz,
			chrsz = chrsz,
			cap_num_peak = cap_num_peak,
			pval_thresh = pval_thresh,
			smooth_win = smooth_win,
			blacklist = blacklist,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,

			mem_mb = macs2_mem_mb,
			disks = macs2_disks,
			time_hr = macs2_time_hr,
		}
	}


	# make peak arrays
	Array[File] true_peaks_ = macs2.npeak
#	Array[Pair[String,Array[File]]] true_peak_pairs = if length(true_peaks_)<=1 then [] else [('rep1-rep2',[true_peaks_[0],true_peaks_[1]])]
	# generate all possible pairs of true replicates (pair: left=prefix, right=[peak1,peak2])
	Array[Pair[String,Array[File]]] true_peak_pairs = if length(true_peaks_)<=1 then [] else if length(true_peaks_)<=2 then  [('rep1-rep2',[true_peaks_[0],true_peaks_[1]])] else if length(true_peaks_)<=3 then [('rep1-rep2',[true_peaks_[0],true_peaks_[1]]), ('rep1-rep3',[true_peaks_[0],true_peaks_[2]]), ('rep2-rep3',[true_peaks_[1],true_peaks_[2]])] else if length(true_peaks_)<=4 then [('rep1-rep2',[true_peaks_[0],true_peaks_[1]]), ('rep1-rep3',[true_peaks_[0],true_peaks_[2]]), ('rep1-rep4',[true_peaks_[0],true_peaks_[3]]), ('rep2-rep3',[true_peaks_[1],true_peaks_[2]]), ('rep2-rep4',[true_peaks_[1],true_peaks_[3]]), ('rep3-rep4',[true_peaks_[2],true_peaks_[3]])] else if length(true_peaks_)<=5 then [('rep1-rep2',[true_peaks_[0],true_peaks_[1]]), ('rep1-rep3',[true_peaks_[0],true_peaks_[2]]), ('rep1-rep4',[true_peaks_[0],true_peaks_[3]]), ('rep1-rep5',[true_peaks_[0],true_peaks_[4]]), ('rep2-rep3',[true_peaks_[1],true_peaks_[2]]), ('rep2-rep4',[true_peaks_[1],true_peaks_[3]]), ('rep2-rep5',[true_peaks_[1],true_peaks_[4]]), ('rep3-rep4',[true_peaks_[2],true_peaks_[3]]), ('rep3-rep5',[true_peaks_[2],true_peaks_[4]]), ('rep4-rep5',[true_peaks_[3],true_peaks_[4]])] else [('rep1-rep2',[true_peaks_[0],true_peaks_[1]]), ('rep1-rep3',[true_peaks_[0],true_peaks_[2]]), ('rep1-rep4',[true_peaks_[0],true_peaks_[3]]), ('rep1-rep5',[true_peaks_[0],true_peaks_[4]]), ('rep1-rep6',[true_peaks_[0],true_peaks_[5]]), ('rep2-rep3',[true_peaks_[1],true_peaks_[2]]), ('rep2-rep4',[true_peaks_[1],true_peaks_[3]]), ('rep2-rep5',[true_peaks_[1],true_peaks_[4]]), ('rep2-rep6',[true_peaks_[1],true_peaks_[5]]), ('rep3-rep4',[true_peaks_[2],true_peaks_[3]]), ('rep3-rep5',[true_peaks_[2],true_peaks_[4]]), ('rep3-rep6',[true_peaks_[2],true_peaks_[5]]), ('rep4-rep5',[true_peaks_[3],true_peaks_[4]]), ('rep4-rep6',[true_peaks_[3],true_peaks_[5]]), ('rep5-rep6',[true_peaks_[4],true_peaks_[5]])]
	
	scatter( pair in true_peak_pairs ) {
		# Naive overlap on every pair of true replicates
		call atac.overlap { input :
			prefix = pair.left,
			peak1 = pair.right[0],
			peak2 = pair.right[1],
			peak_pooled = macs2_pooled.npeak,
			peak_type = peak_type,
			blacklist = blacklist,
			chrsz = chrsz,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
			ta = pool_ta.ta_pooled,
		}
		if ( enable_idr ) {
			# IDR on every pair of true replicates
			call atac.idr { input :
				prefix = pair.left,
				peak1 = pair.right[0],
				peak2 = pair.right[1],
				peak_pooled = macs2_pooled.npeak,
				idr_thresh = idr_thresh,
				peak_type = peak_type,
				rank = idr_rank,
				blacklist = blacklist,
				chrsz = chrsz,
				keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
				ta = pool_ta.ta_pooled,
			}
		}
	}
	Array[File] peaks_pr1_ = macs2_pr1.npeak
	Array[File] peaks_pr2_ = macs2_pr2.npeak

	# overlap on pseudo-replicates (pr1, pr2) for each true replicate
	scatter( i in range(length(tas_)) ) {
		if ( (align_only ==false) && (true_rep_only== false) ) {
			call atac.overlap as overlap_pr { input :
				prefix = "rep"+(i+1)+"-pr",
				peak1 = peaks_pr1_[i],
				peak2 = peaks_pr2_[i],
				peak_pooled = true_peaks_[i],
				peak_type = peak_type,
				blacklist = blacklist,
				chrsz = chrsz,
				keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
				ta = tas_[i],
			}
		}
	}
	scatter( i in range(length(tas_))) {
		if ( (align_only==false) && (true_rep_only== false) && enable_idr ) {
			# IDR on pseduo replicates
			call atac.idr as idr_pr { input :
				prefix = "rep"+(i+1)+"-pr",
				peak1 = peaks_pr1_[i],
				peak2 = peaks_pr2_[i],
				peak_pooled = true_peaks_[i],
				idr_thresh = idr_thresh,
				peak_type = peak_type,
				rank = idr_rank,
				blacklist = blacklist,
				chrsz = chrsz,
				keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
				ta = tas_[i],
			}
		}
	}

	if ( (align_only==false) && (true_rep_only== false) && length(tas_)>1 ) {
		# Naive overlap on pooled pseudo replicates
		call atac.overlap as overlap_ppr { input :
			prefix = "ppr",
			peak1 = macs2_ppr1.npeak, #peak_ppr1_[0]
			peak2 = macs2_ppr2.npeak, #peak_ppr2_[0]
			peak_pooled = macs2_pooled.npeak,
			peak_type = peak_type,
			blacklist = blacklist,
			chrsz = chrsz,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
			ta = pool_ta.ta_pooled,
		}
	}

	if ( (align_only==false) && (true_rep_only== false) && length(tas_)>1  ) {
		# IDR on pooled pseduo replicates
		call atac.idr as idr_ppr { input :
			prefix = "ppr",
			peak1 = macs2_ppr1.npeak, #peak_ppr1_[0]
			peak2 = macs2_ppr2.npeak, #peak_ppr2_[0]
			peak_pooled = macs2_pooled.npeak,
			idr_thresh = idr_thresh,
			peak_type = peak_type,
			rank = idr_rank,
			blacklist = blacklist,
			chrsz = chrsz,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
			ta = pool_ta.ta_pooled,
		}
	}

		# reproducibility QC for overlap/IDR peaks
	if ( (align_only==false) && (true_rep_only== false) ) {
		# reproducibility QC for overlapping peaks
		call atac.reproducibility as reproducibility_overlap { input :
			prefix = 'overlap',
			peaks = overlap.bfilt_overlap_peak,
			peaks_pr = overlap_pr.bfilt_overlap_peak,
			peak_ppr = overlap_ppr.bfilt_overlap_peak,
			peak_type = peak_type,
			chrsz = chrsz,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
		}
	}

	if ( (align_only==false) && (true_rep_only== false) && enable_idr ) {
		# reproducibility QC for IDR peaks
		call atac.reproducibility as reproducibility_idr { input :
			prefix = 'idr',
			peaks = idr.bfilt_idr_peak,
			peaks_pr = idr_pr.bfilt_idr_peak,
			peak_ppr = idr_ppr.bfilt_idr_peak,
			peak_type = peak_type,
			chrsz = chrsz,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
		}
	}

	scatter( i in range(length(tas_))) {
		if ( (disable_ataqc==false) ){
			call atac.ataqc { input :
				paired_end = paired_end,
				read_len_log = bowtie2.read_len_log[i],
				flagstat_qc = bowtie2.flagstat_qc[i],
				bowtie2_log = bowtie2.align_log[i],
				pbc_qc = filter.pbc_qc[i],
				dup_qc = filter.dup_qc[i],
				bam = bowtie2.bam[i],
				nodup_flagstat_qc = filter.flagstat_qc[i],
				mito_dup_log = filter.mito_dup_log[i],
				nodup_bam = filter.nodup_bam[i],
				ta = tas_[i],
				peak = if defined(idr_pr.bfilt_idr_peak[i]) then idr_pr.bfilt_idr_peak[i] else reproducibility_overlap.optimal_peak,
				idr_peak = reproducibility_idr.optimal_peak,
				overlap_peak= reproducibility_overlap.optimal_peak,
				pval_bw = macs2_signal_track.pval_bw[i],
				ref_fa = Ref_fa,
				chrsz = chrsz,
				tss = tss,
				blacklist = blacklist,
				dnase = dnase,
				prom = prom,
				enh = enh,
				reg2map_bed = reg2map_bed,
				reg2map = reg2map,
				roadmap_meta = roadmap_meta,
				mito_chr_name = mito_chr_name,

				mem_mb = ataqc_mem_mb,
				mem_java_mb = ataqc_mem_java_mb,
				time_hr = ataqc_time_hr,
				disks = ataqc_disks,
			}
		}
	}

	# Generate final QC report and JSON
	call atac.qc_report { input :
		pipeline_ver = "v1",
		title = "atac",
		description = "encode_atac",
		genome = basename([Ref_fa]),
		multimapping = multimapping,
		paired_ends = paired_ends,
		pipeline_type = pipeline_type,
		peak_caller = 'macs2',
		macs2_cap_num_peak = cap_num_peak,
		idr_thresh = idr_thresh,
		
		flagstat_qcs = bowtie2.flagstat_qc,
		nodup_flagstat_qcs = filter.flagstat_qc,
		dup_qcs = filter.dup_qc,
		pbc_qcs = filter.pbc_qc,
		xcor_plots = xcor.plot_png,
		xcor_scores = xcor.score,

		frip_macs2_qcs = macs2.frip_qc,
		frip_macs2_qcs_pr1 = macs2_pr1.frip_qc,
		frip_macs2_qcs_pr2 = macs2_pr2.frip_qc,

		frip_macs2_qc_pooled = macs2_pooled.frip_qc,
		frip_macs2_qc_ppr1 = macs2_ppr1.frip_qc,
		frip_macs2_qc_ppr2 = macs2_ppr2.frip_qc,

		idr_plots = idr.idr_plot,
		idr_plots_pr = idr_pr.idr_plot,
		idr_plot_ppr = idr_ppr.idr_plot,
		frip_idr_qcs = idr.frip_qc,
		frip_idr_qcs_pr = idr_pr.frip_qc,
		frip_idr_qc_ppr = idr_ppr.frip_qc,
		frip_overlap_qcs = overlap.frip_qc,
		frip_overlap_qcs_pr = overlap_pr.frip_qc,
		frip_overlap_qc_ppr = overlap_ppr.frip_qc,
		idr_reproducibility_qc = reproducibility_idr.reproducibility_qc,
		overlap_reproducibility_qc = reproducibility_overlap.reproducibility_qc,

		ataqc_txts = ataqc.txt,
		ataqc_htmls = ataqc.html,
	}

	output {
		File report = qc_report.report
		File qc_json = qc_report.qc_json
		Boolean qc_json_ref_match = qc_report.qc_json_ref_match
	}
}
