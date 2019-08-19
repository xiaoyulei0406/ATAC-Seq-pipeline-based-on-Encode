task trim_adapter { 
	Array[File] fastqs_R1 		# [merge_id]
	Array[File] fastqs_R2
	Array[Array[File]] fastqs
	Array[Array[String]] adapters

	String? adapter 	# adapter for all fastqs,
						#	this will override individual adapters in adapters_R1/R2
	Array[String] adapters_R1
	Array[String] adapters_R2
	Boolean paired_end
	Boolean auto_detect_adapter
	String cutadapt_param 
	# resource
	Int cpu
	Int mem_mb
	Int time_hr
	String disks

	# tmp vars
	File? null_f
	Array[Array[File]] tmp_fastqs = if paired_end then transpose([fastqs_R1, fastqs_R2]) else transpose([fastqs_R1])
	Array[Array[String]] tmp_adapters = if paired_end then transpose([adapters_R1, adapters_R2]) else transpose([adapters_R1])
	command {
		python /opt/src/encode_trim_adapter.py \
			${write_tsv(tmp_fastqs)} \
			${"--adapter " + adapter} \
			--adapters ${write_tsv(tmp_adapters)} \
			${if paired_end then "--paired-end" else ""} \
			${if auto_detect_adapter then "--auto-detect-adapter" else ""} \
			--cutadapt-param ' ${cutadapt_param}' \
			${"--nth " + cpu}
	}
	output {
		File trim_merged_fastq_R1 = glob("R1/*.fastq.gz")[0]
		File? trim_merged_fastq_R2 = if paired_end then glob("R2/*.fastq.gz")[0] else null_f
	}
	runtime {
		docker:"encode_atac/ycl:v1"
		cpu:"10"
		memory:"20G"
	}
}

task only_trim {
  File fq1
  File fq2
  String cutadapt_param
  String adapters_R1
  String adapters_R2
  String prefix

  command {
    cutadapt ${cutadapt_param} -a ${adapters_R1} -A ${adapters_R2} ${fq1} ${fq2} -o ${prefix}_1.trim.fastq.gz -p ${prefix}_2.trim.fastq.gz
  }

  runtime {
    docker:"aws_tool:2.1.4"
    cpu:"2"
    memory:"10GB"
  }

  output {
    File trimmed_fq1 = "${prefix}_1.trim.fastq.gz"
    File trimmed_fq2 = "${prefix}_2.trim.fastq.gz"
  }
}

task bowtie2 {
	File idx_tar 		# reference bowtie2 index tar
	File? fastq_R1 		# [read_end_id]
	File? fastq_R2
	Boolean paired_end
	Int multimapping
	String bowtie2_param_se
	String bowtie2_param_pe
	Int cpu
	Int mem_mb
	Int time_hr
	String disks

	command {
		python /opt/src/encode_bowtie2.py \
			${idx_tar} \
			${fastq_R1} ${fastq_R2} \
			${if paired_end then "--paired-end" else ""} \
			${"--multimapping " + multimapping} \
			--bowtie2-param-se ' ${bowtie2_param_se}' \
			--bowtie2-param-pe ' ${bowtie2_param_pe}' \
			${"--nth " + cpu}
	}
	output {
		File bam = glob("*.bam")[0]
		File bai = glob("*.bai")[0]
		File align_log = glob("*.align.log")[0]
		File flagstat_qc = glob("*.flagstat.qc")[0]
		File read_len_log = glob("*.read_length.txt")[0] # read_len
	}
	runtime {
		docker:"encode_atac/ycl:v1"
		cpu:"10"
		memory:"20G"
	}
}

task filter {
	File bam
	Boolean paired_end
#	Int multimapping
	String prefix
	String dup_marker 			# picard.jar MarkDuplicates (picard) or 
								# sambamba markdup (sambamba)
	Int mapq_thresh				# threshold for low MAPQ reads removal
	Boolean no_dup_removal 		# no dupe reads removal when filtering BAM
	String mito_chr_name
	Int cpu
	Int mem_mb
	Int time_hr
	String disks

	command {
		python /opt/src/encode_filter.py \
			${bam} \
			${if paired_end then "--paired-end" else ""} \
			${"--multimapping " + 0} \
			${"--dup-marker " + dup_marker} \
			${"--mapq-thresh " + mapq_thresh} \
			${if no_dup_removal then "--no-dup-removal" else ""} \
			${"--mito-chr-name " + mito_chr_name} \
			${"--nth " + cpu}
	}
	output {
		File nodup_bam = "${prefix}_1.trim.nodup.bam"
		File nodup_bai = glob("*.bai")[0]
		File flagstat_qc = glob("*.flagstat.qc")[0]
		File dup_qc = glob("*.dup.qc")[0]
		File pbc_qc = glob("*.pbc.qc")[0]
		File mito_dup_log = glob("*.mito_dup.txt")[0] # mito_dups, fract_dups_from_mito
	}

	runtime {
		docker:"encode_atac/ycl:v1"
		cpu:"10"
		memory:"20G"
	}
}
task bam2ta {
	File bam
	Boolean paired_end
	Boolean disable_tn5_shift 	# no tn5 shifting (it's for dnase-seq)
	String regex_grep_v_ta   	# Perl-style regular expression pattern 
                        		# to remove matching reads from TAGALIGN
	String mito_chr_name 		# mito chromosome name
	Int subsample 				# number of reads to subsample TAGALIGN
								# this affects all downstream analysis
	Int cpu
	Int mem_mb
	Int time_hr
	String disks

	command {
		python /opt/src/encode_bam2ta.py \
			${bam} \
			${if paired_end then "--paired-end" else ""} \
			${if disable_tn5_shift then "--disable-tn5-shift" else ""} \
			${if regex_grep_v_ta!="" then "--regex-grep-v-ta '"+regex_grep_v_ta+"'" else ""} \
			${"--mito-chr-name " + mito_chr_name} \
			${"--subsample " + subsample} \
			${"--nth " + cpu}
	}
	output {
		File ta = glob("*.tagAlign.gz")[0]
	}

	runtime {
		docker:"encode_atac/ycl:v1"
		cpu:"10"
		memory:"20G"
	}
}

task spr { # make two self pseudo replicates
	File ta
	Boolean paired_end

	Int mem_mb

	command {
		python /opt/src/encode_spr.py \
			${ta} \
			${if paired_end then "--paired-end" else ""}
	}
	output {
		File ta_pr1 = glob("*.pr1.tagAlign.gz")[0]
		File ta_pr2 = glob("*.pr2.tagAlign.gz")[0]
	}

	runtime {
		docker:"encode_atac/ycl:v1"
		cpu:"10"
		memory:"20G"
	}
}

task pool_ta {
	# input variables
	Array[File?] tas 	# TAG-ALIGNs to be merged

	command {
		python /opt/src/encode_pool_ta.py \
			${sep=' ' tas}
	}
	output {
		File ta_pooled = glob("*.tagAlign.gz")[0]
	}

	runtime {
		docker:"encode_atac/ycl:v1"
		cpu:"10"
		memory:"20G"
	}
}

task xcor {
	File ta
	Boolean paired_end
	String mito_chr_name
	Int subsample  # number of reads to subsample TAGALIGN
				# this will be used for xcor only
				# will not affect any downstream analysis
	Int cpu
	Int mem_mb	
	Int time_hr
	String disks

	command {
		python /opt/src/encode_xcor.py \
			${ta} \
			${if paired_end then "--paired-end" else ""} \
			${"--mito-chr-name " + mito_chr_name} \
			${"--subsample " + subsample} \
			--speak=0 \
			${"--nth " + cpu}
	}
	output {
		File plot_pdf = glob("*.cc.plot.pdf")[0]
		File plot_png = glob("*.cc.plot.png")[0]
		File score = glob("*.cc.qc")[0]
		Int fraglen = read_int(glob("*.cc.fraglen.txt")[0])
	}

	runtime {
		docker:"encode_atac/ycl:v1"
		cpu:"10"
		memory:"20G"
	}
}

task count_signal_track {
	File ta 			# tag-align
	File chrsz			# 2-col chromosome sizes file

	command {
		python /opt/src/encode_count_signal_track.py \
			${ta} \
			${"--chrsz " + chrsz}
	}
	output {
		File pos_bw = glob("*.positive.bigwig")[0]
		File neg_bw = glob("*.negative.bigwig")[0]
	}

	runtime {
		docker:"encode_atac/ycl:v1"
		cpu:"10"
		memory:"20G"
	}
}

task macs2 {
	File ta
	String gensz		# Genome size (sum of entries in 2nd column of 
                        # chr. sizes file, or hs for human, ms for mouse)
	File chrsz			# 2-col chromosome sizes file
	Int cap_num_peak	# cap number of raw peaks called from MACS2
	Float pval_thresh  	# p.value threshold
	Int smooth_win 		# size of smoothing window
	File? blacklist 	# blacklist BED to filter raw peaks
	Boolean	keep_irregular_chr_in_bfilt_peak

	String? trim	
	Int mem_mb
	Int time_hr
	String disks

	File? null_f

	command {
		python /opt/src/encode_macs2_atac.py \
			${ta} \
			${"--gensz "+ gensz} \
			${"--chrsz " + chrsz} \
			${"--cap-num-peak " + cap_num_peak} \
			${"--pval-thresh "+ pval_thresh} \
			${"--smooth-win "+ smooth_win} \
			${if keep_irregular_chr_in_bfilt_peak then "--keep-irregular-chr" else ""} \
			${"--blacklist "+ blacklist}
	}
	output {
		File npeak = glob("*.300K.narrowPeak.gz")[0]
		File bfilt_npeak = glob("*.bfilt.narrowPeak.gz")[0]
		File bfilt_npeak_bb = glob("*.bfilt.narrowPeak.bb")[0]
		File bfilt_npeak_hammock = glob("*.bfilt.narrowPeak.hammock.gz*")[0]
		File bfilt_npeak_hammock_tbi = glob("*.bfilt.narrowPeak.hammock.gz*")[1]
		File frip_qc = glob("*.frip.qc")[0]
	}

	runtime {
		docker:"encode_atac/ycl:v1"
		cpu:"10"
		memory:"20G"
	}
}

task macs2_signal_track {
	File ta
	String gensz		# Genome size (sum of entries in 2nd column of 
                        # chr. sizes file, or hs for human, ms for mouse)
	File chrsz			# 2-col chromosome sizes file
	Float pval_thresh  	# p.value threshold
	Int smooth_win 		# size of smoothing window
	
	Int mem_mb
	Int time_hr
	String disks

	command {
		python /opt/src/encode_macs2_signal_track_atac.py \
			${ta} \
			${"--gensz "+ gensz} \
			${"--chrsz " + chrsz} \
			${"--pval-thresh "+ pval_thresh} \
			${"--smooth-win "+ smooth_win}
	}
	output {
		File pval_bw = glob("*.pval.signal.bigwig")[0]
		File fc_bw = glob("*.fc.signal.bigwig")[0]
	}

	runtime {
		docker:"encode_atac/ycl:v1"
		cpu:"10"
		memory:"30G"
	}
}

task idr {
	String prefix 		# prefix for IDR output file
	File peak1 			
	File peak2
	File peak_pooled
	Float idr_thresh
	File? blacklist 	# blacklist BED to filter raw peaks
	Boolean	keep_irregular_chr_in_bfilt_peak
	# parameters to compute FRiP
	File? ta			# to calculate FRiP
	File chrsz			# 2-col chromosome sizes file
	String peak_type
	String rank

	File? null_f

	command {
		${if defined(ta) then "" else "touch null.frip.qc"}
		touch null
		python /opt/src/encode_idr.py \
			${peak1} ${peak2} ${peak_pooled} \
			${"--prefix " + prefix} \
			${"--idr-thresh " + idr_thresh} \
			${"--peak-type " + peak_type} \
			--idr-rank ${rank} \
			${"--chrsz " + chrsz} \
			${"--blacklist "+ blacklist} \
			${if keep_irregular_chr_in_bfilt_peak then "--keep-irregular-chr" else ""} \
			${"--ta " + ta}
	}
	output {
		File idr_peak = "${prefix}.idr.${idr_thresh}.${peak_type}.gz"
		File bfilt_idr_peak = glob("*.bfilt."+peak_type+".gz")[0]
		File bfilt_idr_peak_bb = glob("*.bfilt."+peak_type+".bb")[0]
		File bfilt_idr_peak_hammock = glob("*.bfilt."+peak_type+".hammock.gz*")[0]
		File bfilt_idr_peak_hammock_tbi = glob("*.bfilt."+peak_type+".hammock.gz*")[1]
		File idr_plot = glob("*.txt.png")[0]
		File idr_unthresholded_peak = glob("*.txt.gz")[0]
		File idr_log = glob("*.idr*.log")[0]
		File frip_qc = if defined(ta) then glob("*.frip.qc")[0] else glob("null")[0]
	}

	runtime {
		docker:"encode_atac/ycl:v1"
		cpu:"10"
		memory:"20G"
	}
}

task overlap {
	String prefix 		# prefix for IDR output file
	File peak1
	File peak2
	File peak_pooled
	File? blacklist 	# blacklist BED to filter raw peaks
	Boolean	keep_irregular_chr_in_bfilt_peak
	File? ta		# to calculate FRiP
	File chrsz			# 2-col chromosome sizes file
	String peak_type

	File? null_f

	command {
		${if defined(ta) then "" else "touch null.frip.qc"}
		touch null 
		python /opt/src/encode_naive_overlap.py \
			${peak1} ${peak2} ${peak_pooled} \
			${"--prefix " + prefix} \
			${"--peak-type " + peak_type} \
			${"--chrsz " + chrsz} \
			${"--blacklist "+ blacklist} \
			--nonamecheck \
			${if keep_irregular_chr_in_bfilt_peak then "--keep-irregular-chr" else ""} \
			${"--ta " + ta}
	}
	output {
		File overlap_peak = "${prefix}.overlap.narrowPeak.gz"
		File bfilt_overlap_peak = "${prefix}.overlap.bfilt.narrowPeak.gz"
		File bfilt_overlap_peak_bb = "${prefix}.overlap.bfilt.narrowPeak.bb"
		File bfilt_overlap_peak_hammock = "${prefix}.overlap.bfilt.narrowPeak.hammock.gz"
		File frip_qc = "${prefix}.overlap.bfilt.frip.qc"
	}

	runtime {
		docker:"encode_atac/ycl:v1"
		cpu:"10"
		memory:"20G"
	}
}

task reproducibility {
	String prefix
	Array[File]? peaks # peak files from pair of true replicates
						# in a sorted order. for example of 4 replicates,
						# 1,2 1,3 1,4 2,3 2,4 3,4.
                        # x,y means peak file from rep-x vs rep-y
	Array[File?] peaks_pr	# peak files from pseudo replicates
	File? peak_ppr			# Peak file from pooled pseudo replicate.
	String peak_type
	File chrsz			# 2-col chromosome sizes file
	Boolean	keep_irregular_chr_in_bfilt_peak

	command {
		python /opt/src/encode_reproducibility_qc.py \
			${sep=' ' peaks} \
			--peaks-pr ${sep=' ' peaks_pr} \
			${"--peak-ppr "+ peak_ppr} \
			--prefix ${prefix} \
			${"--peak-type " + peak_type} \
			${if keep_irregular_chr_in_bfilt_peak then "--keep-irregular-chr" else ""} \
			${"--chrsz " + chrsz}
	}
	output {
		File optimal_peak = glob("optimal_peak.*.gz")[0]
		File conservative_peak = glob("conservative_peak.*.gz")[0]
		File optimal_peak_bb = glob("optimal_peak.*.bb")[0]
		File conservative_peak_bb = glob("conservative_peak.*.bb")[0]
		File optimal_peak_hammock = glob("optimal_peak.*.hammock.gz*")[0]
		File optimal_peak_hammock_tbi = glob("optimal_peak.*.hammock.gz*")[1]
		File conservative_peak_hammock = glob("conservative_peak.*.hammock.gz*")[0]
		File conservative_peak_hammock_tbi = glob("conservative_peak.*.hammock.gz*")[1]
		File reproducibility_qc = glob("*reproducibility.qc")[0]
	}

	runtime {
		docker:"encode_atac/ycl:v1"
		cpu:"10"
		memory:"20G"
	}
}

# annotation-based analysis
task ataqc {
	Boolean paired_end
	File? read_len_log
	File? flagstat_qc
	File? bowtie2_log
	File? bam
	File? nodup_flagstat_qc
	File? mito_dup_log
	File? dup_qc
	File? pbc_qc
	File? nodup_bam
	File? ta
	File? peak
	File? idr_peak 
	File? overlap_peak
	File? pval_bw
	# from genome database
	File? ref_fa
	File? chrsz
	File? tss
	File? blacklist
	File? dnase
	File? prom
	File? enh
	File? reg2map_bed
	File? reg2map
	File? roadmap_meta
	String mito_chr_name

	Int mem_mb
	Int mem_java_mb
	Int time_hr
	String disks

	command {
		export _JAVA_OPTIONS="-Xms256M -Xmx${mem_java_mb}M -XX:ParallelGCThreads=1 $_JAVA_OPTIONS"

		python /opt/src/encode_ataqc.py \
			${if paired_end then "--paired-end" else ""} \
			${"--read-len-log " + read_len_log} \
			${"--flagstat-log " + flagstat_qc} \
			${"--bowtie2-log " + bowtie2_log} \
			${"--bam " + bam} \
			${"--nodup-flagstat-log " + nodup_flagstat_qc} \
			${"--mito-dup-log " + mito_dup_log} \
			${"--dup-log " + dup_qc} \
			${"--pbc-log " + pbc_qc} \
			${"--nodup-bam " + nodup_bam} \
			${"--ta " + ta} \
			${"--bigwig " + pval_bw} \
			${"--peak " + peak} \
			${"--idr-peak " + idr_peak} \
			${"--overlap-peak " + overlap_peak} \
			${"--ref-fa " + ref_fa} \
			${"--blacklist " + blacklist} \
			${"--chrsz " + chrsz} \
			${"--dnase " + dnase} \
			${"--tss " + tss} \
			${"--prom " + prom} \
			${"--enh " + enh} \
			${"--reg2map-bed " + reg2map_bed} \
			${"--reg2map " + reg2map} \
			${"--roadmap-meta " + roadmap_meta} \
			${"--mito-chr-name " + mito_chr_name}

	}
	output {
		File html = glob("*_qc.html")[0]
		File txt = glob("*_qc.txt")[0]
	}
	runtime {
		docker:"encode_atac/ycl:v1"
		cpu:"10"
		memory:"20G"
	}
}
# gather all outputs and generate 
# - qc.html		: organized final HTML report
# - qc.json		: all QCs
task qc_report {
	String pipeline_ver
 	String title
	String description
	String? genome	
	# workflow params
	Int multimapping
	Array[Boolean?] paired_ends
	String pipeline_type
	String peak_caller
	Int? macs2_cap_num_peak
	Int? spp_cap_num_peak
	Float idr_thresh
	# QCs
	Array[File?] flagstat_qcs
	Array[File?] nodup_flagstat_qcs
	Array[File?] dup_qcs
	Array[File?] pbc_qcs
	Array[File?] xcor_plots
	Array[File?] xcor_scores
	Array[File]? idr_plots
	Array[File?] idr_plots_pr
	File? idr_plot_ppr
	Array[File?] frip_macs2_qcs
	Array[File?] frip_macs2_qcs_pr1
	Array[File?] frip_macs2_qcs_pr2
	File? frip_macs2_qc_pooled
	File? frip_macs2_qc_ppr1 
	File? frip_macs2_qc_ppr2 
	Array[File]? frip_idr_qcs
	Array[File?] frip_idr_qcs_pr
	File? frip_idr_qc_ppr 
	Array[File?] frip_overlap_qcs
	Array[File?] frip_overlap_qcs_pr
	File? frip_overlap_qc_ppr
	File? idr_reproducibility_qc
	File? overlap_reproducibility_qc
	Array[File?] ataqc_txts
	Array[File?] ataqc_htmls

	File? qc_json_ref

	command {
		python /opt/src/encode_qc_report.py \
			${"--pipeline-ver " + pipeline_ver} \
			${"--title '" + sub(title,"'","_") + "'"} \
			${"--desc '" + sub(description,"'","_") + "'"} \
			${"--genome " + genome} \
			${"--multimapping " + multimapping} \
			--paired-ends ${sep=" " paired_ends} \
			--pipeline-type ${pipeline_type} \
			--peak-caller ${peak_caller} \
			${"--macs2-cap-num-peak " + macs2_cap_num_peak} \
			${"--spp-cap-num-peak " + spp_cap_num_peak} \
			--idr-thresh ${idr_thresh} \
			--flagstat-qcs ${sep="_:_" flagstat_qcs} \
			--nodup-flagstat-qcs ${sep="_:_" nodup_flagstat_qcs} \
			--dup-qcs ${sep="_:_" dup_qcs} \
			--pbc-qcs ${sep="_:_" pbc_qcs} \
			--xcor-plots ${sep="_:_" xcor_plots} \
			--xcor-scores ${sep="_:_" xcor_scores} \
			--idr-plots ${sep="_:_" idr_plots} \
			--idr-plots-pr ${sep="_:_" idr_plots_pr} \
			${"--idr-plot-ppr " + idr_plot_ppr} \
			--frip-macs2-qcs ${sep="_:_" frip_macs2_qcs} \
			--frip-macs2-qcs-pr1 ${sep="_:_" frip_macs2_qcs_pr1} \
			--frip-macs2-qcs-pr2 ${sep="_:_" frip_macs2_qcs_pr2} \
			${"--frip-macs2-qc-pooled " + frip_macs2_qc_pooled} \
			${"--frip-macs2-qc-ppr1 " + frip_macs2_qc_ppr1} \
			${"--frip-macs2-qc-ppr2 " + frip_macs2_qc_ppr2} \
			--frip-idr-qcs ${sep="_:_" frip_idr_qcs} \
			--frip-idr-qcs-pr ${sep="_:_" frip_idr_qcs_pr} \
			${"--frip-idr-qc-ppr " + frip_idr_qc_ppr} \
			--frip-overlap-qcs ${sep="_:_" frip_overlap_qcs} \
			--frip-overlap-qcs-pr ${sep="_:_" frip_overlap_qcs_pr} \
			${"--frip-overlap-qc-ppr " + frip_overlap_qc_ppr} \
			${"--idr-reproducibility-qc " + idr_reproducibility_qc} \
			${"--overlap-reproducibility-qc " + overlap_reproducibility_qc} \
			--ataqc-txts ${sep="_:_" ataqc_txts} \
			--ataqc-htmls ${sep="_:_" ataqc_htmls} \
			--out-qc-html qc.html \
			--out-qc-json qc.json \
			${"--qc-json-ref " + qc_json_ref}
	}
	output {
		File report = glob('*qc.html')[0]
		File qc_json = glob('*qc.json')[0]
		Boolean qc_json_ref_match = read_string("qc_json_ref_match.txt")=="True"
	}
	runtime {
		docker:"encode_atac/ycl:v1"
		cpu:"10"
		memory:"20G"
	}
}
