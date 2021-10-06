version 1.0

import "structs/DNASeqStructs.wdl" as Structs
import "Qc.wdl" as QC
import "BamProcessing.wdl" as Processing
import "Utilities.wdl" as Utils

workflow Alignment {

  input {
    Input inp
    DNASeqSingleSampleReferences references
    PapiSettings papi_settings

    Boolean check_contamination = true
    Boolean check_fingerprints = true

    File contamination_sites_ud
    File contamination_sites_bed
    File contamination_sites_mu

    String cross_check_fingerprints_by
    File haplotype_database_file
    Float lod_threshold
    Boolean hard_clip_reads = false
    
#    Boolean to_cram = false
    String? subset_region

    Boolean bin_base_qualities = true
    Boolean somatic = false
    Boolean dragmap = false
  }

  Float cutoff_for_large_rg_in_gb = 7.5
  Int compression_level = 2

  # We don't partition input BAM for distributed alignment, to avoid extra copying
  # and extra deduplication+merge step with a large overhead.
  # However, if you want shared parallelilsm, consider adding 
  # bazam sharded parallelism like in this pipeline:
  # https://github.com/Oshlack/STRetch/blob/c5345e5dea4adfde790befb9903ec2d81ed5b2c1/pipelines/pipeline_stages.groovy#L101

  if (inp.bam_or_cram_or_fastq1 == sub(inp.bam_or_cram_or_fastq1, ".cram$", "") + ".cram" ||
    inp.bam_or_cram_or_fastq1 == sub(inp.bam_or_cram_or_fastq1, ".bam$", "") + ".bam") {
  
    call BwaBamsormadupFromBamOrCram as BamCramToBam {
      input:
        bam_or_cram = inp.bam_or_cram_or_fastq1,
        bai_or_crai = inp.bai_or_crai_or_fastq2,
        RGID = inp.RGID,
        RGPL = inp.RGPL,
        RGPU = inp.RGPU,
        RGLB = inp.RGLB,
        RGCN = inp.RGCN,
        sample_name = inp.sample_name,
        output_bam_basename = inp.base_file_name,
        reference_fasta = references.reference_fasta,
        preemptible_tries = papi_settings.preemptible_tries,
        duplicate_metrics_fname = inp.base_file_name + ".duplicate_metrics",
        # to_cram = to_cram,
        subset_region = subset_region,
    }
  }

  if (inp.bam_or_cram_or_fastq1 != sub(inp.bam_or_cram_or_fastq1, ".cram$", "") + ".cram" &&
    inp.bam_or_cram_or_fastq1 != sub(inp.bam_or_cram_or_fastq1, ".bam$", "") + ".bam") {

    if(dragmap) {
      call DragmapAndBamsormadup as DragmapFastqToBam {
        input:
          fastq1 = inp.bam_or_cram_or_fastq1,
          fastq2 = inp.bai_or_crai_or_fastq2,
          RGID = inp.RGID,
          RGPL = inp.RGPL,
          RGPU = inp.RGPU,
          RGLB = inp.RGLB,
          RGCN = inp.RGCN,
          sample_name = inp.sample_name,
          output_bam_basename = inp.base_file_name,
          reference_fasta = references.reference_fasta,
          dragmap_reference_dir = references.dragmap_reference_dir,
          preemptible_tries = papi_settings.preemptible_tries,
          duplicate_metrics_fname = inp.base_file_name + ".duplicate_metrics"
      }


      call Dragmap {
        input:
          fastq1 = inp.bam_or_cram_or_fastq1,
          fastq2 = inp.bai_or_crai_or_fastq2,
          RGID = inp.RGID,
          sample_name = inp.sample_name,
          output_bam_basename = inp.base_file_name,
          dragmap_reference_dir = references.dragmap_reference_dir,
          preemptible_tries = papi_settings.preemptible_tries
      }

      call Processing.RenameHeader as RenameHeader {
        input:
          sortedbam = Dragmap.output_file,
          sample_name = inp.sample_name,
          RGID = inp.RGID,
          RGPL = inp.RGPL,
          RGPU = inp.RGPU,
          RGLB = inp.RGLB,
          RGCN = inp.RGCN,
          output_bam_basename = inp.base_file_name,

          preemptible_tries = papi_settings.preemptible_tries
      }

      Float total_bam_size = size(RenameHeader.output_file, "GiB")

      call Processing.MarkDuplicates as MarkDuplicates {
        input:
          input_bams = [RenameHeader.output_file],
          output_bam_basename = inp.base_file_name,
          metrics_filename = inp.base_file_name + ".duplicate_metrics",
          total_input_size = total_bam_size,
          compression_level = 2,
          preemptible_tries = papi_settings.preemptible_tries
      }
    }
    if(!dragmap){
      call BwaAndBamsormadup as FastqToBam {
        input:
          fastq1 = inp.bam_or_cram_or_fastq1,
          fastq2 = inp.bai_or_crai_or_fastq2,
          RGID = inp.RGID,
          RGPL = inp.RGPL,
          RGPU = inp.RGPU,
          RGLB = inp.RGLB,
          RGCN = inp.RGCN,
          sample_name = inp.sample_name,
          output_bam_basename = inp.base_file_name,
          reference_fasta = references.reference_fasta,
          preemptible_tries = papi_settings.preemptible_tries,
          duplicate_metrics_fname = inp.base_file_name + ".duplicate_metrics"
          # to_cram = to_cram
      }
    }
  }
  
  File mapped_file = select_first([BamCramToBam.output_file, FastqToBam.output_file, DragmapFastqToBam.output_file])
  File mapped_indx = select_first([BamCramToBam.output_indx, FastqToBam.output_indx, DragmapFastqToBam.output_indx])
  Float mapped_file_size = size(mapped_file, "GiB")

  # Run BQSR here

  # Create list of sequences for scatter-gather parallelization
  call Utils.CreateSequenceGroupingTSV as CreateSequenceGroupingTSV {
    input:
      ref_dict = references.reference_fasta.ref_dict,
      preemptible_tries = papi_settings.preemptible_tries
  }

  Int num_of_bqsr_scatters = length(CreateSequenceGroupingTSV.sequence_grouping)
  Int potential_bqsr_divisor = num_of_bqsr_scatters - 10
  Int bqsr_divisor = if potential_bqsr_divisor > 1 then potential_bqsr_divisor else 1

  # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel
  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
    # Generate the recalibration model by interval
    call Processing.BaseRecalibrator as BaseRecalibrator {
      input:
        input_bam = mapped_file,
        input_bam_index = mapped_indx,
        recalibration_report_filename = inp.base_file_name + ".recal_data.csv",
        sequence_group_interval = subgroup,
        dbsnp_vcf = references.dbsnp_vcf,
        dbsnp_vcf_index = references.dbsnp_vcf_index,
        known_indels_sites_vcfs = references.known_indels_sites_vcfs,
        known_indels_sites_indices = references.known_indels_sites_indices,
        ref_dict = references.reference_fasta.ref_dict,
        ref_fasta = references.reference_fasta.ref_fasta,
        ref_fasta_index = references.reference_fasta.ref_fasta_index,
        bqsr_scatter = bqsr_divisor,
        preemptible_tries = papi_settings.agg_preemptible_tries
    }
  }

  # Merge the recalibration reports resulting from by-interval recalibration
  # The reports are always the same size
  call Processing.GatherBqsrReports as GatherBqsrReports {
    input:
      input_bqsr_reports = BaseRecalibrator.recalibration_report,
      output_report_filename = inp.base_file_name + ".recal_data.csv",
      preemptible_tries = papi_settings.preemptible_tries
  }

  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping_with_unmapped) {
    # Apply the recalibration model by interval
    call Processing.ApplyBQSR as ApplyBQSR {
      input:
        input_bam = mapped_file,
        input_bam_index = mapped_indx,
        output_bam_basename = inp.base_file_name + ".recalibrated",
        recalibration_report = GatherBqsrReports.output_bqsr_report,
        sequence_group_interval = subgroup,
        ref_dict = references.reference_fasta.ref_dict,
        ref_fasta = references.reference_fasta.ref_fasta,
        ref_fasta_index = references.reference_fasta.ref_fasta_index,
        bqsr_scatter = bqsr_divisor,
        compression_level = compression_level,
        preemptible_tries = papi_settings.agg_preemptible_tries,
        bin_base_qualities = bin_base_qualities,
        somatic = somatic
    }
  }

  # Merge the recalibrated BAM files resulting from by-interval recalibration
  call Processing.GatherSortedBamFiles as GatherBamFiles {
    input:
      input_bams = ApplyBQSR.recalibrated_bam,
      output_bam_basename = inp.base_file_name,
      total_input_size = mapped_file_size,
      compression_level = compression_level,
      preemptible_tries = papi_settings.agg_preemptible_tries
  }

  # Convert to CRAM here
#  if(to_cram){
#    call Utils.ConvertToCram as Crammer {
#      input:
#        input_bam = GatherBamFiles.output_bam,
#        ref_fasta = references.reference_fasta.ref_fasta,
#        ref_fasta_index = references.reference_fasta.ref_fasta_index,
#        output_basename = inp.base_file_name,
#        preemptible_tries = papi_settings.agg_preemptible_tries
#    }
#  }

  if (defined(haplotype_database_file) && check_fingerprints) {
    # Check identity of fingerprints across readgroups
    call QC.CrossCheckFingerprints as CrossCheckFingerprints {
      input:
        input_bams = [ mapped_file ],
        input_bam_indexes = [ mapped_indx ],
        haplotype_database_file = haplotype_database_file,
        metrics_filename = inp.base_file_name + ".crosscheck",
        total_input_size = mapped_file_size,
        lod_threshold = lod_threshold,
        cross_check_by = cross_check_fingerprints_by,
        preemptible_tries = papi_settings.agg_preemptible_tries,
        ref_dict = references.reference_fasta.ref_dict,
        ref_fasta = references.reference_fasta.ref_fasta,
        ref_fasta_index = references.reference_fasta.ref_fasta_index
    }
  }

  if (check_contamination) {
    # Estimate level of cross-sample contamination
    call Processing.CheckContamination as CheckContamination {
      input:
        input_bam = mapped_file,
        input_bam_index = mapped_indx,
        contamination_sites_ud = contamination_sites_ud,
        contamination_sites_bed = contamination_sites_bed,
        contamination_sites_mu = contamination_sites_mu,
        ref_fasta = references.reference_fasta.ref_fasta,
        ref_fasta_index = references.reference_fasta.ref_fasta_index,
        output_prefix = inp.base_file_name,
        preemptible_tries = papi_settings.agg_preemptible_tries,
        contamination_underestimation_factor = 0.75
    }
  }

  # Outputs that will be retained when execution is complete
  output {
    File? cross_check_fingerprints_metrics = CrossCheckFingerprints.cross_check_fingerprints_metrics

    File? selfSM = CheckContamination.selfSM
    Float? contamination = CheckContamination.contamination

    File duplicate_metrics = select_first([BamCramToBam.duplicate_metrics, FastqToBam.duplicate_metrics, DragmapFastqToBam.duplicate_metrics])

    File output_file = GatherBamFiles.output_bam
#    File output_file = mapped_file
    File output_indx = GatherBamFiles.output_bam_index
#    File output_indx = mapped_indx
#    File? output_cram = Crammer.output_cram
#    File? output_cram_index = Crammer.output_cram_index
#    File? output_cram_md5 = Crammer.output_cram_md5
  }
  meta {
    allowNestedInputs: true
  }
}

workflow Crammer {

  input {
    File input_bam
    String base_file_name
    DNASeqSingleSampleReferences references
    PapiSettings papi_settings
  }

  call Utils.ConvertToCram as Crammer {
    input:
      input_bam = input_bam,
      ref_fasta = references.reference_fasta.ref_fasta,
      ref_fasta_index = references.reference_fasta.ref_fasta_index,
      output_basename = base_file_name,
      preemptible_tries = papi_settings.agg_preemptible_tries
    }


  output {
    File output_cram = Crammer.output_cram
    File output_cram_index = Crammer.output_cram_index
    File output_cram_md5 = Crammer.output_cram_md5
  }

  meta {
    allowNestedInputs: true
  }
}

task Dragmap {
  input {
    File fastq1
    File fastq2
    String RGID
    String sample_name
    String output_bam_basename

    DragmapReferenceDir dragmap_reference_dir

    Int preemptible_tries
  }
  Int total_cpu = 32
  String output_file = "~{output_bam_basename}.bam"

  command <<<
    set -o pipefail
    set -ex

    dragen-os -r `dirname ~{dragmap_reference_dir.ht_cfg}` -1 ~{fastq1} -2 ~{fastq2} \
      --RGID ~{RGID} --RGSM ~{sample_name} \
      --num-threads ~{total_cpu} 2> >(tee ~{output_bam_basename}.dragmap.stderr.log >&2) | \
    samtools sort -Obam -o ~{output_file}
  >>>
  runtime {
    docker: "shyrav/dragmap-biobambam2-samtools:0.0"
    preemptible: preemptible_tries
    memory: "120 GiB"
    cpu: total_cpu
    disks: "local-disk " + 300 + " HDD"
  }
  output {
    File output_file = output_file
    File dragmap_stderr_log = "~{output_bam_basename}.dragmap.stderr.log"
  }
}

task DragmapAndBamsormadup {
  input {
    File fastq1
    File fastq2
    String sample_name
    String RGID
    String RGPU
    String RGPL
    String RGCN
    String RGLB
    String output_bam_basename
    String duplicate_metrics_fname

    ReferenceFasta reference_fasta
    DragmapReferenceDir dragmap_reference_dir

    Int preemptible_tries
  }

  Int total_cpu = 32
  String rg_line = "@RG\\tID:~{RGID}\\tSM:~{sample_name}\\tPL:~{RGPL}\\tPU:~{RGPU}\\tLB:~{RGLB}\\tCN:~{RGCN}"
  String output_file = "~{output_bam_basename}.bam"

  command <<<
    set -o pipefail
    set -ex

    (while true; do df -h; pwd; du -sh *; free -m; sleep 300; done) &

    # Align and Markdups
    dragen-os -r `dirname ~{dragmap_reference_dir.ht_cfg}` -1 ~{fastq1} -2 ~{fastq2} \
      --RGID ~{RGID} --RGSM ~{sample_name} \
      --num-threads ~{total_cpu} 2> >(tee ~{output_bam_basename}.dragmap.stderr.log >&2) | \
    bamsormadup threads=~{total_cpu} SO=coordinate inputformat=sam outputformat=bam \
      reference=~{reference_fasta.ref_fasta} \
      M=~{duplicate_metrics_fname} O=tmp.bam > tmp.bam

    # Remove FASTQs to make space to reheader bam file
    rm -f ~{fastq1} ~{fastq2}

    # Changing header of bam to best practices version.
    samtools view -H tmp.bam > header.sam
    oldline=`grep "^@RG" header.sam`
    newline=`echo -e "~{rg_line}"`
    sed -i "s/$oldline/$newline/" header.sam
    samtools reheader header.sam tmp.bam > ~{output_file}
    rm tmp.bam
    samtools index ~{output_file}

    df -h; pwd; du -sh *
  >>>
  runtime {
    docker: "shyrav/dragmap-biobambam2-samtools:0.0"
    preemptible: preemptible_tries
    memory: "120 GiB"
    cpu: total_cpu
    disks: "local-disk " + 400 + " HDD"
  }
  output {
    File output_file = output_file
    File output_indx = "~{output_file}.bai"
    File dragmap_stderr_log = "~{output_bam_basename}.dragmap.stderr.log"
    File duplicate_metrics = "~{duplicate_metrics_fname}"
  }
}

task BwaAndBamsormadup {
  input {
    File fastq1
    File fastq2
    String sample_name
    String RGID
    String RGPU
    String RGPL
    String RGCN
    String RGLB
    String output_bam_basename
    String duplicate_metrics_fname

    # reference_fasta.ref_alt is the .alt file from bwa-kit
    # (https://github.com/lh3/bwa/tree/master/bwakit),
    # listing the reference contigs that are "alternative".
    ReferenceFasta reference_fasta

    Int preemptible_tries
#    Boolean to_cram = false
  }

  #  Int bwa_cpu = 16
  #  Int bamsormadup_cpu = 16
  Int total_cpu = 16

  String rg_line = "@RG\\tID:~{RGID}\\tSM:~{sample_name}\\tPL:~{RGPL}\\tPU:~{RGPU}\\tLB:~{RGLB}\\tCN:~{RGCN}"

  String output_file = "~{output_bam_basename}.bam"
#  String output_file = if to_cram then "~{output_bam_basename}.cram" else "~{output_bam_basename}.bam"
  String output_indx = "~{output_bam_basename}.bai"
#  String output_indx = if to_cram then "~{output_bam_basename}.crai" else "~{output_bam_basename}.bai"

  # BWA command options:
  # -K     process INT input bases in each batch regardless of nThreads (for reproducibility)
  # -t16   threads
  # -R     read group header line such as '@RG\tID:foo\tSM:bar'
  command <<<
    set -o pipefail
    set -ex

    (while true; do df -h; pwd; du -sh *; free -m; sleep 300; done) &

    bwa mem \
      -t 16 \
      -R '~{rg_line}' \
      -K 100000000 \
      -v 3 -Y ~{reference_fasta.ref_fasta} ~{fastq1} ~{fastq2} \
      2> >(tee ~{output_bam_basename}.bwa.stderr.log >&2) | \
    bamsormadup \
      threads=16 \
      inputformat=sam \
      outputformat=bam \
      reference=~{reference_fasta.ref_fasta} \
      M=~{duplicate_metrics_fname} \
      indexfilename=~{output_indx} \
      optminpixeldif=2500 \
      O=~{output_file} > ~{output_file}

    df -h; pwd; du -sh *
  >>>

  runtime {
    # docker: "australia-southeast1-docker.pkg.dev/cpg-common/images/bazam:v2"
    # cromwell doesn't work with artifact registry:
    # java.lang.Exception: Registry australia-southeast1-docker.pkg.dev is not supported
    # docker: "gcr.io/cpg-common/bwa-bazam:v1"
    docker: "gcr.io/pb-dev-312200/biobambam2-samtools-picard-bwa:latest"
    preemptible: preemptible_tries
    memory: "64 GiB"
    cpu: total_cpu
    disks: "local-disk " + 400 + " HDD"
  }
  output {
    File output_file = output_file
    File output_indx = output_indx
    File bwa_stderr_log = "~{output_bam_basename}.bwa.stderr.log"
    File duplicate_metrics = "~{duplicate_metrics_fname}"
  }
}
  
task BwaFromFastq {
  input {
    File fastq1
    File fastq2
    String RGID
    String RGPU
    String RGPL
    String RGCN
    String RGLB
    String sample_name
    String output_bam_basename
    String duplicate_metrics_fname

    # reference_fasta.ref_alt is the .alt file from bwa-kit
    # (https://github.com/lh3/bwa/tree/master/bwakit),
    # listing the reference contigs that are "alternative".
    ReferenceFasta reference_fasta

    Int preemptible_tries
    # Boolean to_cram = false
  }
  
#  String output_format = if to_cram then "cram" else "bam"
  String output_format = "bam"

  Int bwa_cpu = 25
  Int bamsormadup_cpu = 6
  Int total_cpu = bwa_cpu + bamsormadup_cpu

  String rg_line = "@RG\\tID:~{RGID}\\tSM:~{sample_name}\\tPL:~{RGPL}\\tPU:~{RGPU}\\tLB:~{RGLB}\\tCN:~{RGCN}"
  
  String output_file = "~{output_bam_basename}.bam"
#  String output_file = if to_cram then "~{output_bam_basename}.cram" else "~{output_bam_basename}.bam"
  String output_indx = "~{output_bam_basename}.bai"
#  String output_indx = if to_cram then "~{output_bam_basename}.crai" else "~{output_bam_basename}.bai"

  # BWA command options:
  # -K     process INT input bases in each batch regardless of nThreads (for reproducibility)
  # -t16   threads
  # -R     read group header line such as '@RG\tID:foo\tSM:bar'
  command <<<
    set -o pipefail
    set -ex

    (while true; do df -h; pwd; du -sh *; free -m; sleep 300; done) &
    
    bwa mem -K 100000000 -t~{bwa_cpu} -R '~{rg_line}' \
      ~{reference_fasta.ref_fasta} ~{fastq1} ~{fastq2} \
      2> >(tee ~{output_bam_basename}.bwa.stderr.log >&2) | \
    bamsormadup inputformat=sam threads=~{bamsormadup_cpu} SO=coordinate \
      M=~{duplicate_metrics_fname} \
      outputformat=sam | \
    samtools view -T ~{reference_fasta.ref_fasta} \
      -O ~{output_format} \
      -o ~{output_file}
    
    samtools index -@~{total_cpu} ~{output_file} ~{output_indx}

    df -h; pwd; du -sh *
  >>>
  
  runtime {
    # docker: "australia-southeast1-docker.pkg.dev/cpg-common/images/bazam:v2"
    # cromwell doesn't work with artifact registry:
    # java.lang.Exception: Registry australia-southeast1-docker.pkg.dev is not supported
    # docker: "gcr.io/cpg-common/bwa-bazam:v1"
    docker: "gcr.io/pb-dev-312200/biobambam2-samtools-picard-bwa:latest"
    preemptible: preemptible_tries
    memory: "64 GiB"
    cpu: total_cpu
    disks: "local-disk " + 300 + " HDD"
  }
  output {
    File output_file = output_file
    File output_indx = output_indx
    File bwa_stderr_log = "~{output_bam_basename}.bwa.stderr.log"
    File duplicate_metrics = "~{duplicate_metrics_fname}"
  }
}

task BwaBamsormadupFromBamOrCram {
  input {
    File bam_or_cram
    File bai_or_crai
    String RGID
    String RGPU
    String RGPL
    String RGCN
    String RGLB
    String sample_name
    String output_bam_basename
    String duplicate_metrics_fname

    # reference_fasta.ref_alt is the .alt file from bwa-kit
    # (https://github.com/lh3/bwa/tree/master/bwakit),
    # listing the reference contigs that are "alternative".
    ReferenceFasta reference_fasta

    Int preemptible_tries
    # Boolean to_cram = false
    String? subset_region
  }

#  String output_format = if to_cram then "cram" else "bam"
  String output_format = "bam"
  String bazam_regions = if defined(subset_region) then "--regions ~{subset_region} " else ""

  # Not been tested yet
  Int bwa_cpu = 16
  Int bazam_cpu = 16
  Int bamsormadup_cpu = 16
  Int total_cpu = 16

  String rg_line = "@RG\\tID:~{RGID}\\tSM:~{sample_name}\\tPL:~{RGPL}\\tPU:~{RGPU}\\tLB:~{RGLB}\\tCN:~{RGCN}"

  String output_file = "~{output_bam_basename}.bam"
#  String output_file = if to_cram then "~{output_bam_basename}.cram" else "~{output_bam_basename}.bam"
#  String output_indx = if to_cram then "~{output_bam_basename}.crai" else "~{output_bam_basename}.bai"
  String output_indx = "~{output_bam_basename}.bai"

  # BWA command options:
  # -K     process INT input bases in each batch regardless of nThreads (for reproducibility)
  # -p     smart pairing (ignoring in2.fq)
  # -t16   threads
  # -R     read group header line such as '@RG\tID:foo\tSM:bar'
  command <<<
    set -o pipefail
    set -ex

    (while true; do df -h; pwd; du -sh *; free -m; sleep 300; done) &

    bazam -Xmx16g -Dsamjdk.reference_fasta=~{reference_fasta.ref_fasta} \
      ~{bazam_regions} -n~{bazam_cpu} -bam ~{bam_or_cram} | \
    bwa mem -K 100000000 -p -t~{bwa_cpu} -R '~{rg_line}' \
      ~{reference_fasta.ref_fasta} /dev/stdin - \
      2> >(tee ~{output_bam_basename}.bwa.stderr.log >&2) | \
    bamsormadup \
      threads=16 \
      inputformat=sam \
      outputformat=bam \
      reference=~{reference_fasta.ref_fasta} \
      M=~{duplicate_metrics_fname} \
      indexfilename=~{output_indx} \
      optminpixeldif=2500 \
      O=~{output_file} > ~{output_file}

    df -h; pwd; du -sh *
  >>>

  runtime {
    # docker: "australia-southeast1-docker.pkg.dev/cpg-common/images/bazam:v2"
    # cromwell doesn't work with artifact registry:
    # java.lang.Exception: Registry australia-southeast1-docker.pkg.dev is not supported
    # docker: "gcr.io/cpg-common/bwa-bazam:v1"
    docker: "gcr.io/pb-dev-312200/biobambam2-samtools-picard-bwa:latest"
    preemptible: preemptible_tries
    memory: "64 GiB"
    cpu: total_cpu
    disks: "local-disk " + 300 + " HDD"
  }
  output {
    File output_file = output_file
    File output_indx = output_indx
    File bwa_stderr_log = "~{output_bam_basename}.bwa.stderr.log"
    File duplicate_metrics = "~{duplicate_metrics_fname}"
  }
}

task BwaFromBamOrCram {
  input {
    File bam_or_cram
    File bai_or_crai
    String RGID
    String RGPU
    String RGPL
    String RGCN
    String RGLB
    String sample_name
    String output_bam_basename
    String duplicate_metrics_fname

    # reference_fasta.ref_alt is the .alt file from bwa-kit
    # (https://github.com/lh3/bwa/tree/master/bwakit),
    # listing the reference contigs that are "alternative".
    ReferenceFasta reference_fasta

    Int preemptible_tries
#    Boolean to_cram = false
    String? subset_region 
  }
  
#  String output_format = if to_cram then "cram" else "bam"
  String output_format = "bam"
  String bazam_regions = if defined(subset_region) then "--regions ~{subset_region} " else ""
  
  Int bwa_cpu = 20
  Int bazam_cpu = 5
  Int bamsormadup_cpu = 6
  Int total_cpu = bwa_cpu + bazam_cpu + bamsormadup_cpu

  String rg_line = "@RG\\tID:~{RGID}\\tSM:~{sample_name}\\tPL:~{RGPL}\\tPU:~{RGPU}\\tLB:~{RGLB}\\tCN:~{RGCN}"
  
  String output_file = "~{output_bam_basename}.bam"
#  String output_file = if to_cram then "~{output_bam_basename}.cram" else "~{output_bam_basename}.bam"
#  String output_indx = if to_cram then "~{output_bam_basename}.crai" else "~{output_bam_basename}.bai"
  String output_indx = "~{output_bam_basename}.bai"

  # BWA command options:
  # -K     process INT input bases in each batch regardless of nThreads (for reproducibility)
  # -p     smart pairing (ignoring in2.fq)
  # -t16   threads
  # -R     read group header line such as '@RG\tID:foo\tSM:bar'
  command <<<
    set -o pipefail
    set -ex

    (while true; do df -h; pwd; du -sh *; free -m; sleep 300; done) &
    
    bazam -Xmx16g -Dsamjdk.reference_fasta=~{reference_fasta.ref_fasta} \
      ~{bazam_regions} -n~{bazam_cpu} -bam ~{bam_or_cram} | \
    bwa mem -K 100000000 -p -t~{bwa_cpu} -R '~{rg_line}' \
      ~{reference_fasta.ref_fasta} /dev/stdin - \
      2> >(tee ~{output_bam_basename}.bwa.stderr.log >&2) | \
    bamsormadup inputformat=sam threads=~{bamsormadup_cpu} SO=coordinate \
      M=~{duplicate_metrics_fname} \
      outputformat=sam | \
    samtools view -T ~{reference_fasta.ref_fasta} \
      -O ~{output_format} \
      -o ~{output_file}
    
    samtools index -@~{total_cpu} ~{output_file} ~{output_indx}

    df -h; pwd; du -sh *
  >>>
  
  runtime {
    # docker: "australia-southeast1-docker.pkg.dev/cpg-common/images/bazam:v2"
    # cromwell doesn't work with artifact registry:
    # java.lang.Exception: Registry australia-southeast1-docker.pkg.dev is not supported
    # docker: "gcr.io/cpg-common/bwa-bazam:v1"
    docker: "gcr.io/pb-dev-312200/biobambam2-samtools-picard-bwa:latest"
    preemptible: preemptible_tries
    memory: "64 GiB"
    cpu: total_cpu
    disks: "local-disk " + 300 + " HDD"
  }
  output {
    File output_file = output_file
    File output_indx = output_indx
    File bwa_stderr_log = "~{output_bam_basename}.bwa.stderr.log"
    File duplicate_metrics = "~{duplicate_metrics_fname}"
  }
}

task SamSplitter {
  input {
    File input_bam
    Int n_reads
    Int preemptible_tries
    Int compression_level
  }

  Float unmapped_bam_size = size(input_bam, "GiB")
  # Since the output bams are less compressed than the input bam we need a disk multiplier that's larger than 2.
  Float disk_multiplier = 2.5
  Int disk_size = ceil(disk_multiplier * unmapped_bam_size + 20)

  command {
    set -e
    mkdir output_dir

    total_reads=$(samtools view -c ~{input_bam})

    java -Dsamjdk.compression_level=~{compression_level} -Xms3000m -jar /usr/gitc/picard.jar SplitSamByNumberOfReads \
      INPUT=~{input_bam} \
      OUTPUT=output_dir \
      SPLIT_TO_N_READS=~{n_reads} \
      TOTAL_READS_IN_INPUT=$total_reads
  }
  output {
    Array[File] split_bams = glob("output_dir/*.bam")
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
    preemptible: preemptible_tries
    memory: "3.75 GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
}
