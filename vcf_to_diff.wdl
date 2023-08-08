version 1.0

task make_mask_and_diff {
	# This combines the creation of the bed graph histogram mask file and the
	# creation of the diff file into one WDL task. Sometimes, in WDL, it is
	# easier to combine tasks to avoid shenanigans with the Array[File] type.
	# The masking part of this task is based on github.com/aofarrel/mask-by-coverage
	input {
		File bam
		File vcf
		File? tbmf
		Int min_coverage_per_site
		Boolean histograms = false
		Float discard_sample_if_more_than_this_percent_is_low_coverage # for sra, default was 0.05

		# runtime attributes
		Int addldisk = 10
		Int cpu      = 8
		Int retries  = 1
		Int memory   = 16
		Int preempt  = 1
	}
	String basename_bam = basename(bam, ".bam")
	String basename_vcf = basename(vcf, ".vcf")
	Int finalDiskSize = ceil(size(bam, "GB")*2) + ceil(size(vcf, "GB")*2) + addldisk
	
	command <<<
	set -eux pipefail
	start=$(date +%s)

	# We want the mask file the user input, if it exists, to be the mask file, else
	# fall back on a default mask file that exists in the Docker image already.
	# We cannot use WDL built-in select_first, or else the user inputting a mask file
	# will result in WDL looking for the literal gs:// URI rather than while the file
	# is localized. Different WDL executors localize files to different places, so the
	# following workaround, while goofy, seems to be the most robust.
	if [[ "~{tbmf}" = "" ]]
	then
		mask="/mask/R00000039_repregions.bed"
	else
		mask="~{tbmf}"
	fi
	
	echo "Copying bam..."
	cp ~{bam} .
	
	echo "Sorting bam..."
	samtools sort -u ~{basename_bam}.bam > sorted_u_~{basename_bam}.bam
	
	echo "Calculating coverage..."
	bedtools genomecov -ibam sorted_u_~{basename_bam}.bam -bga | \
		awk '$4 < ~{min_coverage_per_site}' > \
		~{basename_bam}_below_~{min_coverage_per_site}x_coverage.bedgraph
	
	if [[ "~{histograms}" = "true" ]]
	then
		echo "Generating histograms..."
		bedtools genomecov -ibam sorted_u_~{basename_bam}.bam > histogram.txt
	fi
	
	echo "Pulling diff script..."
	wget https://raw.githubusercontent.com/aofarrel/parsevcf/1.1.8/vcf_to_diff_script.py
	echo "Running script..."
	python3 vcf_to_diff_script.py -v ~{vcf} \
	-d . \
	-tbmf ${mask} \
	-bed ~{basename_bam}_below_~{min_coverage_per_site}x_coverage.bedgraph \
	-cd ~{min_coverage_per_site}
	
	# if the sample has too many low coverage sites, throw it out
	this_files_info=$(awk -v file_to_check="~{basename_vcf}.diff" '$1 == file_to_check' "~{basename_vcf}.report")
	if [[ ! "$this_files_info" = "" ]]
	then
		# okay, we have information about this file. is it above the removal threshold?
		echo "$this_files_info" > temp
		amount_low_coverage=$(cut -f2 temp)
		percent_low_coverage=$(echo "$amount_low_coverage"*100 | bc)
		echo "$percent_low_coverage percent of ~{basename_vcf} is below ~{min_coverage_per_site}."
		
		# piping an inequality to `bc` will return 0 if false, 1 if true
		is_bigger=$(echo "$amount_low_coverage>~{discard_sample_if_more_than_this_percent_is_low_coverage}" | bc)
		if [[ $is_bigger == 0 ]]
		then
			# amount of low coverage is BELOW the removal threshold and passes
			echo "PASS" >> ERROR
	
		else
			# amount of low coverage is ABOVE the removal threshold and fails
			echo FAILURE - "$percent_low_coverage" is above "~{discard_sample_if_more_than_this_percent_is_low_coverage}" cutoff
			echo VCF2DIFF_"$percent_low_coverage"_PCT_BELOW_"~{min_coverage_per_site}"x_COVERAGE >> ERROR
		fi
	fi
			
	end=$(date +%s)
	seconds=$(echo "$end - $start" | bc)
	minutes=$(echo "$seconds" / 60 | bc)
	echo "Finished in about $minutes minutes ($seconds sec))"
	ls -lha
	>>>

	runtime {
		cpu: cpu
		docker: "ashedpotatoes/sranwrp:1.1.12"
		disks: "local-disk " + finalDiskSize + " HDD"
		maxRetries: "${retries}"
		memory: "${memory} GB"
		preemptible: "${preempt}"
	}

	meta {
		author: "Lily Karim (WDLization by Ash O'Farrell)"
	}

	output {
		File mask_file = basename_bam+"_below_"+min_coverage_per_site+"x_coverage.bedgraph"  # !StringCoercion
		File? diff = basename_vcf+".diff"
		File? report = basename_vcf+".report"
		File? histogram = "histogram.txt"
		String errorcode = read_string("ERROR")
	}
}

task make_diff {
	input {
		File vcf
		File tbmf
		File cf
		Int cd = 10

		# runtime attributes
		Int addldisk = 10
		Int cpu	= 8
		Int retries	= 1
		Int memory = 16
		Int preempt	= 1
	}
	# estimate disk size
	String basename = basename(vcf)
	Int finalDiskSize = 2*ceil(size(vcf, "GB")) + addldisk

	command <<<
		set -eux pipefail
		mkdir outs
		wget https://raw.githubusercontent.com/lilymaryam/parsevcf/1.0.4/vcf_to_diff_script.py
		python3.10 vcf_to_diff_script.py -v ~{vcf} -d ./outs/ -tbmf ~{tbmf} -cf ~{cf} -cd ~{cd}
		ls -lha outs/
	>>>

	runtime {
		cpu: cpu
		disks: "local-disk " + finalDiskSize + " SSD"
		docker: "ashedpotatoes/sranwrp:1.1.6"
		maxRetries: "${retries}"
		memory: "${memory} GB"
		preemptible: "${preempt}"
	}

	meta {
		author: "Lily Karim (WDLization by Ash O'Farrell)"
	}

	output {
		File diff = "outs/"+basename+".diff"
		File report = "outs/"+basename+".report"
	}
}
