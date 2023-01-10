version 1.0

task make_mask_and_diff {
	# This combines the creation of the bed graph histogram mask file and the
	# creation of the diff file into one WDL task. Sometimes, in WDL, it is
	# easier to combine tasks to avoid shenanigans with the Array[File] type.
	# The masking part of this task is based on github.com/aofarrel/mask-by-coverage
	input {
		File bam
		File vcf
		File tbmf
		Int min_coverage

		# runtime attributes
		Int addldisk = 250
		Int cpu      = 16
		Int retries  = 1
		Int memory   = 32
		Int preempt  = 1
	}
	String basestem = basename(bam, ".bam")
	Int finalDiskSize = ceil(size(bam, "GB")) + addldisk
	
	command <<<
	set -eux pipefail
	cp ~{bam} .
	samtools sort -u ~{basestem}.bam > sorted_u_~{basestem}.bam
	bedtools genomecov -ibam sorted_u_~{basestem}.bam -bga | \
		awk '$4 < ~{min_coverage}' > \
		~{basestem}_below_~{min_coverage}x_coverage.bedgraph
	mkdir outs
	wget https://raw.githubusercontent.com/lilymaryam/parsevcf/fix_endstart/vcf_to_diff_script.py
	python3.10 vcf_to_diff_script.py -v ~{vcf} -d ./outs/ -tbmf ~{tbmf} -cf ~{basestem}_below_~{min_coverage}x_coverage.bedgraph -cd ~{min_coverage}
	>>>

	runtime {
		cpu: cpu
		docker: "ashedpotatoes/sranwrp:1.0.7"
		disks: "local-disk " + finalDiskSize + " HDD"
		maxRetries: "${retries}"
		memory: "${memory} GB"
		preemptible: "${preempt}"
	}

	meta {
		author: "Lily Karim (WDLization by Ash O'Farrell)"
	}

	output {
		File diff = glob("outs/*.diff")[0]
		File debug_script = "vcf_to_diff_script.py" # to keep track of what's on main
		File mask_file = glob("*coverage.bedgraph")[0]
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
	Int finalDiskSize = 2*ceil(size(vcf, "GB")) + addldisk

	command <<<
		set -eux pipefail
		mkdir outs
		wget https://raw.githubusercontent.com/lilymaryam/parsevcf/main/vcf_to_diff_script.py
		python3.10 vcf_to_diff_script.py -v ~{vcf} -d ./outs/ -tbmf ~{tbmf} -cf ~{cf} -cd ~{cd}
	>>>

	runtime {
		cpu: cpu
		disks: "local-disk " + finalDiskSize + " SSD"
		docker: "ashedpotatoes/sranwrp:1.0.7"
		maxRetries: "${retries}"
		memory: "${memory} GB"
		preemptible: "${preempt}"
	}

	meta {
		author: "Lily Karim (WDLization by Ash O'Farrell)"
	}

	output {
		File diff = glob("outs/*.diff")[0]
		File debug_script = "vcf_to_diff_script.py" # to keep track of what's on main
	}
}
