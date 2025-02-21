
// Set here a list of files or directories to use. E.g. Channel.fromPath(["/path/*", "/path2/file"], type: 'any')
INPUT = Channel.fromPath(params.input.tokenize())
OUTPUT = file(params.output)
STEPS_FOLDER = file(params.stepsFolder)
ANNOTATIONS = Channel.value(params.annotations)
REGIONS = Channel.value("${params.datasets}/regions/cds.regions.gz")



process ParseInput {
	tag "Parse input ${input}"
	label "core"
	publishDir "${STEPS_FOLDER}/inputs", mode: "copy"
	errorStrategy 'finish'

	input:
		path input from INPUT
		path annotations from ANNOTATIONS

	output:
		path("*.parsed.tsv.gz") into COHORTS

	script:
		if ( input.toRealPath().toFile().isDirectory() || input.endsWith(".bginfo" )) {
			"""
			openvar groupby ${input} --header -g DATASET -q -s 'gzip > \${GROUP_KEY}.parsed.tsv.gz'
			"""
		}
		else {
			cohort = input.baseName.split('\\.')[0]
			"""
			openvar cat ${input} --header | gzip > ${cohort}.parsed.tsv.gz
			"""
		}
}


COHORTS
	.flatten()
	.map{it -> [it.baseName.split('\\.')[0], it]}
	.into{ COHORTS1; COHORTS2; COHORTS3; COHORTS4; COHORTS5 }

process LoadCancer {
	tag "Load cancer type ${cohort}"
	label "core"

	input:
		tuple val(cohort), path(input) from COHORTS1

	output:
		tuple val(cohort), stdout into CANCERS

	script:
		"""
		get_field.sh ${input} CANCER
		"""
}

CANCERS.into { CANCERS1; CANCERS2; CANCERS3 }

process LoadPlatform {
	tag "Load sequencing platform ${cohort}"
	label "core"

	input:
		tuple val(cohort), path(input) from COHORTS2

	output:
		tuple val(cohort), stdout into PLATFORMS

	script:
		"""
		get_field.sh ${input} PLATFORM
		"""
}

PLATFORMS.into { PLATFORMS1; PLATFORMS2; PLATFORMS3 }

process LoadGenome {
	tag "Load reference genome ${cohort}"
	label "core"

	input:
		tuple val(cohort), path(input) from COHORTS3

	output:
		tuple val(cohort), stdout into GENOMES

	script:
		"""
		get_field.sh ${input} GENOMEREF
		"""
}

CUTOFFS = ['WXS': 10, 'WGS': 10]

process ProcessVariants {
	tag "Process variants ${cohort}"
	label "core"
	errorStrategy 'ignore'  // if a cohort does not pass the filters, do not proceed with it
	publishDir "${STEPS_FOLDER}/variants", mode: "copy"

	input:
		tuple val(cohort), path(input), val(platform), val(genome) from COHORTS4.join(PLATFORMS1).join(GENOMES)

	output:
		tuple val(cohort), path(output) into VARIANTS
		tuple val(cohort), path("${output}.stats.json") into STATS_VARIANTS

	script:
		cutoff = CUTOFFS[platform]
		output = "${cohort}.tsv.gz"
		if (cutoff)
			"""
			parse-variants --input ${input} --output ${output} \
				--genome ${genome.toLowerCase()} \
				--cutoff ${cutoff}
			"""
		else
			error "Invalid cutoff. Check platform: $platform"

}

VARIANTS.into { VARIANTS1; VARIANTS2; VARIANTS3; VARIANTS4; VARIANTS5 }

process FormatVEP {
	tag "Prepare for VEP ${cohort}"
	label "core"
	publishDir "${STEPS_FOLDER}/vep", mode: "copy"

	input:
		tuple val(cohort), path(input) from VARIANTS5

	output:
		tuple val(cohort), path(output) into VARIANTS_VEP

	script:
		output = "${cohort}.in.tsv.gz"
		"""
		format-variants --input ${input} --output ${output}.tmp \
			--format vep
		sort -k1V -k2n -k3n ${output}.tmp |gzip > ${output}
		"""

}

process VEP {
	tag "VEP ${cohort}"
	publishDir "${STEPS_FOLDER}/vep", mode: "copy"

    input:
        tuple val(cohort), path(input) from VARIANTS_VEP

    output:
        tuple val(cohort), path(output) into OUT_VEP

	script:
		output = "${cohort}.vep.tsv.gz"
		"""
		vep -i ${input} -o STDOUT --assembly GRCh38 \
			--no_stats --cache --offline --symbol \
			--protein --tab --canonical --mane \
			--dir ${params.datasets}/vep \
			| grep -v ^## | gzip > ${output}
		"""
}

process ProcessVEPoutput {
	tag "Process vep output ${cohort}"
	label "core"
	publishDir "${STEPS_FOLDER}/vep", mode: "copy"

    input:
        tuple val(cohort), path(input) from OUT_VEP

    output:
        tuple val(cohort), path(output) into PARSED_VEP
        tuple val(cohort), path("${output}.stats.json") into STATS_VEP

	script:
		output = "${cohort}.tsv.gz"
		"""
		parse-vep --input ${input} --output ${output}
		"""
}


PARSED_VEP.into { PARSED_VEP1; PARSED_VEP2; PARSED_VEP3; PARSED_VEP4; PARSED_VEP5; PARSED_VEP6; PARSED_VEP7 }