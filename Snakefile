# Snakefile
import os
from urllib.parse import urlparse

configfile: "config.yaml"

# --- Helper function to get filename from URL ---
def get_filename_from_url(url):
    path = urlparse(url).path
    return os.path.basename(path)

# --- Dynamically generate all target files for the 'all' rule ---
ALL_TARGETS = []

# 1. Reference Genomes
# Uses a set to ensure each unique reference is listed only once.
UNIQUE_REFERENCE_TARGETS = set()
for sample_name, sample_info in config["samples"].items():
    ref_key = sample_info.get("reference_key")
    if not ref_key:
        # Fallback to reference_build if reference_key is not specified
        ref_key = sample_info.get("reference_build", "GRCh38") # Default to GRCh38 or make it mandatory
        print(f"Warning: 'reference_key' not specified for sample {sample_name}. Defaulting to '{ref_key}'.")

    if ref_key in config["references"]:
        ref_config = config["references"][ref_key]
        ref_fasta_path = os.path.join(
            config["output_ref_dir"].format(data_dir=config["data_dir"]),
            ref_config["filename"]
        )
        UNIQUE_REFERENCE_TARGETS.add(ref_fasta_path)
    else:
        raise ValueError(
            f"Reference key '{ref_key}' for sample '{sample_name}' not found in 'references' section of config.yaml."
        )
ALL_TARGETS.extend(list(UNIQUE_REFERENCE_TARGETS))


# 2. GIAB Truth Sets (VCF, TBI, BED) per sample
for sample_name, sample_info in config["samples"].items():
    giab_version = sample_info["giab_version"]
    ref_build = sample_info["reference_build"]

    output_giab_dir = config["output_giab_dir_template"].format(
        data_dir=config["data_dir"],
        sample_name=sample_name,
        giab_version=giab_version,
        reference_build=ref_build
    )

    # Construct the basename for truth files for this sample
    truth_basename = config["giab_truth_basename_template"].format(
        sample_name=sample_name,
        reference_build=ref_build,
        giab_version=giab_version
    )

    # Construct full paths for VCF, TBI, BED
    giab_vcf = os.path.join(output_giab_dir, config["giab_vcf_filename_template"].format(giab_truth_basename=truth_basename))
    giab_tbi = os.path.join(output_giab_dir, config["giab_tbi_filename_template"].format(giab_truth_basename=truth_basename))
    giab_bed = os.path.join(output_giab_dir, config["giab_bed_filename_template"].format(giab_truth_basename=truth_basename))

    ALL_TARGETS.extend([giab_vcf, giab_tbi, giab_bed])

# 3. FASTQ files per sample
for sample_name, sample_info in config["samples"].items():
    output_fastq_dir = config["output_fastq_dir_template"].format(
        data_dir=config["data_dir"],
        sample_name=sample_name
    )
    if "fastq_urls" in sample_info and sample_info["fastq_urls"]:
        for fq_url in sample_info["fastq_urls"]:
            if fq_url: # Ensure URL is not empty or None
                fastq_filename = get_filename_from_url(fq_url)
                ALL_TARGETS.append(os.path.join(output_fastq_dir, fastq_filename))
    else:
        print(f"Warning: No FASTQ URLs provided for sample {sample_name}.")


# --- Main Rule ---
rule all:
    input:
        ALL_TARGETS
    shell:
        "echo 'All GIAB data and FASTQs specified in config.yaml downloaded successfully.'"

# --- Download Rules ---

# Rule to download reference genomes
rule download_reference:
    output:
        # Path: {data_dir}/reference/{ref_filename}
        # Wildcards: ref_filename (data_dir is fixed from config)
        os.path.join(config["output_ref_dir"].format(data_dir=config["data_dir"]), "{ref_filename}")
    params:
        # Find the URL for the specific ref_filename from the config
        get_ref_url = lambda wildcards: next(
            (ref_info["url"] for ref_key, ref_info in config["references"].items()
            if ref_info["filename"] == wildcards.ref_filename), None
        )
    shell:
        """
        set -eo pipefail
        if [ -z "{params.get_ref_url}" ]; then
            echo "Error: Could not find URL for reference file {wildcards.ref_filename}" >&2
            exit 1
        fi
        mkdir -p $(dirname {output})
        echo "Downloading reference: {params.get_ref_url}"
        wget -nv -c "{params.get_ref_url}" -O {output}
        """

# Rule to download reference genome index (.fai)
rule download_reference_fai:
    output:
        os.path.join(config["output_ref_dir"].format(data_dir=config["data_dir"]), "{ref_filename}.fai")
    params:
        get_ref_url = lambda wildcards: next(
            (ref_info["url"] for ref_key, ref_info in config["references"].items()
            if ref_info["filename"] == wildcards.ref_filename), None
        )
    shell:
        """
        set -eo pipefail
        if [ -z "{params.get_ref_url}" ]; then
            echo "Error: Could not find URL for reference file {wildcards.ref_filename}" >&2
            exit 1
        fi
        mkdir -p $(dirname {output})
        echo "Downloading reference index: {params.get_ref_url}.fai"
        wget -nv -c "{params.get_ref_url}.fai" -O {output}
        """

# Rule to download reference genome gzip index (.gzi)
rule download_reference_gzi:
    output:
        os.path.join(config["output_ref_dir"].format(data_dir=config["data_dir"]), "{ref_filename}.gzi")
    params:
        get_ref_url = lambda wildcards: next(
            (ref_info["url"] for ref_key, ref_info in config["references"].items()
            if ref_info["filename"] == wildcards.ref_filename), None
        )
    shell:
        """
        set -eo pipefail
        if [ -z "{params.get_ref_url}" ]; then
            echo "Error: Could not find URL for reference file {wildcards.ref_filename}" >&2
            exit 1
        fi
        mkdir -p $(dirname {output})
        echo "Downloading reference gzip index: {params.get_ref_url}.gzi"
        wget -nv -c "{params.get_ref_url}.gzi" -O {output}
        """

# Rule to download GIAB truth VCF files
rule download_giab_vcf:
    output:
        # Path: {data_dir}/giab/{sample_name}/{giab_version}/{reference_build}/{sample_name}_{reference_build}_{giab_version}_benchmark.vcf.gz
        # Wildcards: sample_name, giab_version, reference_build (data_dir is fixed)
        # The output pattern is constructed by joining the dir template and filename template parts, with wildcards.
        os.path.join(
            config["output_giab_dir_template"], # contains {sample_name}, {giab_version}, {reference_build}
            config["giab_vcf_filename_template"].format(
                giab_truth_basename=config["giab_truth_basename_template"] # also contains {sample_name}, {giab_version}, {reference_build}
            )
        ).replace("{data_dir}", config["data_dir"]) # Fix data_dir part
    params:
        get_url = lambda wildcards: os.path.join(
            config["giab_base_url"],
            f"{config['samples'][wildcards.sample_name]['id']}_{wildcards.sample_name}",
            wildcards.giab_version,
            wildcards.reference_build,
            config["giab_vcf_filename_template"].format( # Reconstruct actual filename for URL
                giab_truth_basename=config["giab_truth_basename_template"].format(
                    sample_name=wildcards.sample_name,
                    reference_build=wildcards.reference_build,
                    giab_version=wildcards.giab_version
                )
            )
        )
    shell:
        """
        set -eo pipefail
        mkdir -p $(dirname {output})
        echo "Downloading VCF for {wildcards.sample_name}: {params.get_url}"
        wget -nv -c "{params.get_url}" -O {output}
        """

# Rule to download GIAB truth VCF index (.tbi) files
rule download_giab_tbi:
    output:
        os.path.join(
            config["output_giab_dir_template"],
            config["giab_tbi_filename_template"].format(
                giab_truth_basename=config["giab_truth_basename_template"]
            )
        ).replace("{data_dir}", config["data_dir"])
    params:
        get_url = lambda wildcards: os.path.join(
            config["giab_base_url"],
            f"{config['samples'][wildcards.sample_name]['id']}_{wildcards.sample_name}",
            wildcards.giab_version,
            wildcards.reference_build,
            config["giab_tbi_filename_template"].format(
                giab_truth_basename=config["giab_truth_basename_template"].format(
                    sample_name=wildcards.sample_name,
                    reference_build=wildcards.reference_build,
                    giab_version=wildcards.giab_version
                )
            )
        )
    shell:
        """
        set -eo pipefail
        mkdir -p $(dirname {output})
        echo "Downloading TBI for {wildcards.sample_name}: {params.get_url}"
        echo "DEBUG: URL is {params.get_url}"
        wget -nv -c "{params.get_url}" -O {output}
        """

# Rule to download GIAB high-confidence BED files
rule download_giab_bed:
    output:
        os.path.join(
            config["output_giab_dir_template"],
            config["giab_bed_filename_template"].format(
                giab_truth_basename=config["giab_truth_basename_template"]
            )
        ).replace("{data_dir}", config["data_dir"])
    params:
        get_url=lambda wildcards: config["samples"][wildcards.sample_name].get("bed_url")
    shell:
        """
        set -eo pipefail
        mkdir -p $(dirname {output})
        echo "Downloading BED for {wildcards.sample_name}: {params.get_url}"
        wget -nv -c "{params.get_url}" -O {output}
        """

# Rule to download FASTQ files
rule download_fastq:
    output:
        # Path: {data_dir}/fastq/{sample_name}/{fastq_filename}
        # Wildcards: sample_name, fastq_filename (data_dir is fixed)
        os.path.join(
            config["output_fastq_dir_template"],
            "{fastq_filename}"
        ).replace("{data_dir}", config["data_dir"]) # Fix data_dir part for the pattern
    params:
        # Find the original URL for this specific sample_name and fastq_filename
        get_fastq_url = lambda wildcards: next(
            (url for url in config["samples"][wildcards.sample_name].get("fastq_urls", [])
            if url and get_filename_from_url(url) == wildcards.fastq_filename),
            None # Return None if not found to prevent error, though rule shouldn't trigger
        )
    wildcard_constraints:
        sample_name="|".join(config["samples"].keys()),
    shell:
        """
        set -eo pipefail
        if [ -z "{params.get_fastq_url}" ]; then
            echo "Error: Could not find URL for FASTQ {wildcards.fastq_filename} of sample {wildcards.sample_name}" >&2
            exit 1
        fi
        mkdir -p $(dirname {output})
        echo "Downloading FASTQ for {wildcards.sample_name}: {params.get_fastq_url}"
        echo "DEBUG: URL is {params.get_fastq_url}"
        wget -nv -c "{params.get_fastq_url}" -O {output}
        """
