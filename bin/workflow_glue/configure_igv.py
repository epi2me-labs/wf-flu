"""Create an IGV config file."""

from itertools import zip_longest
import json
import sys

from .util import get_named_logger, wf_parser  # noqa: ABS101


def parse_fnames(fofn):
    """Parse list with filenames and return them grouped as ref-, BAM-, or VCF-related.

    :param fofn: File with list of file names (one per line)
    :return: dict of reference-related filenames (with keys 'ref', 'fai', and '.gzi' and
        `None` as default values); lists of BAM- and VCF-related filenames
    """
    ref_extensions = [".fasta", ".fasta.gz", ".fa", ".fa.gz", ".fna", ".fna.gz"]
    ref_dict = {}
    bams = []
    bam_indices = []
    vcfs = []
    vcf_indices = []
    with open(fofn, "r") as f:
        for line in f:
            fname = line.strip()
            if any(fname.endswith(ext) for ext in ref_extensions):
                ref_dict["ref"] = fname
            elif fname.endswith(".fai"):
                ref_dict["fai"] = fname
            elif fname.endswith(".gzi"):
                ref_dict["gzi"] = fname
            elif fname.endswith(".bam"):
                bams.append(fname)
            elif fname.endswith(".bai"):
                bam_indices.append(fname)
            elif fname.endswith(".vcf") or fname.endswith(".vcf.gz"):
                vcfs.append(fname)
            elif fname.endswith(".csi") or fname.endswith(".tbi"):
                vcf_indices.append(fname)
    # do some sanity checks
    if "ref" not in ref_dict:
        raise ValueError(
            "No reference file (i.e. file ending in one of "
            f"{ref_extensions} was found)."
        )
    ref = ref_dict["ref"]
    if (gzi := ref_dict.get("gzi")) is not None:
        # since we got a '.gzi' index, make sure that the reference is actually
        # compressed
        if not ref_dict["ref"].endswith(".gz"):
            raise ValueError(
                f"Found GZI reference index '{gzi}', but the reference file "
                f"'{ref}' appears not to be compressed."
            )
    if bam_indices:
        if len(bams) != len(bam_indices):
            raise ValueError("Got different number of BAM and BAM index files.")
    if vcf_indices:
        if len(vcfs) != len(vcf_indices):
            raise ValueError("Got different number of VCF and VCF index files.")
    if bams and vcfs:
        if len(bams) != len(vcfs):
            raise ValueError("Got different number of BAM and VCF files.")
    # if we got BAM or VCF indices, pair them up with their corresponding files (and
    # otherwise with `None`)
    bams_with_indices = zip_longest(bams, bam_indices)
    vcfs_with_indices = zip_longest(vcfs, vcf_indices)
    return ref_dict, bams_with_indices, vcfs_with_indices


def get_reference_options(ref, fai=None, gzi=None):
    """Create dict with IGV reference options.

    :param ref: reference file name
    :param fai: name reference `.fai` index file
    :param gzi: name of `.gzi` index file for a compressed reference
    :return: dict with reference options
    """
    # initialise the options dict and add the index attributes later
    ref_opts = {
        "id": "ref",
        "name": "ref",
        "wholeGenomeView": False,
        "fastaURL": ref,
    }
    if fai is not None:
        ref_opts["indexURL"] = fai
    if gzi is not None:
        ref_opts["compressedIndexURL"] = gzi
    return ref_opts


def get_alignment_track(bam, bai=None, extra_opts=None):
    """Create dict with options for IGV alignment track.

    :param bam: name of BAM file to be displayed
    :param bai: name of BAM index file
    :param extra_opts: dict of extra options for the alignment track
    :return: dict with alignment track options
    """
    alignment_track_dict = {
        "name": bam,
        "type": "alignment",
        "format": "bam",
        "url": bam,
    }
    # add the BAM index if present
    if bai is not None:
        alignment_track_dict["indexURL"] = bai
    alignment_track_dict.update(extra_opts or {})
    return alignment_track_dict


def get_variant_track(vcf, index=None, extra_opts=None):
    """Create dict with options for IGV variant track.

    :param vcf: name of VCF file to be displayed
    :param index: name of VCF index file (ending in `.csi` or `.tbi`)
    :param extra_opts: dict of extra options for the variant track
    :return: dict with variant track options
    """
    variant_track_dict = {
        "name": vcf,
        "type": "variant",
        "format": "vcf",
        "url": vcf,
    }
    # add the VCF index if we got an index extension
    if index is not None:
        variant_track_dict["indexURL"] = index
    variant_track_dict.update(extra_opts or {})
    return variant_track_dict


def main(args):
    """Run the entry point."""
    logger = get_named_logger("configIGV")

    # parse the FOFN
    ref_dict, bams_with_indices, vcfs_with_indices = parse_fnames(args.fofn)

    # initialise the IGV options dict with the reference options
    json_dict = {"reference": get_reference_options(**ref_dict)}

    # if we got JSON files with extra options for the alignment / variant tracks, read
    # them
    extra_alignment_opts = {}
    if args.extra_bam_opts is not None:
        with open(args.extra_bam_opts, "r") as f:
            extra_alignment_opts = json.load(f)
    extra_variant_opts = {}
    if args.extra_vcf_opts is not None:
        with open(args.extra_vcf_opts, "r") as f:
            extra_variant_opts = json.load(f)

    # now add the alignment and variant tracks
    json_dict["tracks"] = []
    # we use `zip_longest` to make sure that variant and alignment tracks from the same
    # sample are added after each other
    for (vcf, vcf_index), (bam, bam_index) in zip_longest(
        vcfs_with_indices, bams_with_indices, fillvalue=(None, None)
    ):
        if vcf is not None:
            # add an variant track for the VCF
            json_dict["tracks"].append(
                get_variant_track(vcf, vcf_index, extra_variant_opts)
            )
        if bam is not None:
            # add an alignment track for the BAM
            json_dict["tracks"].append(
                get_alignment_track(bam, bam_index, extra_alignment_opts)
            )

    if args.locus is not None:
        json_dict["locus"] = args.locus

    json.dump(json_dict, sys.stdout, indent=4)

    logger.info("Printed IGV config JSON to STDOUT.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("configure_igv")
    parser.add_argument(
        "--fofn",
        required=True,
        help=(
            "File with list of names of reference / BAM / VCF files and indices "
            "(one filename per line)"
        ),
    )
    parser.add_argument(
        "--locus",
        help="Locus string to set initial genomic coordinates to display in IGV",
    )
    parser.add_argument(
        "--extra-bam-opts",
        help="JSON file with extra options for alignment tracks",
    )
    parser.add_argument(
        "--extra-vcf-opts",
        help="JSON file with extra options for variant tracks",
    )
    return parser
