#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""A module to draw boxplot of genotype.
"""

import argparse

from qmutils import BoxPlot


def get_args():
    """Get CLI arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--phtp-file", dest="phtp_file", required=True, help="The file read the phenotype level from.")
    parser.add_argument("-g", "--gntp-file", dest="gntp_file", required=True, help="The file read the genotype from.")
    parser.add_argument("-i", "--gntp-info-file", dest="gntp_info_file", required=True, help="The file read the genotype information from.")
    parser.add_argument("-c", "--cvrt-file", dest="cvrt_file", help="The file read the covriates from.")
    parser.add_argument("-e", "--effect-allele", dest="effal", default="EffectAllele", help="The allele that encoded as 1. Default: %(default)s")
    parser.add_argument("-P", "--phtp-name", dest="phtp_name", nargs="*", help="The name of phenotype for which the script draws the boxplot.")
    parser.add_argument("-G", "--snps-indx", dest="snps_indx", nargs="*", help="The genotype id for which the script draw boxplot.")
    parser.add_argument("-C", "--cvrt-name", dest="cvrt_name", nargs="*", help="The name of covriates by which the script correct.")
    parser.add_argument("-Y", "--y-info", dest="y_info", nargs="*", help="Information string for y-axis.")
    parser.add_argument("-o", "--output-prefix", dest="output_dir", default="./", help="The prefix for the output aka boxplot")

    return parser


def main():
    """The main entry.
    """
    args = get_args().parse_args()
    phtp_file = args.phtp_file
    phtp_name = args.phtp_name
    y_info = args.y_info
    gntp_file = args.gntp_file
    snps_indx = args.snps_indx
    gntp_info_file = args.gntp_info_file
    cvrt_file = args.cvrt_file
    cvrt_name = args.cvrt_name

    effal = args.effal

    output_dir = args.output_dir

    box_plot = BoxPlot(phtp_file, phtp_name, gntp_file, snps_indx,
                       gntp_info_file, cvrt_file, cvrt_name,
                       effect_allele=effal, y_info=y_info,
                       output_dir=output_dir)

    box_plot.init().draw_boxplot()


if __name__ == "__main__":
    main()
