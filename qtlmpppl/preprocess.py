#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""A script to pre-process traits for QTL mapping analysis.
"""

import argparse

from qmutils import PreProcess


def get_args():
    """Get CLI options.
    """
    parser = argparse.ArgumentParser(description="A script to preprocessing the phenotype data.")
    parser.add_argument("-c", "--conf-file", dest="conf_file", metavar="CONFIG", help="The file from which read configurations. Default: %(default)s")
    parser.add_argument("-p", "--phtp-file", dest="phtp_file", metavar="PHENOTYPE-FILE", help="The file from which read the phenotypes.")
    parser.add_argument("-t", "--phtp-list", dest="phtp_list", nargs="*", help="Traits to be processed. If is None, all numerical variable will be processed")
    parser.add_argument("-r", "--cvrt-file", dest="cvrt_file", metavar="COVARIATE-FILE", help="The file from which read the covariates.")
    parser.add_argument("-v", "--cvrt-list", dest="cvrt_list", nargs="*", metavar="COVARIATE", help="Traits used as covariates.")
    parser.add_argument("-o", "--output-pre", dest="output_pre", default="phenotype", metavar="OUTPUT_PREFIX", help="Prefix of output file after processed. Default: %(default)s")
    parser.add_argument("-s", "--stdv-times", dest="stdv_times", default=3, type=int, metavar="TIMES", help="Times of standard deviation. Default: %(default)s")

    return parser


def main():
    """The main entry of the module.
    """
    args = get_args().parse_args()

    phtp_file = args.phtp_file
    cvrt_file = args.cvrt_file
    conf_file = args.conf_file

    pre_processer = PreProcess(conf_file, phtp_file, cvrt_file)
    pre_processer.init() \
            .encode_sex() \
            .transform() \
            .check_outliers(skip=True) \
            .mask_outliers() \
            .check_dist() \
            .correlate() \
            .save_results()


if __name__ == "__main__":
    main()
