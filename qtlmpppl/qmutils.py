#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""QTL mapping pipeline utils.
"""

import re
import json
import logging

import numpy as np
import pandas as pd
import seaborn as sb
import pingouin as pg
import scipy.stats as stt
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf

try:
    import rpy2.robjects as robj
    import rpy2.robjects.packages as rpk
except:
    pass

from sklearn.decomposition import PCA

# Set up a global logger.
logger = logging.getLogger("matplotlib")
logger.setLevel(logging.ERROR)

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

cs_stream = logging.StreamHandler()
fmt = logging.Formatter("| {levelname: ^8} | {asctime} | {name}: {message}", datefmt="%Y-%m-%d %H:%M:%S %p", style="{")
cs_stream.setFormatter(fmt)
cs_stream.setLevel(logging.INFO)

logger.addHandler(cs_stream)


CHROM_LEN_GRCH37 = {
    'chr1' : 249250621, 'chr2' : 243199373, 'chr3' : 198022430,
    'chr4' : 191154276, 'chr5' : 180915260, 'chr6' : 171115067,
    'chr7' : 159138663, 'chr8' : 146364022, "chr9" : 141213431,
    "chr10": 135534747, "chr11": 135006516, "chr12": 133851895,
    "chr13": 115169878, "chr14": 107349540, "chr15": 102531392,
    "chr16":  90354753, "chr17":  81195210, "chr18":  78077248,
    "chr19":  59128983, "chr20":  63025520, "chr21":  48129895,
    "chr22":  51304566, "chrx" : 155270560, "chry" :  59373566,
}


CELL_SUBSET_DICT = {
    "CD45": "Leukocytes",
    "M": "Monocytes",
    "L": "Lymphocytes",
    "CD4": "CD4+",
    "RApR7p": "CD4+Naive",
    "RAnR7p": "CD4+CM",
    "mTreg": "CD4+mTreg",
    "nTreg": "CD4+nTreg",
    "TEM": "CD4+TEM",
    "RAnR7n": "CD4+EM",
    "RApR7n": "CD4+TEMRA",
    "CD8": "CD8+",
    "Naive_CD8": "CD8+Naive",
    "CM_CD8": "CD8+CM",
    "EM_CD8": "CD8+TEM",
    "RAnR7n_EM_CD8": "CD8+EM",
    "RApR7n_EM_CD8": "CD8+TEMRA",
}



class DataSet:
    """A class to handle data set.
    """
    def __init__(self, file_path=None):
        self.file_path = file_path

    def load_data(self, file_path, fmt='tsv'):
        """Load data from given file path"""
        return self

    def dump_data(self, file_path, fmt='tsv'):
        """Dump processed data into given file path"""
        return self


class BoxPlot:
    """Draw boxplot.
    """
    def __init__(self, phtp_file, phtp_name, gntp_file, snps_indx, gntp_info_file,
                 cvrt_file=None, cvrt_name=None, effect_allele=None,
                 y_info=None, output_dir="./"):
        self.phtp_file = phtp_file
        self.phtp_name = phtp_name

        if y_info is None:
            self.y_info = self.phtp_name
        else:
            self.y_info = y_info

        self.gntp_file = gntp_file
        self.snps_indx = snps_indx
        self.gntp_info_file = gntp_info_file
        self.cvrt_file = cvrt_file
        self.cvrt_name = cvrt_name
        self.effect_allele = effect_allele

        self.output_dir = output_dir

        self.phtp_dtfm = None
        self.cvrt_dtfm = None
        self.gntp_dtfm = None
        self.gntp_info_dtfm = None
        self.dtfm = None

    @staticmethod
    def _pr_clear():
        plt.cla()
        plt.clf()
        plt.close()

    @staticmethod
    def _pr_load_file(file_path, **kwargs):
        sep = kwargs.pop("sep") if "sep" in kwargs else "\t"
        dtfm = pd.read_csv(file_path, sep=sep, **kwargs)
        return dtfm

    def _pr_make_dtfm(self):
        # Load phenotype data
        phtp_dtfm = self._pr_load_file(self.phtp_file, index_col=0)

        if self.phtp_name is None:
            self.phtp_dtfm = phtp_dtfm
            self.phtp_name = phtp_dtfm.index
        else:
            self.phtp_name = [phtp for phtp in self.phtp_name if phtp in phtp_dtfm.index]
            self.phtp_dtfm = phtp_dtfm.loc[self.phtp_name]

        # Load genotype data
        gntp_dtfm = self._pr_load_file(self.gntp_file, index_col=0)
        gntp_info_dtfm = self._pr_load_file(self.gntp_info_file, index_col=0)

        if self.snps_indx is None:
            self.gntp_dtfm = gntp_dtfm
            self.snps_indx = gntp_dtfm.index
            self.gntp_info_dtfm = gntp_info_dtfm
        else:
            self.snps_indx = [
                snps_indx for snps_indx in self.snps_indx
                if snps_indx in gntp_dtfm.index and snps_indx in gntp_info_dtfm.index]

            self.gntp_dtfm = gntp_dtfm.loc[self.snps_indx]
            self.gntp_info_dtfm = gntp_info_dtfm.loc[self.snps_indx]


        if len(self.snps_indx) > 16:
            raise ValueError("More than 16 snps were selected. Unsupported!!!")

        # Load covariates
        if self.cvrt_file is not None:
            cvrt_dtfm = self._pr_load_file(self.cvrt_file, index_col=0)

            if self.cvrt_name is None:
                self.cvrt_dtfm = cvrt_dtfm
                self.cvrt_name = cvrt_dtfm.index
            else:
                self.cvrt_name = [cvrt for cvrt in self.cvrt_name if cvrt in cvrt_dtfm.index]
                self.cvrt_dtfm = cvrt_dtfm.loc[self.cvrt_name]

        # Merge the previous three dataframe
        self.dtfm = pd.concat([self.phtp_dtfm, self.gntp_dtfm, self.cvrt_dtfm], join="inner")

    def _pr_encode_gntp(self):
        for gntp in self.snps_indx:
            eff, alt = self.gntp_info_dtfm.loc[gntp, ["EffectAllele", "AlternativeAllele"]]

            dosage2code_dict = {0: eff + eff, 1: eff + alt, 2: alt + alt}
            if self.effect_allele == "EffectAllele":
                dosage2code_dict = {0: alt + alt, 1: alt + eff, 2: eff + eff}

            self.dtfm.loc[gntp + "_code"] = self.dtfm \
                    .loc[gntp, :] \
                    .apply(lambda x: dosage2code_dict[round(x)])

    def _pr_correct_cvrt(self):
        for phtp in self.phtp_name:
            dtfm = self.dtfm.loc[[phtp]]
            smf.ols('', data=dtfm).fit()

    def _pr_draw_boxplots(self, svfmt="pdf", **kwargs):
        sb.set(style="ticks")
        width = kwargs.pop("width") if "width" in kwargs else 6
        height = kwargs.pop("height") if "height" in kwargs else 8

        for phtp, y_info in zip(self.phtp_name, self.y_info):
            for gntp in self.snps_indx:
                gntp_code = gntp + "_code"
                dtfm = self.dtfm.loc[[phtp, gntp_code]]
                chrom, pos, eff, alt, *_ = self.gntp_info_dtfm.loc[gntp, ]

                if self.effect_allele == "EffectAllele":
                    allele_pair = [alt, eff]
                    allele_order = [alt + alt, alt + eff, eff + eff]
                else:
                    allele_pair = [eff, alt]
                    allele_order = [eff + eff, eff + alt, alt + alt]

                allele_count = (dtfm.loc[gntp_code]
                                .value_counts()
                                .reindex(allele_order)
                                .fillna(0))

                xtick_labels = ['{} ({})'.format(aa, cc)
                                for aa, cc in zip(allele_order, allele_count)]

                order_dict = dict(zip(allele_order, range(len(allele_count))))
                dtfm = dtfm.transpose().sort_values(by=gntp_code, key=lambda x: x.map(order_dict))

                axes = sb.boxplot(x=gntp_code, y=phtp, data=dtfm, width=0.4, color='0.9')#, order=allele_code)
                axes = sb.swarmplot(x=gntp_code, y=phtp, data=dtfm, color="0", alpha=0.8)#, order=allele_code)

                axes.spines["top"].set_visible(False)
                axes.spines["right"].set_visible(False)

                axes.set_title("Boxplot for {}".format(gntp))
                axes.set_xlabel("Genptype (chr:{}, pos:{}, allele:{}>{})".format(chrom, pos, *allele_pair))
                axes.set_ylabel(y_info)
                axes.set_xticklabels(xtick_labels)

                fig_name = ".".join(["boxplot", phtp, gntp, svfmt])
                fig_name = self.output_dir + "/" + fig_name.replace('+', 'p')

                fig = axes.get_figure()
                fig.set_figheight(height)
                fig.set_figwidth(width)
                fig.savefig(fig_name)

                self._pr_clear()

    def init(self,):
        """Init."""
        self._pr_make_dtfm()
        self._pr_encode_gntp()
        return self

    def stats(self):
        """Do statistics."""
        self._pr_correct_cvrt()
        return self

    def draw_boxplot(self):
        """Draw boxplot."""
        self._pr_draw_boxplots()
        return self


class PreProcess:
    """A class to pre-process phenotypes.
    """
    def __init__(self, conf_file=None, phtp_file=None, cvrt_file=None):
        self.conf_file = conf_file
        self.phtp_file = phtp_file
        self.cvrt_file = cvrt_file

        self.configs = None   # Configurations parsed from self.conf_file
        self.dataframe = None # The data frame loaded from self.phtp_file

        self.phtp_list = None
        self.phtp_dtfm = None
        self.phtp_trans_dict = None
        self.phtp_transed_name_list = None

        self.smpl_idx = None
        self.smpl_id_col = None
        self.smpl_list = None
        self.smpl_blck_list = None

        self.cvrt_list = None
        self.cvrt_dtfm = None

        self.sex_col = None

        self.pcor_dtfm = None

        self.transpose = False
        self.stdv_times = None
        self.output_prefix = None

    def _pr_load_config(self, **kwargs):
        with open(self.conf_file, "r") as config_file_handle:
            self.configs = json.load(config_file_handle, **kwargs)

    def _pr_parse_phtp_list(self):
        self.phtp_list = self.configs["preprocess"]["phtp_list"]

    def _pr_parse_cvrt_list(self):
        self.cvrt_list = self.configs["preprocess"]["cvrt_list"]

    def _pr_parse_pntp_tran_dict(self):
        phtp_trans_dict = self.configs["preprocess"]["trans_dict"]
        _kept_cols = phtp_trans_dict.pop("orig")

        col_to_func_dict = {}
        transformed_cols_vec = []
        for key, val_vec in phtp_trans_dict.items():
            if key in ["log2", "log10"]:
                col_to_func_dict.update({x: key for x in val_vec})
            elif key in ["ivrk"]:
                col_to_func_dict.update({x: self._inverse_rank for x in val_vec})
            else:
                logger.error("Only log2, log10, and inverse rank are supported")
                continue

            transformed_cols_vec.extend([x + "_" + key for x in val_vec])
        self.phtp_trans_dict = col_to_func_dict

    def _pr_parse_misc(self):
        if self.phtp_file is None:
            self.phtp_file = self.configs["preprocess"]["phtp_file"]

        if self.cvrt_file is None:
            self.cvrt_file = self.configs["preprocess"]["cvrt_file"]

        self.transpose = self.configs["preprocess"]["transpose"]
        self.output_prefix = self.configs["preprocess"]["output_prefix"]
        self.stdv_times = self.configs["preprocess"]["stdv_times"]
        self.smpl_id_col = self.configs["preprocess"]["smpl_id_col"]
        self.sex_col = self.configs["preprocess"]["sex_col"]

    def _pr_parse_sample_id(self):
        self.smpl_idx = self.configs["preprocess"]["smpl_idx"]

        if self.smpl_idx is not None:
            smpl_idx = []
            for idx_range in self.smpl_idx.split(","):
                if idx_range:
                    if "-" in idx_range:
                        start, stop = idx_range.split("-")
                        smpl_idx.extend(list(range(int(start), int(stop))))
                    else:
                        smpl_idx.append(int(idx_range))
                else:
                    logger.debug("Empty range of sample index")
            self.smpl_idx = smpl_idx

        self.smpl_list = self.configs["preprocess"]["smpl_list"]
        self.smpl_blck_list = self.configs["preprocess"]["smpl_blck_list"]

    def _pr_load_dtfm(self, **kwargs):
        sep = kwargs.pop("sep") if "sep" in kwargs else "\t"

        phtp_dtfm = pd.read_csv(self.phtp_file, sep=sep, **kwargs)
        if self.cvrt_file is None or self.cvrt_file == "":
            self.dataframe = phtp_dtfm
        else:
            cvrt_dtfm = pd.read_csv(self.cvrt_file, sep=sep, **kwargs)
            self.dataframe = pd.merge(phtp_dtfm, cvrt_dtfm, how="outer", on=self.smpl_id_col)

        if self.smpl_idx is not None and self.smpl_idx != "":
            self.dataframe = self.dataframe.loc[self.smpl_idx, :]
            self.dataframe.index = range(len(self.smpl_idx))

        smpl_list = self.smpl_list
        smpl_id_col = self.smpl_id_col
        smpl_blck_list = self.smpl_blck_list
        if self.smpl_list is not None:
            self.dataframe = self.dataframe.loc[self.dataframe[smpl_id_col].isin(smpl_list), :]

        if self.smpl_blck_list is not None:
            self.dataframe = self.dataframe.loc[self.dataframe[smpl_id_col].isin(smpl_blck_list) == False, :]

    @staticmethod
    def _inverse_rank(raw_list):
        # Inverse rank transformer.
        return stt.norm.ppf(stt.rankdata(raw_list) / (len(raw_list) + 1))

    def init(self, **kwargs):
        """Initialize the processor.
        """
        if self.conf_file:
            self._pr_load_config()
        else:
            logger.info("No configuration file.")

        self._pr_parse_misc()
        self._pr_parse_phtp_list()
        self._pr_parse_cvrt_list()
        self._pr_parse_sample_id()
        self._pr_parse_pntp_tran_dict()

        sep = kwargs.get("sep") if "sep" in kwargs else ","
        self._pr_load_dtfm(sep=sep)
        return self

    def encode_sex(self, sex_col=None, mapping=None, skip=False):
        """Encode gender from string into binary"""
        if skip:
            return self

        if sex_col is None:
            if self.sex_col is None:
                return self
            sex_col = self.sex_col

        mapping = {"Male": 0, "Female": 1} if mapping is None else mapping
        self.dataframe[sex_col] = self.dataframe[sex_col].apply(lambda x: mapping[x])

        return self

    def transform(self, skip=False):
        """Transform the target columns using given math function.

        Note:
            1. The `pandas.DataFrame.transform()` method could be a good choice,
            as it allow to specify the function to apply on each column by a
            axis labels->functions implementation.
        """
        if skip:
            return self

        phtp_trans_dict = self.phtp_trans_dict
        phtp_tobe_transed = phtp_trans_dict.keys()
        if phtp_tobe_transed:
            phtp_transed_dtfm = (self.dataframe
                                 .loc[:, phtp_tobe_transed]
                                 .transform(phtp_trans_dict, axis=0))

            phtp_transed_name_list = []
            for col_name in phtp_transed_dtfm.columns:
                tran_func = phtp_trans_dict[col_name]
                func_name = tran_func if isinstance(tran_func, str) else "ivrk"
                phtp_transed_name_list.append(col_name + "_" + func_name)

            phtp_transed_dtfm.columns = phtp_transed_name_list
            self.dataframe[phtp_transed_name_list] = phtp_transed_dtfm
            self.phtp_transed_name_list = phtp_transed_name_list
        else:
            self.phtp_transed_name_list = []

        self.dataframe.replace((np.inf, -np.inf), np.nan, inplace=True)

        return self

    def check_outliers(self, figfmt="png", skip=False):
        """Perform a PCA analysis and show the results.

        Note:
            1. It works for dataset with at least 2 traits.
        """
        if skip:
            return self

        pca = PCA()
        transformed_values = pca.fit_transform(self.dataframe.loc[:, self.phtp_list])

        fig, axes = plt.subplots()
        axes.scatter(transformed_values[:, 0], transformed_values[:, 1], s=0.5)

        pca_save_name = ".".join([self.output_prefix, "pca", figfmt])
        fig.savefig(pca_save_name)

        return self

    def mask_outliers(self, skip=False):
        """Remove outliers.
        """
        if skip:
            return self

        stdv_times = self.stdv_times
        def _mask_outliers(vec: pd.Series, stdv_times):
            vec_mean = vec.mean()
            vec_stdv = vec.std()
            upper = vec_mean + vec_stdv * stdv_times
            lower = vec_mean - vec_stdv * stdv_times
            vec[((lower > vec) | (vec > upper))] = np.nan

            return vec

        phtp_list = self.phtp_list
        self.dataframe.loc[:, phtp_list] = self.dataframe.loc[:, phtp_list].transform(_mask_outliers, stdv_times=stdv_times)

        return self

    def check_dist(self, figfmt="png", skip=False):
        """Draw figures to show the distribution of the data.
        """
        if skip:
            return self

        if self.phtp_transed_name_list is None:
            pntp_check_dist_list = self.phtp_list
        else:
            pntp_check_dist_list = self.phtp_list + self.phtp_transed_name_list

        for col_name in pntp_check_dist_list:
            fig, axes = plt.subplots()
            self.dataframe.loc[:, col_name].plot(ax=axes, kind="hist")

            hist_opt_name = re.sub("[()/]", "", col_name.replace(" ", "_"))
            hist_opt_name = ".".join([self.output_prefix, hist_opt_name, figfmt])
            fig.savefig(hist_opt_name)

            plt.cla()
            plt.clf()
            plt.close()
        return self

    def correlate(self, x_list=None, y_list=None, c_list=None, method="spearman", figfmt="png", skip=False):
        """Correlate variables.
        """
        if skip:
            return self

        x_list = self.phtp_list if x_list is None else x_list
        y_list = self.phtp_list if y_list is None else y_list
        c_list = self.cvrt_list if c_list is None else c_list

        x_len, y_len = len(x_list), len(y_list)
        pcorr_mtrx_np = np.eye(x_len, y_len)
        for idx_x, pntp_x in enumerate(x_list):
            for idx_y, pntp_y in enumerate(y_list):
                if pntp_y != pntp_x:
                    try:
                        _pcorr_dtfm = pg.partial_corr(self.dataframe, pntp_x, pntp_y, self.cvrt_list, method=method)
                        pcorr_mtrx_np[idx_x, idx_y] = _pcorr_dtfm['r']
                    except AssertionError as err:
                        pass

        self.pcor_dtfm = pd.DataFrame(pcorr_mtrx_np, index=self.phtp_list, columns=self.phtp_list)
        pcorr_htmp_name = ".".join([self.output_prefix, "correlation_heatmap", figfmt])
        ctmp_grid = sb.clustermap(self.pcor_dtfm, col_cluster=True, row_cluster=True, cmap="Greens")
        ctmp_grid.fig.savefig(pcorr_htmp_name)

        plt.cla()
        plt.clf()
        plt.close()

        return self

    def save_results(self, dtfm_fmt="tsv", skip=False):
        """Save the processed results into disk.
        """
        if skip:
            return self

        if dtfm_fmt == "tsv":
            sep = "\t"
        elif dtfm_fmt == "csv":
            sep = ","
        elif dtfm_fmt == "ssv":
            sep = " "
        else:
            logger.info("Only ssv, tsv and csv are supported, use tsv by defult")
            dtfm_fmt, sep = "tsv", "\t"

        if self.smpl_id_col:
            self.dataframe.index = self.dataframe.loc[:, self.smpl_id_col]

        phtp_tosave_list = [x for x in self.phtp_list if x not in self.phtp_trans_dict] + self.phtp_transed_name_list
        self.phtp_dtfm = self.dataframe.loc[:, phtp_tosave_list]
        phtp_dtfm_name = ".".join([self.output_prefix, "proc_phtp", dtfm_fmt])

        if self.transpose:
            self.phtp_dtfm = self.phtp_dtfm.transpose()

        self.phtp_dtfm.to_csv(phtp_dtfm_name, sep=sep)

        if self.cvrt_list:
            cvrt_cols_to_save_list = self.cvrt_list
            self.cvrt_dtfm = self.dataframe.loc[:, cvrt_cols_to_save_list]

            if self.transpose:
                self.cvrt_dtfm = self.cvrt_dtfm.transpose()

            cvrt_dtfm_name = ".".join([self.output_prefix, "proc_cvrt", dtfm_fmt])
            self.cvrt_dtfm.to_csv(cvrt_dtfm_name, sep=sep)

        if self.pcor_dtfm is not None:
            pcorr_mtrx_name = ".".join([self.output_prefix, "pcorr_phtp", dtfm_fmt])
            self.pcor_dtfm.to_csv(pcorr_mtrx_name, sep=sep, index=True)

        return self


class QTLMapping:
    '''A wrapper of R package MatrixEQTL.
    '''
    def __init__(self, phtp=None, gntp=None, cvrt=None,
                 output_pref="./"):
        rpk.importr('MatrixEQTL')

        self.gntp = self._make_sliced_data(gntp)
        self.phtp = self._make_sliced_data(phtp)
        self.cvrt = self._make_sliced_data(cvrt)

        self.output_pref = output_pref

    def _load_dtfms(self):
        pass

    def _make_sliced_data(self, data):
        return ''

    def _matrix_eqtl_engine(self, optfn=None, pv_opt_thrd=None, use_model=None):
        if use_model is None:
            use_model = robj.r['modelLINEAR']

        if pv_opt_thrd is None:
            pv_opt_thrd = 0.05

        err_cov = robj.r['numeric']()
        matrix_eqtl_engine = robj.r['Matrix_eQTL_engine']
        matrix_eqtl_engine(
            snps=self.gntp, gene=self.phtp, cvrt=self.cvrt,
            output_file_name=optfn, pvOutputThreshold=pv_opt_thrd,
            useModel=use_model, errorCovariance=err_cov, verbose=False,
            pvalue_hist=False, min_pv_by_genesnp=True)

    def _estimate_dispersoin(self):
        pass

    def map_qtls(self, optfn=None, pv_opt_thrd=0.05, use_model=None):
        self._matrix_eqtl_engine(None, None)

        return self

    def draw_qq_plot(self):
        return self

    def draw_manhattan_plot(self):
        return self

    def save_results(self, output_path):
        return self


class LocusZoom:
    '''Draw locus zoom plot.'''
    def __init__(self):
        pass


class CollectAndReportSNP:
    """Collect and report SNPs.
    """
    def __init__(self, input_file_pl):
        self.input_file_pl = input_file_pl

    def _load_snps(self, file_path=None):
        if file_path is None:
            file_path = self.input_file_pl

        if isinstance(file_path, list):
            for flp in file_path:
                pd.read_csv(flp)
        elif isinstance(file_path, str):
            pd.read_csv(file_path)


class PreapreDataForCircos:
    """Preapre QTL mapping results for the Circos plot.
    """
    def __init__(self, fppl):
        self.fppl = fppl
        self.snps_dfpl = None
        self.snps_dtfm = None

    def _load_snps(self, flpt, **kwargs):
        pass

    @staticmethod
    def _filter(dtfm: pd.DataFrame):
        return dtfm

    @staticmethod
    def _merge(dtfm_pl: list):
        return dtfm_pl

    def prepare(self, **kwargs):
        '''Preapre dataset.'''
        self.snps_dtfm = self._merge([self._filter(self._load_snps(flpt, **kwargs)) for flpt in self.fppl])
        return self

    def write(self, opt_path, **kwargs):
        '''Write the dataset to disk.
        '''
        self.snps_dtfm.to_csv(opt_path, **kwargs)
        return self


if __name__ == '__main__':
    logger.error("Import only module.")
