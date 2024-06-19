import vcfETLExceptions
from vcfETL import vcfETL

import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import r
import numpy as np
from scipy.stats import chi2_contingency
import xlsxwriter
import pickle


def get_called_samples_gt(samples):
    """
    Returns a datafram with the sample ID and its genotype for the samples for which there is calling. Not valid for
    multiallelic loci.

    """

    def get_sample_gt(sdata):
        """ Returns a sample's genotype. If there is no call, it returns nocall; if the genotype is not valid, it returns no_valid_gt). """
        
        if sdata.get('GT')[0].startswith('.') or sdata.get('DP')[0] == 0:
            return 'nocall'
        else:
            alelo1 = sdata.get('GT')[0]
            alelo2 = sdata.get('GT')[2]

            if alelo1 == '0' and alelo2 == '0':
                return 0
            elif (alelo1 == '0' and alelo2 == '1') or (alelo1 == '1' and alelo2 == '0'):
                return 1
            elif alelo1 == '1' and alelo2 == '1':
                return 2
            else:
                raise vcfETLExceptions.UnexpectedGenotypeError(sdata.get('GT'))
                return 'no_valid_gt'

    return pd.DataFrame({sample: gt for sample, sdata in samples.items()
                         if (gt := get_sample_gt(sdata)) not in ('nocall', 'no_valid_gt')}.items(),
                        columns=['ID', 'GT'])


def filter_dataframe(df, headers, values):
    """
    Filters a dataframe based on the values specified in the corresponding columns.

    Args:
        df (pd.DataFrame): original DataFrame.
        headers (list): list of header names of the columns that must be filtered.
        values (list): list of values that must coincide in the correspoding columns.

    Returns:
        pd.DataFrame: filtered DataFrame.
    """

    filtered_df = df.copy()  # Make a copy of the original DataFrame

    for header, value in zip(headers, values):
        if isinstance(value, list) or isinstance(value, tuple):
            value_mask = filtered_df[header].isin(value)
        else:
            value_mask = (filtered_df[header] == value)

        filtered_df = filtered_df.loc[value_mask]

    return filtered_df


def segregating_families(df, casos=None, controles=None):
    segregating_families = []
    non_segregating_families = []

    for family, group in df.groupby('FAMILIA'):
        casos_gt = group[group['GRUPO'].isin(casos)]['GT'].unique()
        controles_gt = group[group['GRUPO'].isin(controles)]['GT'].unique()
        # If there is at least one case and one control in the family
        if len(casos_gt) >= 1 and len(controles_gt) >= 1:
            # If there are no common genotypes between cases and controls
            if not set(casos_gt).intersection(controles_gt):
                segregating_families.append(family)
            else:
                non_segregating_families.append(family)

    return segregating_families, non_segregating_families


def at_least_2x2(observed_matrix):
    """
    Check if the contingency matrix has at least 2 rows and 2 columns
    """
    return observed_matrix.shape[0] > 1 and observed_matrix.shape[1] > 1


def do_i_fisher_test(observed_table):
    # Verify the sample size in each cell
    min_expected_count = 5
    return any(observed_table < min_expected_count)


def do_fisher_test(observed_matrix):
    # Convert the pandas DataFrame to an R matrix
    observed_matrix_r = r.matrix(robjects.FloatVector(observed_matrix.values.ravel()),
                                 nrow=observed_matrix.shape[0], ncol=observed_matrix.shape[1])

    # Convert row and column names to strings
    rownames = robjects.StrVector(list(map(str, observed_matrix.index)))
    colnames = robjects.StrVector(list(map(str, observed_matrix.columns)))

    # Assign row and column names to the R matrix
    robjects.r.assign("observed_misc_vs_kd_r", observed_matrix_r)
    robjects.r("rownames(observed_misc_vs_kd_r) <- %s" % rownames.r_repr())
    robjects.r("colnames(observed_misc_vs_kd_r) <- %s" % colnames.r_repr())

    # Perform Fisher's Exact Test using R's fisher.test function
    return stats.fisher_test(robjects.r["observed_misc_vs_kd_r"])


def do_chi_square_test(observed_matrix):
    # Perform chi-square test using scipy's chi2_contingency function
    chi2, p_value, dof, expected = chi2_contingency(observed_matrix)

    return chi2, p_value

def read_gwas_catalog(path):
    # Use a context manager to open the file you want to read from
    with open(path, 'rb') as f:
        # Use pickle.load() to read the dictionary from the file
        data = pickle.load(f)
    return data

if __name__ == '__main__':

    gwas_catalog_pickle = 'Path/To/File/gwas_catalog.pickle'
    with open(gwas_catalog_pickle, "rb") as file:
        gwas_catalog = pickle.load(file)

    out_vars = dict()

    # Import R's "stats" package
    stats = importr("stats")

    try:
        # Load GWAS Catalog data
        gwas_catalog = pickle.load(
            open('Path/To/File/gwas_catalog_selected_dedup.pickle', 'rb'))
        # Load phenotype data
        phenodata = pd.read_csv('Path/To/File/pheno_data.csv')
        # Load VCF data
        vcf = vcfETL('Path/To/File/CoKid.norm.vcf.gz',
                     normalized=True)
        
        for var_id, variant in vcf.parse_line():
            # Consider only biallelic variants
            if len(variant['ALT']) == 1:
                snp_id_vcf = f"{variant['CHROM']}:{variant['POS']}"
                if traits := gwas_catalog.get(snp_id_vcf, False):
                    # It is a variant with traits of interest
                    samples_called = get_called_samples_gt(variant['SAMPLES'])
                    df = pd.merge(samples_called, phenodata.dataset, on=['ID'])
                    df_misc_vs_kd = filter_dataframe(df, ['GRUPO'], [('KD', 'MISC')])
                    
                    observed_misc_vs_kd = pd.crosstab(df_misc_vs_kd.GRUPO, df_misc_vs_kd.GT)
                    if at_least_2x2(observed_misc_vs_kd):
                        if do_i_fisher_test(observed_misc_vs_kd):
                            result = do_fisher_test(observed_misc_vs_kd)
                            test = 'fisher'
                            p_value = result[0][0]
                        else:
                            test = 'chi2'
                            chi2, p_value = do_chi_square_test(observed_misc_vs_kd)

                        # Genotype and phenotype (KD or MIS-C) are not independent
                        if p_value < 0.05:
                            segregating_kd, no_segregating_kd = segregating_families(df, controles=['C'], casos=['KD'])
                            segregating_misc, no_segregating_misc = segregating_families(df, controles=['C'],
                                                                                         casos=['MISC'])

                            # Generate one output line per each trait
                            for trait, gwas_data in traits.items():
                                out_vars[var_id] = {
                                    'var_id': var_id,
                                    'dbSNP': ','.join(str(x) for x in variant.get('ID', [])),
                                    'test': test,
                                    'p_value': p_value,
                                    'trait': trait,
                                    'general_category': gwas_data.get('General category', ''),
                                    'trait_group': gwas_data.get('Group', ''),
                                    'trait_disease': gwas_data.get('Disease', ''),
                                    'trait_notes': gwas_data.get('Notes', ''),
                                    'num_segr_kd': len(segregating_kd),
                                    'segr_kd': ','.join(str(x) for x in segregating_kd),
                                    'num_no_segr_kd': len(no_segregating_kd),
                                    'no_segr_kd': ','.join(str(x) for x in no_segregating_kd),
                                    'num_segr_misc': len(segregating_misc),
                                    'segr_misc': ','.join(str(x) for x in segregating_misc),
                                    'num_no_segr_misc': len(no_segregating_misc),
                                    'no_segr_misc': ','.join(str(x) for x in no_segregating_misc),
                                }
                                pass
            else:
                print(f"WARNING: Multiallelic variant {var_id} found. Skipping...")

    except vcfETLExceptions.UnexpectedGenotypeError as e:
        print('Unexpected genotype: {}'.format(e))
    except vcfETLExceptions.NonNormalizedVCFError as e:
        print('Error in the VCF: {}'.format(e))
    except RuntimeError as e:
        print('Error in the Fisher test for {}'.format(var_id))

    # Output: CSV file
    # Convert out_vars dictionary to a pandas DataFrame
    out_vars_df = pd.DataFrame.from_dict(out_vars, orient='index')
    out_vars_df.to_csv('Path/To/File/out_selected_variants_information.csv', index=False)
