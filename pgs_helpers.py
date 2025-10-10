import pandas as pd
import numpy as np


def split_rows_with_multiple_alleles(df, samples):
    # Create an empty list to store the new rows
    new_rows = []
    # Iterate through each row
    for index, row in df.iterrows():
        alt_alleles = row['alt'].split(',')
        # Check if there are multiple alleles in the ALT field
        if len(alt_alleles) > 1:
            for allele_num, allele in enumerate(alt_alleles):
                # Create a new row for each allele
                new_row = row.copy()
                new_row['alt'] = allele
                # Update genotype fields
                for col in samples:
                    genotype = row[col]
                    # Split the genotype and process it
                    if genotype != './.':
                        gt_alleles = genotype.split('/')
                        new_gt = ['0' if (int(gt) != allele_num + 1 and gt != '0') else gt for gt in gt_alleles]
                        new_row[col] = '/'.join(new_gt)
                new_rows.append(new_row)
        else:
            new_rows.append(row)
    
    new_df = pd.DataFrame(new_rows).reset_index(drop=True)
    return new_df


def convert_genotype_to_alt_allele_count(df, samples):
    nalt_df = df.copy()
    # Iterate through each row
    for index, row in df.iterrows():
        # Update genotype fields
        for col in samples:
                genotype = row[col]
                if genotype != './.':
                    # Split the genotype and count non-zero alleles
                    alleles = genotype.split('/')
                    alt_allele_count = sum([1 for allele in alleles if allele != '0'])
                    nalt_df.at[index, col] = alt_allele_count
                else:
                    nalt_df.at[index, col] = np.nan

    return nalt_df

def calculate_pgs(df_snp, df_eff, df_amp_samples, sig_only=False, standardise=True):
    
    if sig_only:
        selected_snps = df_eff.query("sig == True").index.tolist()
        
    PGS = _calculate_pgs(genotypes=df_snp.values, 
                        odds_ratios=df_eff.odds_ratio, 
                        standardise=standardise)
    
    pgs_df = pd.DataFrame({'sample_id': df_snp.columns,
                           'pgs':PGS})

    return pgs_df

def _calculate_pgs(genotypes, odds_ratios, standardise):
    """
    Calculate Polygenic Scores (PGS) for samples.
    
    Parameters:
    genotypes : ndarray
        3D array of shape (n_snps, n_samples, ploidy)
    odds_ratios : ndarray
        1D array of odds ratios for each SNP
    standardise : bool, optional
        Whether to standardise the final scores
    
    Returns:
    ndarray
        1D array of PGS for each sample
    """
    # Ensure odds_ratios is a 1D array
    odds_ratios = np.asarray(odds_ratios).flatten().round(3)
    
    # Calculate Betas / log ORs (effect sizes)
    effect_sizes = np.log(odds_ratios)
    
    # Reshape effect sizes to broadcast correctly
    effect_sizes = effect_sizes.reshape(-1, 1, 1)
    
    # Calculate PRS
    pgs = np.sum(genotypes * effect_sizes, axis=(0, 1))
    
    # Standardise if requested
    if standardise:
        pgs = (pgs - np.mean(pgs)) / np.std(pgs)
    
    return pgs













################ LD

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def compute_snp_correlations(df, correlation_threshold=0.5):
    corr_matrix = df.corr()
    mask = np.triu(np.ones_like(corr_matrix, dtype=bool))
   
    high_corr_pairs = []
    for i in range(len(corr_matrix.columns)):
        for j in range(i+1, len(corr_matrix.columns)):
            if abs(corr_matrix.iloc[i, j]) > correlation_threshold:
                high_corr_pairs.append((corr_matrix.columns[i], corr_matrix.columns[j], corr_matrix.iloc[i, j]))
    
    print(f"Total SNPs analyzed: {len(df.columns)}")
    print(f"Number of highly correlated SNP pairs (|r| > {correlation_threshold}): {len(high_corr_pairs)}")
    
    return corr_matrix, high_corr_pairs

def plot_correlation_heatmap(corr_matrix):
    plt.figure(figsize=(12, 10))
    sns.heatmap(corr_matrix, annot=False, cmap='coolwarm', vmin=-1, vmax=1, center=0)
    plt.title('SNP Correlation Heatmap')
    plt.show()

def select_representative_snps(df, high_corr_pairs):
    to_remove = set()
    pos_corr_removed = 0
    neg_corr_removed = 0
    
    for snp1, snp2, corr in high_corr_pairs:
        if snp1 not in to_remove and snp2 not in to_remove:
            if corr < 0:  # Negative correlation
                if df.loc[snp1, 'odds_ratio'] < 1:
                    to_remove.add(snp1)
                    neg_corr_removed += 1
                elif df.loc[snp2, 'odds_ratio'] < 1:
                    to_remove.add(snp2)
                    neg_corr_removed += 1
                else:
                    if df.loc[snp1, 'fdr'] < df.loc[snp2, 'fdr']:
                        to_remove.add(snp2)
                        neg_corr_removed += 1
                    else:
                        to_remove.add(snp1)
                        neg_corr_removed += 1
            else:  # Positive correlation
                if df.loc[snp1, 'fdr'] < df.loc[snp2, 'fdr']:
                    to_remove.add(snp2)
                    pos_corr_removed += 1
                else:
                    to_remove.add(snp1)
                    pos_corr_removed += 1
    
    selected_snps = np.array(list(set(df.index) - to_remove))
    
    print(f"Total SNPs before selection: {len(df)}")
    print(f"SNPs removed due to positive correlation: {pos_corr_removed}")
    print(f"SNPs removed due to negative correlation: {neg_corr_removed}")
    print(f"Total SNPs removed: {len(to_remove)}")
    print(f"SNPs retained after selection: {len(selected_snps)}")
    
    return selected_snps

def linked(df_genos, df_eff, correlation_threshold=0.5):
    corr_matrix, high_corr_pairs = compute_snp_correlations(df_genos.filter(like='snp_'), correlation_threshold)
    plot_correlation_heatmap(corr_matrix)
    selected_snps = select_representative_snps(df_eff, high_corr_pairs)
    return np.sort(selected_snps)