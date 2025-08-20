#name: calculateMSR
#description: Calculates Minimum Significant Ratio (MSR) for compounds based on IC50 values, run dates, assay names, and target entities. 
#environment: channels: [Conda-forge], dependencies: [python=3.12, {pip: [statsmodels]}]
#language: python
#tags: Transform
#top-menu: Curves | Calculate MSR
#input: dataframe table [Input data table]
#input: column ic50Column [IC50 column]
#input: column compoundIdColumn [Compound ID column]
#input: column runDateColumn [Run date column]
#input: column assayNameColumn [Assay name column]
#input: column targetEntityColumn [Target entity column]
#output: dataframe res {action:join(table)}

import pandas as pd
import numpy as np
import warnings
from statsmodels.regression.mixed_linear_model import MixedLM
from typing import Dict, Tuple, Optional
from scipy import stats


def prepare_data(df, ic50_col, compound_col, 
                run_date_col, assay_col, target_col):
    """Prepare and clean the data for MSR calculation."""
    
    # Make a copy and select relevant columns
    df_work = df[[ic50_col, compound_col, run_date_col, assay_col, target_col]].copy()
    
    # Remove rows with missing values
    df_work = df_work.dropna()
    # remove rows where reported_ic50 is not a valid float or is zero
    df_work = df_work[df[ic50_col].apply(lambda x: isinstance(x, (int, float)) and x > 0)]
    # Remove non-positive IC50 values
    df_work = df_work[df_work[ic50_col] > 0]
    
    # Calculate log10 IC50
    df_work['log10_ic50'] = np.log10(df_work[ic50_col])
    
    # Convert run_date to datetime if it's not already
    if not pd.api.types.is_datetime64_any_dtype(df_work[run_date_col]):
        df_work[run_date_col] = pd.to_datetime(df_work[run_date_col])
    
    # Create numeric run_date for modeling (days since first run per group)
    min_date = df_work[run_date_col].min()
    df_work['run_date_numeric'] = (df_work[run_date_col] - min_date).dt.days
    
    return df_work

def calculate_msr_mixed_model(group, compound, assay, target,
                             n_runs, n_measurements):
    """Calculate MSR using mixed effects model for compounds with sufficient runs."""
    
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            # For a single compound, we don't need fixed effects for compound
            # We just need to model the variance components
            # Random effect: run_date, Fixed effectercept only
            
            # Create design matrix (just intercept)
            exog = np.ones((len(group), 1))
            
            model = MixedLM(
                endog=group['log10_ic50'].values,
                exog=exog,
                groups=group['run_date_numeric'].values
            )
            
            result = model.fit(method='lbfgs', maxiter=1000)
            # Extract variance components
            residual_variance = result.scale  # Within-run variance
            random_effect_variance = float(result.cov_re[0][0]) if len(result.cov_re) > 0 and len(result.cov_re[0]) > 0 else 0  # Between-run variance
            
            # Calculate s as square root of sum of variance components
            s = np.sqrt(random_effect_variance + residual_variance)
            
            # Calculate MSR = 10^(2*sqrt(2)*s)
            msr = 10**(2 * np.sqrt(2) * s)
            
            return {
                'compound_id': compound,
                'assay_name': assay,
                'target_entity': target,
                'database_msr': msr,
                'residual_variance': residual_variance,
                'run_date_variance': random_effect_variance,
                'total_variance': random_effect_variance + residual_variance,
                's_value': s,
                'n_runs': n_runs,
                'n_measurements': n_measurements,
                'method': 'mixed_model',
                'mean_log10_ic50': group['log10_ic50'].mean(),
                'geometric_mean_ic50': 10**group['log10_ic50'].mean()
            }
            
    except Exception as e:
        print(e)
        # Fallback to simple method if mixed model fails
        return calculate_msr_simple(group, compound, assay, target, n_runs, n_measurements)

def calculate_msr_simple(group, compound, assay, target,
                        n_runs, n_measurements):
    """Calculate MSR using simple variance decomposition for compounds with fewer runs."""
    
    # Calculate run means
    run_means = group.groupby('run_date_numeric')['log10_ic50'].agg(['mean', 'count']).reset_index()
    run_means.columns = ['run_date_numeric', 'run_mean', 'n_per_run']
    
    # Calculate overall mean
    overall_mean = group['log10_ic50'].mean()
    
    # Calculate between-run variance (if multiple runs)
    if len(run_means) > 1:
        between_run_variance = np.var(run_means['run_mean'], ddof=1)
    else:
        between_run_variance = 0
    
    # Calculate within-run variance
    within_run_variances = []
    for run_date in group['run_date_numeric'].unique():
        run_data = group[group['run_date_numeric'] == run_date]['log10_ic50']
        if len(run_data) > 1:
            within_run_variances.append(np.var(run_data, ddof=1))
    
    if within_run_variances:
        within_run_variance = np.mean(within_run_variances)
    else:
        # If no within-run replicates, estimate from total variance
        total_variance = np.var(group['log10_ic50'], ddof=1)
        within_run_variance = max(0, total_variance - between_run_variance)
    
    # Calculate s as square root of sum of variance components
    s = np.sqrt(between_run_variance + within_run_variance)
    
    # Calculate MSR = 10^(2*sqrt(2)*s)
    msr = 10**(2 * np.sqrt(2) * s)
    
    return {
        'compound_id': compound,
        'assay_name': assay,
        'target_entity': target,
        'database_msr': msr,
        'residual_variance': within_run_variance,
        'run_date_variance': between_run_variance,
        'total_variance': between_run_variance + within_run_variance,
        's_value': s,
        'n_runs': n_runs,
        'n_measurements': n_measurements,
        'method': 'simple_variance',
        'mean_log10_ic50': overall_mean,
        'geometric_mean_ic50': 10**overall_mean
    }

def calculate_msr_for_compound(group, ic50_col, compound_col, 
                              run_date_col, min_runs, min_measurements,
                              compound, assay, target):
    """Calculate MSR for a specific compound/assay/target combination."""
    
    # Check minimum requirements
    n_runs = group['run_date_numeric'].nunique()
    n_measurements = len(group)
    
    if n_measurements < min_measurements:
        return None
        
    if n_runs < 2:  # Need at least 2 different runs to estimate between-run variance
        return None
    
    try:
        # Method 1: Use mixed effects model if we have enough runs
        if n_runs >= min_runs:
            return calculate_msr_mixed_model(group, compound, assay, target, n_runs, n_measurements)
        else:
            # Method 2: Simple variance decomposition for fewer runs
            return calculate_msr_simple(group, compound, assay, target, n_runs, n_measurements)
            
    except Exception as e:
        print("Model fitting failed for {compound}/{assay}/{target}: {ee}".format(
            compound=compound, assay=assay, target=target, ee=str(e)
        ))
        # If both methods fail, return None
        return None

def calculate_compound_database_msr(df, 
                                   ic50_col = 'reported_ic50',
                                   compound_col = 'compound_id',
                                   run_date_col = 'run_date',
                                   assay_col = 'assay_name',
                                   target_col = 'target_entity',
                                   min_runs = 6,
                                   min_measurements = 2):
    """
    Calculate Database MSR for each unique combination of compound_id, assay_name, and target_entity.
    Each compound must have multiple measurements across different runs.
    
    Parameters:
    -----------
    df 
        Input dataframe with assay data
    ic50_col 
        Column name for IC50 values
    compound_col 
        Column name for compound identifiers
    run_date_col 
        Column name for run dates
    assay_col 
        Column name for assay names
    target_col 
        Column name for target entities
    min_runs 
        Minimum number of different runs required for a compound
    min_measurements 
        Minimum number of measurements required for a compound
        
    Returns:
    --------
    pd.DataFrame
        DataFrame with MSR values for each compound/assay/target combination
    """
    
    # Prepare the data
    df_clean = prepare_data(df, ic50_col, compound_col, run_date_col, assay_col, target_col)
    
    # Group by compound, assay, and target
    results = []
    
    for (compound, assay, target), group in df_clean.groupby([compound_col, assay_col, target_col]):
        try:
            msr_result = calculate_msr_for_compound(
                group, ic50_col, compound_col, run_date_col, 
                min_runs, min_measurements, compound, assay, target
            )
            if msr_result is not None:
                results.append(msr_result)
        except Exception as e:
            print("Error calculating MSR for {compound}/{assay}/{target}: {ee}".format(
                compound=compound, assay=assay, target=target, ee=str(e)
            ))
            # Continue to next group if there's an error
            continue
    
    if not results:
        return pd.DataFrame()
    
    return pd.DataFrame(results)

def calculate_control_compound_msr(df, control_compounds,
                                  ic50_col = 'reported_ic50',
                                  compound_col = 'compound_id',
                                  run_date_col = 'run_date',
                                  assay_col = 'assay_name',
                                  target_col = 'target_entity'):
    """
    Calculate Control Compound MSR for specified control compounds.
    This is useful for assay performance monitoring.
    """
    
    control_df = df[df[compound_col].isin(control_compounds)]
    return calculate_compound_database_msr(
        control_df, ic50_col, compound_col, run_date_col, assay_col, target_col,
        min_runs=6, min_measurements=6
    )

def generate_msr_summary(msr_results):
    """Generate a summary of MSR results with additional statistics."""
    
    if msr_results.empty:
        return msr_results
    
    summary = msr_results.copy()
    
    # Add interpretation categories
    summary['msr_category'] = pd.cut(
        summary['database_msr'], 
        bins=[0, 3, 5, 10, float('inf')], 
        labels=['Excellent (<3)', 'Good (3-5)', 'Fair (5-10)', 'Poor (>10)']
    )
    
    # Add percentage of variance from each source
    summary['pct_run_date_variance'] = (
        summary['run_date_variance'] / summary['total_variance'] * 100
    ).fillna(0)
    summary['pct_residual_variance'] = (
        summary['residual_variance'] / summary['total_variance'] * 100
    ).fillna(0)
    
    # Add 95% confidence interval for the compound potency
    # CI width = MSR gives the range where we expect 95% of measurements
    summary['ic50_lower_95ci'] = summary['geometric_mean_ic50'] / summary['database_msr']
    summary['ic50_upper_95ci'] = summary['geometric_mean_ic50'] * summary['database_msr']
    
    return summary

def get_compounds_for_msr(df, 
                         compound_col = 'compound_id',
                         run_date_col = 'run_date',
                         assay_col = 'assay_name',
                         target_col = 'target_entity',
                         min_measurements = 2):
    """
    Identify which compounds have sufficient data for MSR calculation.
    """
    
    # Count measurements per compound/assay/target combination
    compound_counts = (df.groupby([compound_col, assay_col, target_col])
                      .agg({
                          run_date_col: 'nunique',
                          compound_col: 'count'  # total measurements
                      })
                      .rename(columns={
                          run_date_col: 'n_runs',
                          compound_col: 'n_measurements'
                      })
                      .reset_index())
    
    # Filter compounds with sufficient data
    sufficient_data = compound_counts[
        (compound_counts['n_measurements'] >= min_measurements) & 
        (compound_counts['n_runs'] >= 2)
    ]
    
    return sufficient_data.sort_values('n_measurements', ascending=False)


np.random.seed(42)
odf = table
# odf = pd.read_csv("assay.csv")
# compoundIdColumn = 'compound_id'
# runDateColumn = 'run_date'
# assayNameColumn = 'assay_name'
# targetEntityColumn = 'cell_line'
# ic50Column = 'reported_ic50'
sufficient_compounds = get_compounds_for_msr(odf, compound_col=compoundIdColumn, run_date_col=runDateColumn, assay_col=assayNameColumn, target_col=targetEntityColumn)
msr_results = calculate_compound_database_msr(odf, min_runs=2, min_measurements=2, ic50_col=ic50Column, compound_col=compoundIdColumn, run_date_col=runDateColumn, assay_col=assayNameColumn, target_col=targetEntityColumn)

# generate a dataset which is the same length as the input data, and for each corresponding row, have the results from the msr_results dataset except the compound_id, assay_name and target_entity. use these three for matching between the two
# res = msr_results

msr_results = msr_results.rename(columns=lambda x: f"{x} (msr)" if x not in ['compound_id', 'assay_name', 'target_entity'] else x) # give postfix
# Merge the results back to the original dataframe
res = odf.merge(
    msr_results,
    left_on=[compoundIdColumn, assayNameColumn, targetEntityColumn],
    right_on=['compound_id', 'assay_name', 'target_entity'],
    how='left',
    suffixes=('', ' (msr)')
)

orig_columns_names = odf.columns.tolist()

res = res.drop(columns=orig_columns_names)