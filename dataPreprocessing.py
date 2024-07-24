"""This module contains functions for importing and modifying SEC-MS data.
More specifically, this script was based on the output from Spectromine for TMT-18 labelled data."""

#import libraries
import click
import pandas as pd
import numpy as np
import sqlite3
from sklearn.linear_model import LinearRegression
from scipy.signal import filtfilt 
import matplotlib.pyplot as plt




#Functions
class SEC_Data_Manager:
    """A class for managing SEC-MS data. This class contains functions for reshaping, normalizing, and exporting dataframes."""

    def __init__(self, input_df_list, labelling_annotation, mix_norm_annotation, uniprot_gene_table, sec_table=None):
        """Initialize the SEC_Data_Manager object with input dataframes and annotation tables.
        
        Args: 
            input_df_list (list): A list of dataframes containing the SEC-MS data.
            labelling_annotation (pd.DataFrame): A dataframe containing the labelling annotation for the SEC-MS data (see reference for details).
            mix_norm_annotation (pd.DataFrame): A dataframe containing the mix normalization annotation for the SEC-MS data (see reference for details).
            uniprot_gene_table (pd.DataFrame): A dataframe containing the mapping between Uniprot IDs and gene names.
            sec_table (pd.DataFrame): A dataframe containing the SEC molecular weight calibration data. Currently only used for exporting to SECAT. Default is None.
        """

        def read_input_files(input_files):
            """Read the input CSV files and return a list of dataframes.
            
            Args:
                input_files (list): A list of input CSV files.
            
            Returns:
                list: A list of dataframes containing the input data.
            """
            input_df_list = []
            for file in input_files:
                input_df = pd.read_csv(file)

                if 'PP.Phospho (STY)' in input_df.columns:
                    input_df = input_df[~input_df['PP.Phospho (STY)'].isna()]

                input_df_list.append(input_df)
            return input_df_list

        self.input_df_list = read_input_files(input_df_list)
        self.df_list = self.input_df_list
        self.labelling_annotation = pd.read_csv(labelling_annotation, delimiter = "\t")
        self.mix_norm_annotation = pd.read_csv(mix_norm_annotation, delimiter = "\t")
        self.uniprot_gene_table = pd.read_csv(uniprot_gene_table)
        self.sec_table = pd.read_csv(sec_table) if sec_table is not None else None

    def reshape_tmt_df(self, id_col_list, tmt_col_list):
        """Reshape the input dataframes based on the column nomenclature and labelling annotation.
        
        Args:
            id_col_list (list): A list of column names that contain the unique identifiers for the dataframes.
            tmt_col_list (list): A list of column names that contain the TMT channel information for the dataframes.
        
        Returns:
            None (Updates the df_list attribute of the SEC_Data_Manager object).
        """

        #reshaping function for a single dataframe
        def reshape_single_df(input_df, id_columns, tmt_col_name, fraction_annotation):
            """Reshape a single dataframe based on the column nomenclature and labelling annotation.
            
            Args:
                input_df (pd.DataFrame): The input dataframe to be reshaped.
                id_columns (list): A list of column names that contain the unique identifiers for the dataframe.
                tmt_col_name (str): The column name that contains the TMT channel information for the dataframe.
                fraction_annotation (pd.DataFrame): A dataframe containing the fraction annotation for the SEC-MS data. (aka labelling annotation)
                
            Returns:
                pd.DataFrame: The reshaped dataframe with additional columns for experiment, condition, replicate, fraction, and sample.
            """

            #melt df by channel
            tmt_columns = [i for i in input_df.columns if tmt_col_name in i]
            keep_columns = id_columns + tmt_columns
            input_df_drop = input_df[keep_columns].drop_duplicates(keep='first', ignore_index=True)      #input_df.drop(columns=drop_cols, inplace = False).drop_duplicates(keep='first', ignore_index=True)
            input_df_melt = input_df_drop.melt(id_columns, var_name='channel', value_name='intensity')

            #edit column values
            input_df_melt['R.FileName'] = input_df_melt['R.FileName'].str.replace('.raw', '', regex=True)
            input_df_melt['channel'] = input_df_melt['channel'].str.replace(tmt_col_name, '', regex=True)

            #add experiment, condition, and replicate columns
            split_file_name = input_df_melt["R.FileName"].str.split("_", n = 2, expand = True, regex=True)
            input_df_melt["ptm"]= split_file_name[0]
            input_df_melt["replicate"]= split_file_name[1]
            input_df_melt["mix"]= split_file_name[2]
            input_df_melt.drop(columns = ['R.FileName'], inplace = True)

            #add fractions
            input_df_fraction = input_df_melt.merge(fraction_annotation, how="left") #, left_on=['mix', 'channel'], right_on=['mix', 'channels']
            input_df_fraction['sample'] = input_df_fraction['ptm'] + '_' + input_df_fraction['replicate'] + '_' + input_df_fraction['condition'] 
            input_df_fraction = input_df_fraction[~pd.isnull(input_df_fraction['intensity'])].drop_duplicates(keep='first', ignore_index=True)
            input_df_fraction = input_df_fraction.rename(columns={'PG.UniprotIds':'protein_id'})

            if 'PEP.StrippedSequence' in input_df_fraction.columns:
                input_df_fraction = input_df_fraction.rename(columns={'PEP.StrippedSequence':'peptide_id'})

            return input_df_fraction
        
        result_df_list = []
        for i, df in enumerate(self.input_df_list):
            result_df = reshape_single_df(df, id_col_list[i], tmt_col_list[i], self.labelling_annotation)
            result_df_list.append(result_df)
        self.df_list = result_df_list
 
    def prot_id_selection(self, sep=';'):
        """Select a single protein ID from the protein groups based on the intersection of protein IDs between datasets.
        
        Args:
            sep (str): The separator used to split the protein IDs in the protein groups. Default is ';'.
        
        Returns:
            None (Updates the df_list attribute of the SEC_Data_Manager object).
        """
        def split_list_intesect(df_list, sep): #list_1, list_2
            """Get the intersection of protein IDs between datasets.
            
            Args:
                df_list (list): A list of dataframes containing the SEC-MS data.
                sep (str): The separator used to split the protein IDs in the protein groups.
            
            Returns:
                list: A list of common protein IDs between datasets.
            """
            #get common uniprot ID from protein groups
            for i, df in enumerate(df_list):
                list_items = np.unique(df['protein_id'])
                list_splitItems = np.concatenate([s.split(sep) for s in list_items])
                if i == 0:
                    IntersectProtIDs = list_splitItems
                else:
                    IntersectProtIDs = np.intersect1d(IntersectProtIDs, list_splitItems)
            return IntersectProtIDs #list(set(list_1) & set(list_2))

        intersect_proteins = split_list_intesect(self.df_list, sep)

        def get_first_item_intersect(original_series, selected_items, sep):
            """Select a single protein ID from the protein groups based on the intersection of protein IDs between datasets.
            
            Args:
                original_series (pd.series): An array containing the protein IDs from the protein groups.
                selected_items (list): A list of common protein IDs between datasets.
            
            Returns:
                np.array: An array of selected protein IDs.
            """
            result = []
            for item in original_series:
                split_items = item.split(sep)
                result_overlap = [item for item in split_items if item in selected_items]
                if len(result_overlap) == 0:
                    result.append(split_items[0])
                else:
                    result.append(result_overlap[0])
            final_array = np.array(result)
            return final_array
        
        result_df_list = []
        for i, df in enumerate(self.df_list):
            df['protein_id'] = get_first_item_intersect(df['protein_id'], intersect_proteins, sep)
            result_df_list.append(df)
        
        self.df_list = result_df_list

    def mix_normalization(self):
        """Normalize the intensities between TMT mixes based on the mix normalization annotation.
        
        Args: 
            None (Uses the mix_norm_annotation attribute of the SEC_Data_Manager object.)

        Returns:
            None (Updates the df_list attribute of the SEC_Data_Manager object).
        """
        intensity_df_list = self.df_list.copy()
        mix_overlap_table = self.mix_norm_annotation.copy()
        result_df_list = []

        #iterate through each intensity DF in the list of samples
        for i, intensity_df in enumerate(intensity_df_list):

            #iterate through each PTM (hence ..._p) in the mix overlap table (e.g. phospho, global)
            for ptm in np.unique(mix_overlap_table['ptm']):
                mix_overlap_table_p = mix_overlap_table[mix_overlap_table['ptm'] == ptm]
                intensity_df_p = intensity_df[intensity_df['ptm'] == ptm]

                #iterate through each condition (hence ..._pc) in the mix overlap table (e.g. HEK, HCT)
                for cond in np.unique(mix_overlap_table_p['condition']):
                    mix_overlap_table_pc = mix_overlap_table_p[mix_overlap_table_p['condition'] == cond]
                    intensity_df_pc = intensity_df_p[intensity_df_p['condition'] == cond]
                    
                    #iterate through each replicate (hence ..._pcr) in the mix overlap table (e.g. rep01, rep02)
                    for rep in np.unique(mix_overlap_table_pc['replicate']):
                        mix_overlap_table_pcr = mix_overlap_table_pc[mix_overlap_table_pc['replicate'] == rep]
                        intensity_df_pcr = intensity_df_pc[intensity_df_pc['replicate'] == rep]
                        
                        #iterate through each mix in the mix overlap table (e.g. mix01, mix02)
                        for index, row in mix_overlap_table_pcr.iterrows():
                            #get the standard and normalization mix (will normalize the normalization mix to the standard mix)
                            intensity_df_pcr_std = intensity_df_pcr[intensity_df_pcr['mix'] == row['std_mix']]
                            intensity_df_pcr_norm = intensity_df_pcr[intensity_df_pcr['mix'] == row['norm_mix']]
                            
                            #get the median intensity for each fraction in the standard and normalization mix
                            fraction_intensity_pcr_std = intensity_df_pcr_std[['fraction', 'intensity']].groupby(['fraction']).median()
                            fraction_intensity_pcr_norm = intensity_df_pcr_norm[['fraction', 'intensity']].groupby(['fraction']).median()
                            
                            #calculate the normalization factor as a mean of the ratio of the median intensities for each overlapping fraction
                            fraction_intensity_median = fraction_intensity_pcr_std.merge(fraction_intensity_pcr_norm, how='inner', on='fraction', suffixes=('_std', '_norm'))
                            fraction_intensity_median['std_norm_ratio'] = fraction_intensity_median['intensity_std'] / fraction_intensity_median['intensity_norm']
                            normalization_factor = np.mean(fraction_intensity_median['std_norm_ratio'])
                            

                            intensity_df_pcr.loc[intensity_df_pcr_norm.index, 'intensity'] = intensity_df_pcr_norm['intensity'] * normalization_factor

                        intensity_df_pc.loc[intensity_df_pcr.index, 'intensity'] = intensity_df_pcr['intensity']
                        
                    intensity_df_p.loc[intensity_df_pc.index, 'intensity'] = intensity_df_pc['intensity']

                intensity_df.loc[intensity_df_p.index, 'intensity'] = intensity_df_p['intensity']

            result_df_list.append(intensity_df)
        self.df_list = result_df_list

    def collapse_on_fraction(self):
        """Collapse the intensities by averaging the fractions that were measured more than once.
        
        Args:
            None (Uses the df_list attribute of the SEC_Data_Manager object.)
        
        Returns:
            None (Updates the df_list attribute of the SEC_Data_Manager object.)"""
        intensity_df_list = self.df_list.copy()
        result_df_list = []

        #iterate through each intensity DF in the list of samples
        for i, intensity_df in enumerate(intensity_df_list):
            mix_channel_df = intensity_df.copy()
            fraction_df = mix_channel_df.drop(columns=['mix', 'channel'])
            id_cols = list(fraction_df.drop(columns='intensity').columns)

            #take the mean of the fractions that were measured more than once
            fraction_df = fraction_df.groupby(id_cols, as_index=False).mean()
            result_df_list.append(fraction_df)
        self.df_list = result_df_list

    def median_normalization(self, category):
        """Median normalize the intensities between conditions or PTMs.

        Args:
            category (str): The category to normalize the intensities between. Can be 'condition' or 'ptm'.
        
        Returns:
            None (Updates the df_list attribute of the SEC_Data_Manager object.)
        """
        intensity_df_list = self.df_list.copy()
        result_df_list = []

        #iterate through each intensity DF in the list of samples
        for i, input_df in enumerate(intensity_df_list):
            result_df = input_df.copy()

            #iterate through each constant (condition or ptm), and use other as the category to normalize
            if category == 'ptm':
                constant = 'condition'
            else:
                constant = 'ptm'
            for c in np.unique(result_df[constant]):
                result_df_c = result_df[result_df[constant] == c]
                
                #iterate through each replicate
                replicates = np.unique(result_df['replicate'])
                for r in replicates:
                    result_df_cr = result_df_c[result_df_c['replicate'] == r]

                    #iterate through each category (condition or PTM) to normalize
                    categories = np.unique(result_df_cr[category])
                    for i in range(0, len(categories)-1):
                        std_cond = categories[i]
                        norm_cond = categories[i+1]

                        #get the median intensity for the standard and normalization condition
                        std_df = result_df_cr[result_df_cr[category] == std_cond]
                        std_median_intensity = np.nanmedian(std_df['intensity'])

                        norm_df = result_df_cr[result_df_cr[category] == norm_cond]
                        norm_median_intensity = np.nanmedian(norm_df['intensity'])

                        #normalize the normalization condition to the standard condition
                        result_df_cr.loc[norm_df.index, 'intensity'] = (norm_df['intensity'] * (std_median_intensity/norm_median_intensity))
                        #print(c + '-' + r + ' ' + ' ' + norm_cond + 'norm factor: ' + str(round(std_median_intensity/norm_median_intensity), 2))
                    
                    result_df_c.loc[result_df_cr.index, 'intensity'] = result_df_cr['intensity']

                result_df.loc[result_df_c.index, 'intensity'] = result_df_c['intensity']
        
            result_df_list.append(result_df)
        self.df_list = result_df_list

    def combine_df_list(self, average_reps=False, two_rep_mode=False):
        """Combine the dataframes from the list of samples into a single dataframe.
        
        Args:
            two_rep_mode (bool): A boolean indicating whether to combine the dataframes from two replicates (a specific exception made for our data). Default is False.
        
        Returns:
            None (Updates the df attribute of the SEC_Data_Manager object.)
        """
        combined_df = pd.concat(self.df_list).reset_index(drop=True)
        
        if two_rep_mode == True:
            combined_df = combined_df[combined_df['replicate'] != 'rep03'].reset_index(drop=True)
            combined_df = combined_df.copy().replace('rep23', 'rep02')
            combined_df['sample'] = combined_df['ptm'] + '_' + combined_df['replicate'] + '_' + combined_df['condition']

        if combined_df.columns.isin(['peptide_id']).any():
            combined_df['peptide_id'] = combined_df['peptide_id'].replace(np.nan, '0')


        if average_reps == True:
            print('Averaging replicates...')
            group_cols = ['protein_id', 'peptide_id', 'sample', 'replicate', 'ptm', 'condition', 'fraction']
            if 'peptide_id' not in combined_df.columns:
                group_cols.remove('peptide_id')
            combined_df['sample'] = combined_df['ptm'] + '_' + combined_df['condition']
            combined_df['replicate'] = 'mean'
            combined_df = combined_df.groupby(group_cols)['intensity'].mean().reset_index()

        self.df = combined_df

    def apply_smoothing(self, id_cols, input_wide=False, return_wide=False):
        """Apply smoothing to the intensities using a low-pass filter (see reference). 
        The filtfilt linear digital filter was chosen over the more commonly used LOESS method for the following reasons:
            1. Preserving the integrity of signal features without introducing phase shifts.
            2. To avoid the potential for overfitting (we found that LOESS may smooth out smaller features or introduce artifacts).
            3. More computationally efficient for large datasets.
        
        It should be noted, however, that LOESS may be more appropriate given the complex non-linear nature of the data. 
        However, we found that rigorous fine-tuning of the smoothing parameters was required to avoid overfitting.
        
        
        Args:
            id_cols (list): A list of column names that contain the unique identifiers for the dataframe.
            input_wide (bool): A boolean indicating whether the input dataframe is in wide format (see reference). Default is False.
            return_wide (bool): A boolean indicating whether to return the dataframe in wide format (see reference). Default is False.
        
        Returns:
            None (Updates the df attribute of the SEC_Data_Manager object.)
        
        Reference:
            - https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.filtfilt.html
            - wide format: fractions as individual columns.
            - long format: fractions as a single column (often called 'fraction').
        """
        input_df = self.df.copy()
        if input_wide == True:
            long_df = pd.melt(input_df.copy(), id_vars=id_cols, var_name='fraction', value_name='intensity')
        else:
            long_df = input_df.copy()

        long_df['fraction'] = long_df['fraction'].astype(int)
        wide_df = pd.pivot_table(long_df.copy(), index=id_cols, columns='fraction', values='intensity', fill_value=0).reset_index()

        #get intensity column names
        intensity_columns = wide_df.copy()[np.unique(long_df['fraction'])].columns

        #smoothing function
        def smooth_intensities(row_intensities): #check code#
            row_intensities = list(row_intensities)
            #the higher the integer smoothing factor, the more intense the smoothing (2 is a good starting point)
            smoothing_factor = 2
            b = [1.0 / smoothing_factor] * smoothing_factor
            a = 1
            y = filtfilt(b, a, row_intensities)
            return pd.Series(y)

        
        #apply smoothing function to each entry
        wide_df[intensity_columns]  = wide_df.copy()[intensity_columns].apply(smooth_intensities, axis=1)

        if return_wide == True:
            self.df = wide_df
        else:
            long_df_norm = pd.melt(wide_df, id_vars=id_cols, var_name='fraction', value_name='intensity')
            long_df_norm = long_df_norm[long_df_norm['intensity'] != 0].reset_index(drop=True)
            self.df = long_df_norm

    def add_gene_names(self, id_col='protein_id', gene_col='Gene.Names', save_missing=False):
        """Add gene names to the dataframe based on the mapping between Uniprot IDs and gene names.

        Args:
            gene_name_df (pd.DataFrame): A dataframe containing the mapping between Uniprot IDs and gene names.
            id_col (str): The column name that contains the unique identifiers for the dataframe. Default is 'protein_id'.
            gene_col (str): The column name that contains the gene names. Default is 'Gene.Names'.
            save_missing (bool): A boolean indicating whether to save the missing gene names to a file (for troubleshooting and filling in data). Default is False.
        
        Returns:
            None (Updates the df attribute of the SEC_Data_Manager object.)
        """
        output_df = self.df.copy()
        gene_name_df = self.uniprot_gene_table.copy()

        #get complimentary gene names and insert column next to id_col
        gene_values = output_df.merge(gene_name_df, left_on=id_col, right_on='Entry', how='left')[gene_col]
        output_df.insert(output_df.columns.get_loc(id_col) + 1, "gene", gene_values, allow_duplicates=True)

        if save_missing == True:
            pd.Series(pd.unique(output_df[output_df['gene'].isnull()][id_col])).to_csv('missing_genes_list.tsv', sep='\t', index=False)

        self.df = output_df
    
    def reorder_columns(self, col_order):
        """A convenience function to reorder the columns in the dataframe based on the input list of column names.
        
        Args:
            col_order (list): A list of column names in the desired order.
        
        Returns:
            None (Updates the df attribute of the SEC_Data_Manager object.)
        """
        output_df = self.df.copy()
        col_order_present = [value for value in col_order if value in output_df.columns]
        output_df = output_df[col_order_present]
        col_order_present.remove('intensity')
        output_df = output_df.sort_values(by=col_order_present, axis=0, ascending=True, ignore_index=True)
        self.df = output_df
    
    def export_to_sql(self, filename, table_name):
        """Export the dataframe to an SQLite database.
        
        Args:
            filename (str): The filename for the SQLite database.
            table_name (str): The table name for the dataframe.
        
        Returns:
            None (Exports the dataframe to an SQLite database.)"""
        conn = sqlite3.connect(filename)
        self.df.to_sql(table_name, conn, index=False, if_exists='replace')
        conn.close()

    def export_to_csv(self, filename):
        """A completely unnecessary function to export the dataframe to a CSV file.
        
        Args:
            filename (str): The filename for the CSV file.
        
        Returns:
            None (Exports the dataframe to a CSV file.)"""
        self.df.to_csv(filename, index=False)

    def export_for_secat(self, intensity_filename, sec_filename, fraction_col='Fraction', first_fraction=0, combination_factor=1, plot_model=False):
        """Export the dataframe to a format compatible with SECAT (see reference).

        Args:
            intensity_filename (str): The filename for the intensity table.
            sec_filename (str): The filename for the SEC table.
            fraction_col (str): The column name that contains the fraction information. Default is 'Fraction'.
            first_fraction (int): The first fraction number (often SEC fractions are renumbered between standard and MS analysis). Default is 0.
            combination_factor (int): The factor to combine fractions (e.g. 2 for combining each two fractions). Default is 1.
            plot_model (bool): A boolean indicating whether to plot the molecular weight calibration model. Default is False.
            
        Returns:
            None (Exports the dataframes compatible with SECAT - intensity (tsv) and SEC annotation (csv).)
        
        Reference:
            - https://github.com/grosenberger/secat
        """

        #make and save the intensity table
        secat_intensity_table = pd.DataFrame({
            'run_id': self.df['condition'].str.cat(self.df['replicate'], sep = "_").str.cat(self.df['fraction'].astype('string'), sep = "_"),
            'protein_id': self.df['protein_id'],
            'peptide_id': self.df['peptide_id'],
            'peptide_intensity': self.df['intensity']
        }).drop_duplicates(ignore_index=True).sort_values(by=['protein_id', 'peptide_id', 'run_id'], axis=0, ignore_index=True)

        secat_intensity_table.to_csv(intensity_filename, index=False, sep='\t')

        def modify_sec_standard(input_sec_table, fraction_col, first_fraction=0, combination_factor=1):
            """Modify the SEC standard table based on the first fraction and combination factor.
            
            Args:
                input_sec_table (pd.DataFrame): The input SEC standard table.
                fraction_col (str): The column name that contains the fraction information.
                first_fraction (int): The first fraction number (often SEC fractions are renumbered between standard and MS analysis). Default is 0.
                combination_factor (int): The factor to combine fractions (e.g. 2 for combining each two fractions). Default is 1.
            
            Returns:
                pd.DataFrame: The modified SEC standard table.
            """
            sec_table_modified = input_sec_table.copy()
            sec_table_modified[fraction_col] = sec_table_modified[fraction_col] - first_fraction
            sec_table_modified[fraction_col] = sec_table_modified[fraction_col] / combination_factor
            return sec_table_modified
        
        sec_table = self.sec_table.copy()

        if (first_fraction != 0) or (combination_factor != 1):
            sec_table = modify_sec_standard(sec_table, fraction_col, first_fraction, combination_factor)


        #make and save the SEC table
        secat_sec_table = pd.DataFrame({
            'run_id': self.df['condition'].str.cat(self.df['replicate'], sep = "_").str.cat(self.df['fraction'].astype('string'), sep = "_"),
            'sec_id': self.df["fraction"],
            'condition_id': self.df["condition"],
            'replicate_id': self.df["replicate"]
        }).drop_duplicates(ignore_index=True).sort_values(by=['condition_id', 'replicate_id', 'sec_id'], axis=0, ignore_index=True)
        
        #MW calibration 
        mw_standard = sec_table.copy()

        #get MW linear model
        x = np.array(mw_standard['Fraction'] - first_fraction).reshape((-1, 1)) 
        y = np.log10(np.array(mw_standard['MW']))
        model = LinearRegression().fit(x, y) 

        if plot_model == True:
            r_sq = model.score(x, y)
            plt.scatter(mw_standard['Fraction'] - first_fraction, np.log10(mw_standard['MW']))
            plt.plot(mw_standard['Fraction'] - first_fraction, model.predict(x),  color='black')
            plt.text(0.05, 0.95, 'R^2 = ' + str(round(r_sq, 3)), horizontalalignment='left', verticalalignment='center', transform=plt.gca().transAxes)
            plt.xlabel('Fraction')
            plt.ylabel('Log10 Molecular Weight (kDa)')
            plt.title('Molecular Weight Calibration')
            plt.savefig(sec_filename + '_MWcalibration.pdf', bbox_inches='tight', format='pdf')
            plt.close()


        sec_fraction_original = np.array(secat_sec_table['sec_id']).reshape((-1, 1)) 
        secat_sec_table['sec_mw'] = np.power(10, model.predict(sec_fraction_original))
        
        secat_sec_table.to_csv(sec_filename, index=False)


##Begin data processing##
#click command line interface
@click.command()
@click.option('--input_files', '-i', multiple=True, default=['sample_input/hekhct_global_2reps_20230815.csv', 'sample_input/hekhct_phospho_2reps_20231106.csv'], help='The input file(s) containing the SEC-MS data. Default is from sample data.')
@click.option('--labelling_annotation', '-l', default='sample_input/TMTinfo_HEKHCT_TMT18.txt', type=click.Path(exists=True), help='The labelling annotation file for the SEC-MS data.')
@click.option('--mix_norm_annotation', '-m', default='sample_input/TMTmixOverlap_HEKHCT_TMT18.txt', type=click.Path(exists=True), help='The mix normalization annotation file for the SEC-MS data.')
@click.option('--uniprot_gene_table', '-u', default='sample_input/uniprotIDtoGENE_human_allGenes_20230719.csv', type=click.Path(exists=True), help='The Uniprot ID to gene name mapping table.')
@click.option('--sec_table', '-s', default='sample_input/SEC_MW_input_hekhct.txt', type=click.Path(exists=True), help='The SEC molecular weight calibration table.')
@click.option('--output_sql', '-o', required=True, help='The output file for the processed SEC-MS data.')
@click.option('--output_csv', '-c', default=None, help='The output file for the processed SEC-MS data.')
@click.option('--output_intensity', '-int', default=None, help='The output file for the intensity table compatible with SECAT.')
@click.option('--output_sec', '-sec', default=None, help='The output file for the SEC table compatible with SECAT.')
@click.option('--fraction_col', '-f', default='Fraction', help='The column name that contains the fraction information.')
@click.option('--first_fraction', '-ff', default=0, help='The first fraction number (often SEC fractions are renumbered between standard and MS analysis).')
@click.option('--combination_factor', '-cf', default=1, help='The factor to combine fractions (e.g. 2 for combining each two fractions).')
@click.option('--plot_model', '-p', default=False, help='A boolean indicating whether to plot the molecular weight calibration model.')
@click.option('--avg_reps', '-ar', is_flag=True, help='A boolean indicating whether to average replicates.')
def main(input_files, labelling_annotation, mix_norm_annotation, uniprot_gene_table, sec_table, output_sql, output_csv, output_intensity, output_sec, fraction_col, first_fraction, combination_factor, plot_model, avg_reps):
    """Main function to perform SEC-MS data preprocessing."""

    #1. Create SEC_Data_Manager object - an object that contains all the dataframes and information
    click.echo('Creating SEC_Data_Manager object...')
    data_manager = SEC_Data_Manager(input_files, 
                                    labelling_annotation, 
                                    mix_norm_annotation, 
                                    uniprot_gene_table, 
                                    sec_table)

    
    #2. Reshape dataframe based on input column nomenclature 
    click.echo('Reshaping TMT data...')
    data_manager.reshape_tmt_df([['R.FileName', 'PG.UniprotIds'], 
                                ['R.FileName', 'PG.UniprotIds', 'PEP.StrippedSequence']], 
                                ['PG.TMT18_', 'PEP.TMT18_'])

    
    #3. Find intersect of protein groups (separated by ';') between datasets & select a single protein ID
    click.echo('Selecting protein IDs...')
    data_manager.prot_id_selection()

    #4. Perform normalization between TMT mixes, collapsing by average of fractions, and median normalize between conditions (or other categorical columns)
    click.echo('Normalizing intensities...')
    data_manager.mix_normalization()
    data_manager.collapse_on_fraction()
    data_manager.median_normalization('condition')

    #5. Combine all data into a single dataframe
    click.echo('Combining data...')
    data_manager.combine_df_list(average_reps=avg_reps) #two_rep_mode=True

    #6. Add gene names - readding in ensures that all common proteins have the same gene names
    click.echo('Adding gene names...')
    data_manager.add_gene_names(id_col='protein_id', 
                                gene_col='Gene.Names', 
                                save_missing=False)

    #7. Reorder columns - according to the 'col_order' list
    click.echo('Reordering columns...')
    col_order = ['sample', 'ptm', 'condition', 'replicate', 'protein_id', 'gene', 'peptide_id', 'PTM.ModificationTitle', 'PTM.Multiplicity', 'PTM.SiteAA', 'PTM.SiteLocation', 'PTM.SiteProbability', 'fraction', 'intensity']
    data_manager.reorder_columns(col_order)

    #8. export to sql
    click.echo('Exporting to SQL...')
    data_manager.export_to_sql(output_sql, 'initial_intensity_table')


    #9. export to csv
    if output_csv is not None:
        click.echo('Exporting to CSV...')
        data_manager.export_to_csv(output_csv)


    #10. export for secat
    if output_intensity is not None and output_sec is not None:
        click.echo('Exporting for SECAT...')
        data_manager.export_for_secat(output_intensity, 
                                    output_sec, 
                                    fraction_col=fraction_col, 
                                    first_fraction=first_fraction, 
                                    combination_factor=combination_factor, 
                                    plot_model=plot_model)

if __name__ == '__main__':
    main()