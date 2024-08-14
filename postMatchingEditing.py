"""This module is used to edit the alignment/matched table.
This version only incorporates alignment of 2 conditions."""

#import libraries
import pandas as pd
import numpy as np
import os
import sqlite3
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import click
import warnings

class SQLiteMatchingManager:
    """Class to manage the SQLite database containing the matched table."""

    def __init__(self, db_filename):
        """Initialize the SQLiteMatchingManager object.
        
        Args:
            db_filename (str): The path to the SQLite database file.
        """
        self.db_filename = db_filename
        self.db_folder_path = os.path.dirname(db_filename)
        self.conn = sqlite3.connect(db_filename)

    def load_dataframes(self, peak_table_name='peak_table', match_table_name='alignment_table'):
        """Load the matched table from the SQLite database.
        
        Args:
            peak_table_name (str): The name of the table containing the peak table (Default: 'peak_table').
            match_table_name (str): The name of the table containing the matched table (Default: 'alignment_table').
            
        Outputs:
            - the peak table is stored in the 'peak_table' attribute of the object.
            - the matched/alignment table is stored in the 'match_table' attribute of the object.
        """
        self.peak_table = pd.read_sql_query(f"SELECT * FROM {peak_table_name};", self.conn)
        self.match_table = pd.read_sql_query(f"SELECT * FROM {match_table_name};", self.conn)

        def convert_string_to_list(cell):
            """Convert a string to a list if the string contains the '|' character. -> previously used to store multiple values in a single cell."""
            if isinstance(cell, str) and '|' in cell:
                return cell.split('|')
            elif pd.isna(cell):
                return cell
            else:
                return cell

        # Apply the function to each element in the DataFrame
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=FutureWarning)
            self.match_table = self.match_table.applymap(convert_string_to_list)
    
    def save_to_SQL(self, df, table_name):
        """Save a DataFrame to a SQLite database.
        
        Args:
            df (DataFrame): The DataFrame to save.
            table_name (str): The name of the table to save the DataFrame.
        
        Returns:
            None (the DataFrame is saved to the SQLite database)."""
        
        def convert_lists_to_string(cell):
            """Convert a list to a string if the cell contains a list. -> used to store multiple values in a single cell."""
            if isinstance(cell, list):
                return '|'.join(map(str, cell))
            elif pd.isna(cell):
                return cell
            else:
                return str(cell)
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=FutureWarning)
            df = df.applymap(convert_lists_to_string)

        # Register adapters for other non-supported types if needed
        sqlite3.register_adapter(np.int64, lambda val: int(val))
        sqlite3.register_adapter(np.int32, lambda val: int(val))

        # Save to SQLite
        df.to_sql(table_name, self.conn, index=False, if_exists='replace')
        self.conn.commit()

    def close_connection(self):
        """Close the connection to the SQLite database."""
        self.conn.close()


class matchTableEditor:
    """Class to edit the matched table."""

    def __init__(self, match_table, peak_table):
        """Initialize the matchTableEditor object and match_table."""
        self.match_table = match_table
        self.peak_table = peak_table

    def add_mean_apex(self):
        """Updates matched table with a column containing the mean of the aligned peak apex values."""

        def find_mean_value(lst):
            #cleaned_list removes NaN values and converts the list to floats
            cleaned_list = [float(value) for value in lst if not pd.isna(float(value))]
            return np.mean(cleaned_list)
        
        self.match_table['mean_peak_apex'] = self.match_table['peak_apex'].apply(lambda x: find_mean_value(x))

    def add_fold_changes(self):
        """Updates matched table with a column containing the log2 fold changes of the peak heights."""

        def fold_changes(lst):
            #cleaned_list removes NaN values and converts the list to floats
            cleaned_list = [float(value) for value in lst if not pd.isna(value)]

            #calculate fold changes if two values are present. Currently, only two values are supported.
            if len(cleaned_list) == 2:
                return np.log2(cleaned_list[0] / cleaned_list[1])
            else:
                return np.nan
            
        self.match_table['height_log2fx'] = self.match_table['peak_height'].apply(lambda x: fold_changes(x))

    def add_fold_change_eval(self, threshold=np.log2(1.5)):
        """Updates matched table with a column containing the fold change evaluation based on a threshold.
        
        Args:
            threshold (float): The threshold to evaluate the fold change (Default: np.log2(1.5)).
        """

        def fold_change_eval(row, threshold):
            if pd.isna(row['height_log2fx']):
                return None
            elif row['height_log2fx'] > threshold:
                return row['condition'][0]
            elif row['height_log2fx']  < -threshold:
                return row['condition'][1]
            else:
                return 'NS'
        self.match_table['enrichment'] = self.match_table.apply(fold_change_eval, args=(threshold,), axis=1)

    def add_monomer_MW(self, uniprot_df):
        """Updates matched table with a column containing the monomer mass of the protein.
        
        Args:
            uniprot_df (DataFrame): The DataFrame containing the Uniprot data.
        """
        uniprot_merge = uniprot_df.copy()[['Entry', 'Mass']]
        self.match_table = pd.merge(self.match_table.copy(), uniprot_merge.rename({'Entry':'protein_id', 'Mass':'monomer_mass'}, axis=1), how='left', on='protein_id')

    def add_estimated_MW(self, mw_standard, first_fraction, plot=False, plot_filename=None):
        """Updates matched table with a column containing the estimated molecular weight of the protein.
        
        Args:
            mw_standard (DataFrame): The DataFrame containing the molecular weight standard data.
            first_fraction (int): The first fraction number.
            plot (bool): Whether to plot the model (Default: False).
            plot_filename (str): The filename to save the plot (Default: None).
        
        Outputs:
            - MW model plot (if plot=True).
            - Updates the matched table with the estimated molecular weight.
        """
        
        #get MW log-linear model
        x = np.array(mw_standard['Fraction'] - first_fraction).reshape((-1, 1)) 
        y = np.log10(np.array(mw_standard['MW']))
        model = LinearRegression().fit(x, y) 

        #plot model
        if plot:
            r_sq = model.score(x, y)
            plt.scatter(mw_standard['Fraction'] - first_fraction, np.log10(mw_standard['MW']))
            plt.plot(mw_standard['Fraction'] - first_fraction, model.predict(x),  color='black')
            plt.text(0.8, 0.95, 'R^2 = ' + str(round(r_sq, 3)), horizontalalignment='left', verticalalignment='center', transform=plt.gca().transAxes)
            plt.xlabel('Fraction')
            plt.ylabel('Log10 Molecular Weight (kDa)')
            plt.title('Molecular Weight vs Fraction')
            plt.savefig(plot_filename, bbox_inches='tight', format='pdf')
            plt.close()

        #Estimate apex MW
        apex_fractions = np.array(self.match_table['mean_peak_apex']).reshape((-1, 1)) 
        self.match_table['estimated_MW'] = np.power(10, model.predict(apex_fractions))

    def add_region(self, mass_factor=1.5):
        """Updates matched table with a column containing the region classification based on the estimated molecular weight."""
        match_table = self.match_table.copy()
        self.match_table['region'] = np.where(match_table['monomer_mass']*mass_factor > match_table['estimated_MW'], 'monomer', 'complex')

    def classify_col(self, col):
        """Updates matched table with a column evaluating the log2fx values of the input column."""

        def classification_fucntion(row, col):
            #get condition list and evaluation list
            condition_list = [value for value in row['condition']]
            evaluation_list = [float(value) for value in row[col]]

            #get result list
            result_i = []
            for i in range(0, len(condition_list)):
                if evaluation_list[i] > 0:
                    result_i.append(condition_list[i])

            #concat result string
            if len(result_i) == 0:
                result_i = None
            elif len(result_i) == 1:
                result_i = result_i[0]
            else:   
                result_i = 'Both'
            return result_i
        
        self.match_table[col+'_eval'] = self.match_table.apply(classification_fucntion, args=(col,), axis=1)

    def add_ptm_log2fx(self):
        """Updates matched table with a column evaluating the log2fx values of the PTM peaks."""

        def ptm_log2fx(row):
            """Function to calculate the log2 fold change of the PTM peaks.
            
            Args:
                row (Series): The row of the matched DataFrame.
            
            Returns:
                log2fx (float): The log2 fold change of the PTM peaks."""

            #get peak id lists
            peak_id_ptm = row['peak_id_ptm']
            peak_id_ptm_a = peak_id_ptm[0]
            peak_id_ptm_b = peak_id_ptm[1]

            #split peak id lists
            peak_id_ptm_a_split = peak_id_ptm_a.split(';')
            peak_id_ptm_b_split = peak_id_ptm_b.split(';')

            #get sum of peak heights (multiple ptms may match to a single peak)
            ptm_a_sum_height = 0
            for i_a in peak_id_ptm_a_split:
                ptm_a_peak_height = peak_table[peak_table['peak_id'] == int(i_a)]['peak_height'].iloc[0]
                ptm_a_sum_height += ptm_a_peak_height

            ptm_b_sum_height = 0
            for i_b in peak_id_ptm_b_split:
                ptm_b_peak_height = peak_table[peak_table['peak_id'] == int(i_b)]['peak_height'].iloc[0]
                ptm_b_sum_height += ptm_b_peak_height

            #calculate log2 fold change
            log2fx = np.log2(ptm_a_sum_height / ptm_b_sum_height)
            return log2fx

        #get tables and filter match table
        peak_table = self.peak_table.copy()
        match_table = self.match_table.copy()
        match_table_filter = match_table[(match_table['global_eval'] == 'Both') & (match_table['phospho_eval'] == 'Both')].copy().reset_index(drop=True)

        #apply ptm_log2fx function
        match_table_filter['ptm_height_log2fx_sum'] = match_table_filter.apply(ptm_log2fx, axis=1)

        #if cluster column values are a list then merge the tables
        if type(match_table['cluster'][0]) == list:
            self.match_table = match_table.merge(match_table_filter[['protein_id', 'mean_peak_apex', 'ptm_height_log2fx_sum']], on=['protein_id', 'mean_peak_apex'], how='left')
        else:
            self.match_table = match_table.merge(match_table_filter[['protein_id', 'cluster', 'mean_peak_apex', 'ptm_height_log2fx_sum']], on=['protein_id', 'cluster', 'mean_peak_apex'], how='left')
    
    def other_edits(self):
        """Updates match/alignment table to move peak parameter columns to the end."""
        self.match_table = self.match_table[self.match_table.columns.drop(list(self.match_table.filter(regex='^peak', axis=1))).tolist() + self.match_table.filter(regex='^peak', axis=1).columns.tolist()]


@click.command()
@click.option('--sql_db', '-s', help='The path to the SQLite database file.', required=True)
@click.option('--peak_table_name', '-p', default='peak_table', help='The name of the table containing the peak table.')
@click.option('--alignment_table_name', '-a', default='alignment_table', help='The name of the table containing the matched table.')
@click.option('--fold_change_threshold', '-t', default=1, help='The fold change threshold.', type=float)
@click.option('--uniprot', '-u', default='sample_input/uniprotk_9606_2023_09_05.tsv', help='The path to the Uniprot data.')
@click.option('--sec_mw', '-m', default='sample_input/SEC_MW_input_hekhct.txt', help='The path to the SEC MW input data.')
@click.option('--first_fraction', '-f', default=8, help='The first fraction number.', type=int)
@click.option('--mass_factor', '-mf', default=1.5, help='The mass factor.', type=float)
@click.option('--output_table_name', '-o', default='alignment_table_modified', help='The name of the output table.')
def main(sql_db, peak_table_name, alignment_table_name, fold_change_threshold, uniprot, sec_mw, first_fraction, mass_factor, output_table_name):
    """Main function to edit the alignment/matched table."""

    #1. Load match/alignment data, uniprot ref, and MW standard and other inputs
    click.echo('Loading data...')
    DB = SQLiteMatchingManager(sql_db)
    DB.load_dataframes(peak_table_name, alignment_table_name)
    match_table= DB.match_table.copy()
    peak_table = DB.peak_table.copy()

    uniprot = pd.read_csv(uniprot, sep='\t')

    mw_standard = pd.read_csv(sec_mw, sep='\t')

    #2. Edit match/alignment data
    click.echo('Editing data:')
    match_table_editor = matchTableEditor(match_table, peak_table)

    #3. Find the mean of the peak apexes of the aligned peaks
    click.echo('Adding mean peak apex...')
    match_table_editor.add_mean_apex()

    #4. Find the fold changes of the peak heights, and evaluate them based on a cutoff.
    click.echo('Adding fold changes...')
    match_table_editor.add_fold_changes()
    match_table_editor.add_fold_change_eval(threshold=fold_change_threshold)

    #5. Add monomer mass and estimated molecular weight - evaluate the region of the peaks apexes.
    click.echo('Adding monomer mass and estimated molecular weight...') 
    match_table_editor.add_monomer_MW(uniprot)
    match_table_editor.add_estimated_MW(mw_standard, first_fraction, plot=True, plot_filename=os.path.dirname(sql_db)+'/MW_Model_Plot.pdf')
    match_table_editor.add_region(mass_factor)
    match_table_editor.classify_col('global')
    match_table_editor.classify_col('phospho')

    #6. Add PTM log2fx by summing the peak heights of the PTM peaks
    click.echo('Adding PTM log2fx...')
    match_table_editor.add_ptm_log2fx()
    match_table_editor.other_edits()

    #7. Save updated alignment table to SQL
    click.echo('Saving data.')
    DB.save_to_SQL(match_table_editor.match_table.copy(), output_table_name)
    DB.close_connection()

if __name__ == '__main__':
    main()
