"""This module contains the peak picking algorithm for the peak picking of the intensity table. 
The peak picking algorithm is based on the intensity table or the smoothed intensity table (can perform smoothing as well)."""

#import libraries
import pandas as pd
import numpy as np
from scipy.signal import find_peaks, peak_widths
import time
import os
import sqlite3
import matplotlib.pyplot as plt
from scipy.signal import filtfilt
from collections import Counter

class SQLitePeakPickingManager:
    """Class to manage the peak picking of the intensity table. The class contains the following functions:
    1. load_dataframe: Load the intensity table from the SQL database.
    2. apply_smoothing: Apply smoothing to the intensity table (optional).
    3. apply_peakpicking: Apply peak picking to the intensity table.
    4. plot_multiplicity_distributions: Plot the multiplicity distributions of the peaks.
    5. save_to_SQL: Save the peak table to the SQL database.
    6. close_connection: Close the connection to the SQL database.
    """

    def __init__(self, db_filename):
        """Initializes the SQLitePeakPickingManager class with the database filename."""
        self.db_filename = db_filename
        self.db_folder_path = os.path.dirname(db_filename)
        self.conn = sqlite3.connect(db_filename)
    
    def load_dataframe(self, table_name):
        """Load the intensity table from the SQL database.
        
        Args:
            table_name (str): The name of the table in the SQL database.
        
        Returns:
            None (Creates the intensity_table and intensity_table_pivot attributes)."""
        self.table_name = table_name
        query = f"SELECT * FROM {table_name};"
        self.intensity_table = pd.read_sql_query(query, self.conn)

        #renames the peptide_id column to cluster and fills in 0 for protein level entries
        if 'peptide_id' in self.intensity_table.columns:
            self.intensity_table = self.intensity_table.rename(columns={'peptide_id':'cluster'})
            self.intensity_table['cluster'] = self.intensity_table['cluster'].fillna('0')

        self.intensity_table_pivot = pd.pivot_table(self.intensity_table.copy(), index=['protein_id', 'gene', 'cluster', 'sample', 'ptm', 'condition', 'replicate'], columns='fraction', values='intensity', fill_value=0).reset_index()

    def apply_smoothing(self, save_to_SQL=True):
        """Apply smoothing to the intensities using a low-pass filter (see reference). See data preprocessing for more details.
        
        Args:
            save_to_SQL (bool): A boolean indicating whether to save the smoothed intensity table to the SQL database.
        
        Returns:
            None (Updates the intensity_table_pivot attribute of the SQLitePeakPickingManager object.)
        
        Reference:
            - https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.filtfilt.html
        """
        #get intensity column names
        intensity_columns = self.intensity_table_pivot.copy()[np.unique(self.intensity_table['fraction'])].columns
        self.intensity_columns = intensity_columns

        #smoothing function
        def smooth_intensities(row_intensities):
            row_intensities = list(row_intensities)

            #the higher the integer smoothing factor, the more intense the smoothing (2 is a good starting point)
            smoothing_factor = 3
            b = [1.0 / smoothing_factor] * smoothing_factor
            a = 1
            y = filtfilt(b, a, row_intensities)
            return pd.Series(y)
    
        #apply smoothing function to each entry
        smoothed_intensity_table = self.intensity_table_pivot.copy()
        smoothed_intensity_table[intensity_columns] = smoothed_intensity_table[intensity_columns].apply(smooth_intensities, axis=1)

        #save to SQL for future use
        if save_to_SQL==True:
            sqlite3.register_adapter(np.int64, lambda val: int(val))
            sqlite3.register_adapter(np.int32, lambda val: int(val))
            smoothed_table_name = self.table_name + '_smoothed'
            smoothed_intensity_table.to_sql(smoothed_table_name, self.conn, index=False, if_exists='replace')

            self.conn.commit()

        self.intensity_table_pivot = smoothed_intensity_table
        return smoothed_intensity_table

    def apply_peakpicking(self):
        """Apply peak picking to the smoothed intensity table. 
        First identifies the peak locations, and then calculates the peak widths and other attributes. 
        This function return a peak table, splitting each peak into a separate row.
        
        The peak picking algorithm uses the following parameters:
        - prominence: The prominence of the peaks (default=0.05 * maximum intensity per profile).
        - distance: The minimum distance between peaks (default=3).
        - height: The minimum height (intensity) of the peaks (default=1000).
        - width: The minimum width of the peaks (default=1).
        
        Returns:
            pd.DataFrame: peak table (also updates the peak_table attribute of the SQLitePeakPickingManager object.)"""

        def peakpicking(row):
            #get intensities
            intensities = np.array(row[self.intensity_columns])
            rel_threshold = np.max(intensities) * 0.05

            #get peaks
            pep_peaks = find_peaks(intensities, prominence=rel_threshold, distance=3, height=1000, width=1)[0]
            result = [[row['protein_id']]*len(pep_peaks), [row['gene']]*len(pep_peaks), [row['cluster']]*len(pep_peaks), [row['sample']]*len(pep_peaks), [row['ptm']]*len(pep_peaks), [row['condition']]*len(pep_peaks), [row['replicate']]*len(pep_peaks), pep_peaks, intensities[pep_peaks]]
            
            #get peak widths
            result.extend(list(peak_widths(intensities, pep_peaks, rel_height=0.3)))

            #edit start/end fractions
            result[12] = result[12] + 2
            result.append(result[11].round())
            result.append(result[12].round())
            return result 
        
        #apply peakpicking algorithm
        intensity_df = self.intensity_table_pivot
        self.intensity_columns = np.unique(self.intensity_table['fraction'])
        df_result = intensity_df.apply(peakpicking, axis=1, result_type='expand')

        #reshape peakpicking results into a peak_table
        df_result.columns = ['protein_id', 'genes', 'cluster', 'sample', 'ptm', 'condition', 'replicate', 'peak_apex', "peak_height","peak_width", "peak_width_height", "peak_left_ips", "peak_right_ips", 'peak_left', 'peak_right']
        df_expanded = pd.DataFrame()
        for col in df_result.columns:
            df_expanded[col] = df_result[col].explode()
        df_expanded = df_expanded.dropna(axis=0, how='all').reset_index(drop=True)
        
        df_expanded.index = df_expanded.index.set_names(['peak_id'])
        df_expanded = df_expanded.reset_index()

        self.peak_table = df_expanded
        return df_expanded
    
    def plot_multiplicity_distributions(self, file_location, file_string_add):
        """Plot the multiplicity distributions of the peaks - the number of peaks per protein per sample.
        
        Args:
            file_location (str): The location to save the multiplicity distributions.
            file_string_add (str): A string to add to the file name.
        
        Output:
            Multiplicity distributions for each sample (saved as a PDF file).
        """
        peak_table = self.peak_table.copy()
        peak_table_select = peak_table[['protein_id', 'sample', 'cluster', 'peak_apex']]
        peak_table_select['peak_apex'] = 1
        peak_table_summary = peak_table_select.groupby(['protein_id', 'sample', 'cluster']).sum().reset_index()

        for s in peak_table_summary['sample'].unique():
            peak_table_summary_s = peak_table_summary[peak_table_summary['sample']==s]
            peakMultiplicity_array = np.array(peak_table_summary_s['peak_apex'])
            mean_value = np.mean(peakMultiplicity_array)
            counter = Counter(peakMultiplicity_array) #peak_table_summary['peak_apex']
            x_values = sorted(counter.keys())
            y_values = [counter[x] for x in x_values]

            # Create bar plot
            plt.figure(figsize=(10, 6))
            plt.bar(x_values, y_values, color='steelblue')
            plt.axvline(x=mean_value, color='red', linestyle='--', label=f'Mean: {mean_value:.2f}')
            plt.xlabel('Number of Peaks')
            plt.ylabel('Protein Frequency')
            plt.title('Peak Multiplicity Histogram: ' + s)
            plt.legend()
            plt.savefig(file_location + "multiplicity_histogram_" + s + "_" + file_string_add + ".pdf", bbox_inches='tight', format='pdf')
            plt.close()
    
    def save_to_SQL(self, peak_table_name='peak_table'):
        """Save the peak table to the SQL database."""
        sqlite3.register_adapter(np.int64, lambda val: int(val))
        sqlite3.register_adapter(np.int32, lambda val: int(val))

        self.peak_table.to_sql(peak_table_name, self.conn, index=False, if_exists='replace')

        self.conn.commit()

    def close_connection(self):
        """Close the connection to the SQL database."""
        self.conn.close()

#start timer - for testing purposes
t0 = time.time()

#1. load SQL file
DB = SQLitePeakPickingManager('20240409_PeakMatching_Package/sample_outputs/DataRepository_Test.db')

#2. load initial intensity table
DB.load_dataframe('initial_intensity_table')

#3. apply smoothing (optional)
DB.apply_smoothing(save_to_SQL=True)

#3. peak picking
DB.apply_peakpicking()

#plot multiplicity distributions
file_string_add = 'multiplicity_distributions_test'
file_location = "20240409_PeakMatching_Package/sample_outputs/"
DB.plot_multiplicity_distributions(file_location, file_string_add)

#4. save results to SQL
DB.save_to_SQL(peak_table_name='peak_table')

#5. close sql connection
DB.close_connection()

#end timer
t1 = time.time()
time_elapsed = t1-t0
print(time_elapsed)

