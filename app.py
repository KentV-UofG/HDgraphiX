#!/usr/bin/env python
# coding: utf-8

# In[3]:


from flask import Flask, flash, render_template, request, send_file, make_response, Response
from io import BytesIO
import logging
import io
import zipfile
import os
import csv
import re
import shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import matplotlib.patheffects as pe
import math
import seaborn as sns
import time
import itertools
import scipy
import tempfile
from pathlib import Path
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import font_manager as fm, rcParams
from matplotlib.colors import to_hex
from scipy.stats import ttest_ind_from_stats
from io import StringIO
from datetime import datetime
from azure.storage.blob import BlobServiceClient, BlobClient, ContainerClient
from azure.core.exceptions import ResourceExistsError

pd.options.display.width = 200

from IPython.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))

app = Flask(__name__, template_folder='templates')
app.config['UPLOAD_FOLDER'] = "Name_Of_An_Azure_Blob_Folder"  # Ensure this folder exists
app.secret_key = 'Key_For_Website_Would_Go_Here'


@app.route('/pdf')
def pdf_view():
    pdf_path = 'images/Logsaso.pdf'
    return send_file(pdf_path, as_attachment=False)

@app.route('/')
def index():
    return render_template('index.html')

# Set up Azure Blob Storage credentials
AZURE_STORAGE_CONNECTION_STRING = "ThisIsObtainedOnAzure"
CONTAINER_NAME = "AContainerForTemporaryOutputStorage" #to ensure file is output properly, deleted almost immediately and never accessed by developers

blob_service_client = BlobServiceClient.from_connection_string(AZURE_STORAGE_CONNECTION_STRING) 
container_client = blob_service_client.get_container_client(CONTAINER_NAME)


output_bitmap_h_count = 0
output_bitmap_v_count = 0
obj_id = ''
seg_id = ''
chain_id = ''
chain_dict = {}    

@app.route('/generate_pdf', methods=['POST'])
def generate_pdf():
    #try:
    request.environ['wsgi.input_terminated'] = 600  # Increase the timeout limit
    timestamp = datetime.now().strftime("%Y.%m.%d_%H.%M.%S")
 
    AZURE_STORAGE_CONNECTION_STRING = "ThisIsObtainedOnAzure"
    CONTAINER_NAME = "AContainerForTemporaryOutputStorage"
    blob_service_client = BlobServiceClient.from_connection_string(AZURE_STORAGE_CONNECTION_STRING) 
    container_client = blob_service_client.get_container_client(CONTAINER_NAME)
    FONT_CONTAINER_NAME = 'Name_Of_An_Azure_Blob_Folder'
    
    fontchoice = float(request.form.get('fontchoiceu', 1))
    if fontchoice == 1:
        FONT_BLOB_NAME = 'Arial.ttf'
        font_path = '/tmp/Arial.ttf'
    elif fontchoice == 2: 
        FONT_BLOB_NAME = 'helvetica.ttf'
        font_path = '/tmp/helvetica.ttf'
    elif fontchoice == 3:
        FONT_BLOB_NAME = 'Times New Roman.ttf'
        font_path = '/tmp/Times New Roman.ttf'
        
    font_container_client = blob_service_client.get_container_client(CONTAINER_NAME)
    blob_client = blob_service_client.get_blob_client(container=FONT_CONTAINER_NAME, blob=FONT_BLOB_NAME)
    font_data = blob_client.download_blob()
    font_bytes = io.BytesIO(font_data.readall())

    with open(font_path, 'wb') as f:
        f.write(font_bytes.getvalue())
    fm.fontManager.addfont(font_path)

    if fontchoice == 1:
        plt.rcParams['font.family'] = 'Arial'
    elif fontchoice == 2:
        plt.rcParams['font.family'] = 'helvetica'
    elif fontchoice == 3:
        plt.rcParams['font.family'] = 'Times New Roman'

    
    ###########################################################################################################################

    
    # Constants and default values. Do not change this section.
    # All of these variables can be set to desired values in the next section.

    output_buffer = []
    heatmap_buffer = []
    buffer = []
    
    # Set up Azure Blob Storage credentials
    AZURE_STORAGE_CONNECTION_STRING = "ThisIsObtainedOnAzure"
    CONTAINER_NAME = "AContainerForTemporaryOutputStorage"
    
    blob_service_client = BlobServiceClient.from_connection_string(AZURE_STORAGE_CONNECTION_STRING) 
    container_client = blob_service_client.get_container_client(CONTAINER_NAME)
    
    # proton mass
    mp = 1.00727647
    
    output_csv_file  = 'HDX processed data.csv'
    output_pdf_file  = 'HDX heatmap.pdf'
    
    mutation_msg = 1

    p_threshold = 0.05
    Dd_threshold = 0.5    
    significant_only = 1
    scatter_plot = 0
    scatter_dir = 'scatter_plots'
    
    split_outp_by_prot = ''
    split_outp_chunks = 0
    split_outp_by_row = []
    
    plot_h = 1
    plot_v = 0

    pymol_print = 1
    pymol_dir = 'pymol_macros'

    def download_blob_as_bytes(container_client, blob_name):
        blob_client = container_client.get_blob_client(blob_name)
        stream = BytesIO()
        blob_client.download_blob().download_to_stream(stream)
        stream.seek(0)
        return stream
    
    def generate_unique_filename(base_name, extension, timestamp):
        unique_filename = f"{base_name}_{timestamp}.{extension}"
        return unique_filename

    def upload_to_blob_storage(container_client, file_content, blob_name):
        blob_client = container_client.get_blob_client(blob_name)
        blob_client.upload_blob(file_content, overwrite=True)
        return blob_name

    
    ###########################################################################################################################
    # --------------------------------------------------------------------------------------------------------
    
    
    #Reset Any Changed Values
    neg_col = None
    pos_col = None
    custom_colors = None
    num_col = None
    num_shades = None
    num_inputs = None
    bound1 = None
    bound2 = None
    bound3 = None
    bound4 = None
    bound5 = None
    bound6 = None
    bound7 = None
    bound8 = None
    custom_bounds = None
    custb = None
    h_or_v =  None
    font_size = 36
    font_size_title = 48
    ao = None
    dtime = None
    drop_times = None
    dpept = None
    drop_pept = None
    dpeppro = None
    dpepst = None
    dpepend = None
    dprot = None
    drop_prot = None
    dpro = None
    renumdict = None
    key1 = None
    value1 = None
    key2 = None
    value2 = None
    renumbering_dict = None
    iddict = None
    mutidWT = None
    mutidMUT = None
    mut_id_dict = None
    mutdict = None
    mutdictres = None
    mutdictwt = None
    mutation_dict = None
    slist = None
    state1_list = None
    state2_list = None
    s1 = None
    s2 = None
    state_list = None
    buffer = None
    input_csv_file = None
    pept_tick_labels = None
    time_tick_labels = None
    f_pymol = None
    zip_file = None
    download_pymol = None
    obj_id = ''
    seg_id = ''
    chain_id = ''
    chain_dict = {}
    output_bitmap_h_count = 0
    output_bitmap_v_count = 0
    output_bitmap_file = None
    relativeUptakeCalc = None
    PulseLabellingNow = None
    separate_plots_pls = 0
    funkybound = None
    globalmax_delta = None
    max_range = None
    altvolc = None
    zerobound = None
    plot_w = None
    color_sections = []
    all_in_min = 0
    what_t_unit = 0
    annotateHM = 0
    addvalHM = False
    editing_figtitle = 0
    output_bitmap = 1
    output_bitmap_name = 'HDX heatmap'
    output_bitmap_dpi = 100
    woodscolthick = 3
    woodsbackthick = 0.5
    wpadthick = 10
    blocktextposwoods = 0.3 #for w_plot
    blockthickwoods = 0.2
    absoluteUptakeValues = 0
    SpecifyStateProt = 0
    PerformAbsoluteUptake = 0
    
    #################################################################################################################################################################
    #################################################################################################################################################################
    #############################################################            Define Functions            ############################################################
    #################################################################################################################################################################
    #################################################################################################################################################################
    
    
    # Define a few functions.
    # from https://stackoverflow.com/questions/24005221/ipython-notebook-early-exit-from-cell
    class StopExecution(Exception):
        def _render_traceback_(self):
            pass
    
    # from https://stackoverflow.com/questions/25668828/how-to-create-colour-gradient-in-python
    def colorFader(c1,c2,mix=0):   #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
        c1 = c1.strip()
        c1=np.array(mc.to_rgb(c1))
        #print('color fader c1:', c1)
        c2 = c2.strip()
        c2=np.array(mc.to_rgb(c2))
        #print('color fader c2:', c2)
        return mc.to_hex((1-mix)*c1 + mix*c2)
    
    def mk_pymol(a,b,c,d):
        chain_ids = chain_dict.get(a, 'default_chain') 
        chain_selection = '+'.join(chain_ids.split(','))
        b = int(b)
        c = int(c)
        return f'alter ///{chain_selection} and resi {b}-{c}, b={d:.3f}'
        
    def mk_chimerax(chain_tocol, start, end, bfactor):
        nstart = int(start)
        nend = int(end)
        command = f"setattr {chain_tocol}:{nstart}-{nend} r colorval {bfactor:.3f} create true"
        return command


    ######THE FOLLOWING FUNCTIONS ARE FOR CONVERTING HDExaminer FILES TO DYNAMX FORMAT#####

    
    def clean_and_divide_entry(entry):
        if entry == 'Full-D':
            divided_entry = 1000
            return divided_entry
        else:
            cleaned_entry = re.sub(r'[a-zA-Z]', '', entry)
            if cleaned_entry:
                try:
                    divided_entry = round(int(cleaned_entry) / 60,3)
                    return divided_entry
                except ValueError:
                    return None
            return None 
    
    def process_first_line(input_file, encoding='utf-8'):
        with open(input_file, 'r', encoding=encoding) as file:
            reader = csv.reader(file)
            first_line = next(reader)
            cleaned_and_divided_entries = [clean_and_divide_entry(entry) for entry in first_line if clean_and_divide_entry(entry) is not None]
            cleaned_and_divided_entries.insert(0, 0)
            #print(cleaned_and_divided_entries)
            num_entries = len(cleaned_and_divided_entries)
            #print(num_entries)
            return num_entries, cleaned_and_divided_entries
    
    def count_rows_after_second(input_file, encoding='utf-8'):
        with open(input_file, 'r', encoding=encoding) as file:
            reader = csv.reader(file)
            next(reader)
            next(reader)
            row_count = sum(1 for row in reader if any(row))
            return row_count
    
    def is_invalid_value(value):
        return value == 0 or value == '' or (isinstance(value, float) and math.isnan(value))
    
    def write_output_csv(output_file, first_row_values, data_rows, repeated_statevalues):
        with open(output_file, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(first_row_values)
            writer.writerow([]) 
            for i in range(len(data_rows)):
                if is_invalid_value(center[i]) and is_invalid_value(RTresult[i]):
                    continue 
                row = [''] * 15 
                if i < len(repeated_statevalues):
                    row[0] = repeated_protlist[i] 
                    row[1] = repeated_startvalues[i]
                    row[2] = repeated_endvalues[i]
                    row[3] = repeated_seqvalues[i]
                    row[6] = repeated_MaxUpvalues[i]
                    row[7] = center[i]
                    row[8] = repeated_statevalues[i]
                    row[9] = repeated_exposurevalues[i]
                    row[10] = repeated_filelist[i]  
                    row[11] = repeated_charge[i]  
                    row[12] = RTresult[i]
                    row[14] = center[i] 
                writer.writerow(row)
                writer.writerow([]) 

    def write_output_csv_summary(output_file, first_row_values, num_rows_summary):
        with open(output_file, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(first_row_values)
            writer.writerow([])
            for i in range(num_rows_summary):
                row = [''] * 15 
                row[0] = repeated_protlist[i]
                row[1] = grab_start_val_sum[i]
                row[2] = grab_end_val_sum[i]
                row[3] = grab_seq_sum[i]
                #rows 4 & 5 (modification and fragment) not necessary
                row[6] = grab_maxup_sum[i]
                row[7] = grab_MHP_sum[i]
                row[8] = grab_state_sum[i]
                row[9] = grab_expo_sum[i]
                row[10] = repeated_filelist[i] #this returns a repeated gen_file list for the file column (likely not necessary)
                row[11] = repeated_charge[i] #Because HDExaminer returns processed D uptake, want Als script to "ignore" charge, treat all as 1
                row[12] = grab_RT_sum[i]
                #row 13 is inten in DynamX files. There is no direct comparison in the HDExaminer file format except maybe score, however the scale is vastly different
                row[14] = grab_center_val_sum[i]
                writer.writerow(row)
                writer.writerow([])

    def grab_columns_summary(csv_file_path, column_index_summary):
        values = []
        with open(csv_file_path, 'r') as file:
            reader = csv.reader(file)
            next(reader, None)
            for row in reader:
                if len(row) > column_index_summary:  
                    value = row[column_index_summary].strip() 
                    check_all_val = [row[i].strip() for i in [0, 2, 3, 4, 5, 6, 7, 8, 10] if len(row) > i] 
                    if len(check_all_val) == 9 and all(check_all_val):
                        if value:  
                            if column_index_summary == grab_expo_sum_col:
                                value = float(value) / 60
                                values.append(value)
                            else:
                                values.append(value)
        return values
    
    def gather_and_condense_columns(input_file, num_entries=None):
        with open(input_file, 'r') as f:
            f.seek(0)
            df = pd.read_csv(f)
        first_non_empty_row_idx = None
        search_value = 'Start RT'
        search_value_2 = 'Search RT'
        for idx in range(len(df)):
            if pd.notna(df.iloc[idx, 0]) and df.iloc[idx, 0] != '.':
                first_non_empty_row_idx = idx
                break
        if first_non_empty_row_idx is None:
            #print("No non-empty cell found in the first column.")
            return []
        for idx in range(len(df)):
            if any(df.iloc[idx].astype(str).str.contains(search_value, na=False)):
                row_with_d_idx = idx
                break
        last_non_empty_row_idx = None
        for idx in range(len(df) - 1, first_non_empty_row_idx - 1, -1):
            if pd.notna(df.iloc[idx, 0]):
                last_non_empty_row_idx = idx
                break
        row_to_search = df.iloc[row_with_d_idx]
        column_indices = [idx for idx, value in enumerate(row_to_search) if value == search_value or value == search_value_2]
        #print(f"First non-empty row index: {first_non_empty_row_idx}")
        #print(f"Last non-empty row index: {last_non_empty_row_idx}")
        #print(f"Column indices with '{search_value}' or '{search_value_2}' in row {first_non_empty_row_idx + 1}: {column_indices}")
        gathered_values = []
        for idx in range(first_non_empty_row_idx + 1, last_non_empty_row_idx + 1):
            if df.iloc[idx].isna().all():
                continue
            row_values = [df.iloc[idx, col_idx] if pd.notna(df.iloc[idx, col_idx]) else 0 for col_idx in column_indices]
            gathered_values.extend(row_values)
        return gathered_values
    
    def replace_nan_with_average(df, idx, col_idx, column_indices):
        exposure_time = df.iloc[2, col_idx] 
        similar_values = []
        for compare_col_idx in column_indices:
            compare_exposure_time = df.iloc[2, compare_col_idx]
            compare_d_value = df.iloc[idx, compare_col_idx]
            if compare_exposure_time == exposure_time and pd.notna(compare_d_value):
                similar_values.append(compare_d_value)
        if similar_values:
            return sum(similar_values) / len(similar_values)
        else:
            return 0 
    
    def gather_D(input_file, num_entries=None):
        df = pd.read_csv(input_file)
        row_with_d_idx = None
        search_value = '#D'
        for idx in range(len(df)):
            if any(df.iloc[idx].astype(str).str.contains(search_value, na=False)): 
                row_with_d_idx = idx
                break
        if row_with_d_idx is None:
            #print(f"No row contains '{search_value}' in the DataFrame.")
            return []
        #print(f"Row containing '{search_value}': {row_with_d_idx + 1}")
        row_to_search = df.iloc[row_with_d_idx]
        column_indices = [idx for idx, value in enumerate(row_to_search) if value == search_value]
        if not column_indices:
            #print(f"No columns found with '{search_value}' in row {row_with_d_idx + 1}.")
            return []
        #print(f"Column indices with '{search_value}' in row {row_with_d_idx + 1}: {column_indices}")
        last_non_empty_row_idx = None
        for idx in range(len(df) - 1, row_with_d_idx, -1):
            if pd.notna(df.iloc[idx, 0]):  
                last_non_empty_row_idx = idx
                break
        #print(f"Last non-empty row index: {last_non_empty_row_idx + 1}")
        gathered_values = []
        for idx in range(row_with_d_idx + 1, last_non_empty_row_idx + 1):
            if df.iloc[idx].isna().all():
                continue
            row_values = [0]
            for col_idx in column_indices:
                d_value = df.iloc[idx, col_idx]
                if pd.isna(d_value) or d_value == '':
                    d_value = replace_nan_with_average(df, idx, col_idx, column_indices)
                row_values.append(d_value)
            gathered_values.extend(row_values) 
        return gathered_values
    
    def create_repeated_string_list(repeated_statevalues, string_to_repeat):
        return [string_to_repeat] * len(repeated_statevalues)

    def create_repeated_string_list_sum(num_rows_summary, string_to_repeat):
        return [string_to_repeat] * num_rows_summary
    
    def count_empty_cells_in_column_a(input_file):
        empty_cell_count = 1
        with open(input_file, 'r') as file:
            reader = csv.reader(file)
            for row in reader:
                if len(row) > 0 and row[0] == '':
                    empty_cell_count += 1
                elif len(row) > 0 and row[0] != '':
                    break 
        return empty_cell_count
    
    def count_empty_cells_between_first_two_filled_cells_first_row(input_file):
        empty_cells_count = 0
        found_first_filled = False
        with open(input_file, 'r') as file:
            csv_reader = csv.reader(file)
            for row in csv_reader:
                for cell_value in row:
                    if cell_value != '':
                        if not found_first_filled:
                            found_first_filled = True
                        else:
                            return empty_cells_count
                    elif found_first_filled:
                        empty_cells_count += 1
        return empty_cells_count

    def count_empty_cells_from_csv(input_file):
        col_empty_cell_count = 0
        with open(input_file, mode='r') as file:
            reader = csv.reader(file)
            data = list(reader)
        first_filled_index = None
        second_filled_index = None
        for row_idx, row in enumerate(data):
            if len(row) > 0 and row[0] != "":
                if first_filled_index is None:
                    first_filled_index = row_idx
                elif second_filled_index is None:
                    second_filled_index = row_idx
                    break 
        for i in range(first_filled_index, second_filled_index):
            if len(data[i]) > 0 and data[i][0] == "":
                col_empty_cell_count += 1
        return col_empty_cell_count

    def gather_column_values(input_file_path, search_term, num_entries):
        collected_values = []
        with open(input_file_path, 'r') as file:
            reader = csv.reader(file)
            first_filled_row = None
            for row in reader:
                if row and row[0].strip(): 
                    first_filled_row = row
                    break
            if first_filled_row is None:
                raise ValueError("No non-empty cell found in the first column.")
            try:
                column_index = first_filled_row.index(search_term)
            except ValueError:
                raise ValueError(f"'{search_term}' not found in the first non-empty row.")
            for row in reader:
                if row: 
                    value = row[column_index].strip()
                    if value: 
                        collected_values.append(value)
            collected_values = list(filter(None, collected_values))
        repeated_values = []
        for value in collected_values:
            repeated_values.extend([value] * num_entries)
        return repeated_values

    def write_output_csv_WB(output_file, df_filtered, column_mapping):
        data_dict = {output_col: df_filtered[original_col].tolist() if original_col in df_filtered.columns else [''] * len(df_filtered)
                     for output_col, original_col in column_mapping.items()}
        with open(output_file, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(column_mapping.keys()) 
            writer.writerow([])  
            num_rows = len(df_filtered)
            for i in range(num_rows):
                row = [data_dict[col][i] for col in column_mapping.keys()]
                writer.writerow(row)
                writer.writerow([]) 
        print(f"CSV file '{output_file}' written successfully.")
    
    def ensure_hex_color(color):
        return color if color.startswith('#') else f'#{color}'

    def save_to_azure_blob(data_avg, container_name, blob_name, connection_string):
        try:
            csv_buffer = StringIO()
            data_avg.to_csv(csv_buffer, index=False)
            csv_data = csv_buffer.getvalue()
            blob_service_client = BlobServiceClient.from_connection_string(connection_string)
            container_client = blob_service_client.get_container_client(container_name)
            blob_client = container_client.get_blob_client(blob_name)
            blob_client.upload_blob(csv_data, overwrite=True)
            print(f"File {blob_name} successfully uploaded to container {container_name}.")
        except Exception as e:
            print(f"An error occurred: {e}")
    
    
    #################################################################################################################################################################
    #################################################################################################################################################################
    #############################################################            Read Input Data            #############################################################
    #################################################################################################################################################################
    #################################################################################################################################################################

    
    PDFgeneration = float(request.form.get('gen_pdf',0))
    file_typeDoH = float(request.form['file_type'])
    if file_typeDoH == 1: # HDeXaminer Input (either Pool or Summary, automatically determined)
        input_file = request.files['csv_file']
        upload_folder = os.path.join('/home', 'site', 'wwwroot', 'heatmap1vahidi')
        os.makedirs(upload_folder, exist_ok=True)
        input_file_path = os.path.join(upload_folder, input_file.filename)
        input_file.save(input_file_path) 
        input_csv_file = 'output_file.csv'
        with open(input_file_path, 'r') as file:
            reader = csv.reader(file)
            first_rowread = next(reader, None) 
            if first_rowread is not None and len(first_rowread) > 0:
                first_cellread = first_rowread[0]
                if first_cellread == '' or first_cellread == '.':  # Check if the first cell is blank
                    resultread = 1  
                else:
                    resultread = 2 
            else:
                resultread = 0 
        if resultread == 1: # POOL
            try:
                num_entries, cleaned_and_divided_entries = process_first_line(input_file_path, encoding='utf-8')
                rows_after_second = count_rows_after_second(input_file_path, encoding='utf-8')
            except UnicodeDecodeError:
                print("UnicodeDecodeError encountered with 'utf-8' encoding, trying 'iso-8859-1' encoding")
                num_entries, cleaned_and_divided_entries = process_first_line(input_file_path, encoding='iso-8859-1')
                rows_after_second = count_rows_after_second(input_file_path, encoding='iso-8859-1')
            empty_cells_count = count_empty_cells_in_column_a(input_file_path)
            #print(f"Number of empty cells in column A: {empty_cells_count}")
            result_empty_cells_count = count_empty_cells_between_first_two_filled_cells_first_row(input_file_path)
            #print(f"Number of empty cells between the first two filled cells in the first row: {result_empty_cells_count}")
            col_empty_cell_count = count_empty_cells_from_csv(input_file_path)
            #print(f"Number of empty cells between the first two filled cells in the first column: {col_empty_cell_count}")
            # Gather and repeat values for different search terms
            repeated_startvalues = gather_column_values(input_file_path, 'Start', num_entries)
            #print('# Start:', len(repeated_startvalues))
            repeated_statevalues = gather_column_values(input_file_path, 'State', num_entries)
            #print('# State:', len(repeated_statevalues))
            repeated_endvalues = gather_column_values(input_file_path, 'End', num_entries)
            #print('# End:', len(repeated_endvalues))
            repeated_seqvalues = gather_column_values(input_file_path, 'Sequence', num_entries)
            #print('# Sequence:', len(repeated_seqvalues))
            repeated_MaxUpvalues = gather_column_values(input_file_path, 'Max D', num_entries)
            #print('# MaxUp:', len(repeated_MaxUpvalues))
            repeated_chargevalues = gather_column_values(input_file_path, 'Charge', num_entries)
            #print('# Charge:', len(repeated_chargevalues))
            repeated_protlist = gather_column_values(input_file_path, 'Protein', num_entries)
            if len(repeated_protlist) != len(repeated_chargevalues):
                prot_empty_pool = True
                string_to_repeat = 'gen_prot'
                repeated_protlist = create_repeated_string_list(repeated_statevalues, string_to_repeat)
            else:
                prot_empty_pool = False  
            #print('# Protein:', len(repeated_protlist))
            input_statevalues = []
            with open(input_file_path, 'r') as file:
                reader = csv.reader(file)
                first_filled_row = None
                for i, row in enumerate(reader):
                    if row and row[0]:
                        first_filled_row = row
                        break
                if first_filled_row is None:
                    print("No filled cells found in column one.")
                else:
                    try:
                        column_index = first_filled_row.index('State')
                    except ValueError:
                        print("The term 'State' was not found in the first filled row.")
                        column_index = -1
                    if column_index != -1: 
                        file.seek(0)  
                        for row in reader:
                            if row and len(row) > column_index:
                                input_statevalues.append(row[column_index])
            input_statevalues = list(filter(None, input_statevalues))
            repeated_exposurevalues = []
            for i in range(len(input_statevalues)):
                repeated_exposurevalues.extend(cleaned_and_divided_entries)
            #print('# Exposure:', len(repeated_exposurevalues))
            first_row_values = ['Protein', 'Start', 'End', 'Sequence', 'Modification',
                               'Fragment', 'MaxUptake', 'MHP', 'State', 'Exposure',
                               'File', 'z', 'RT', 'Inten', 'Center']
            data_rows = [[''] * 15 for _ in range(len(repeated_statevalues))]
            RTresult = gather_and_condense_columns(input_file_path, num_entries)
            #print(RTresult)
            #print('# RTresult:', len(RTresult))
            center = gather_D(input_file_path, num_entries)
            #print(center)
            #print('# Center:', len(center))
            string_to_repeat = 'gen_file'
            repeated_filelist = create_repeated_string_list(repeated_statevalues, string_to_repeat)
            string_to_repeat = '1'
            repeated_charge = create_repeated_string_list(repeated_statevalues, string_to_repeat)
            #write_output_csv(output_file, first_row_values, data_rows, repeated_statevalues)
            #print(f"Number of empty cells in column A: {empty_cells_count}")
            #print(f'Number of exposures: {num_entries}')
            #print(f'Exposure Times (converted to min): {cleaned_and_divided_entries}')
            #print(f'Number of rows after the second row (peptides): {rows_after_second}')
            csv_filename = "output_file.csv"
            csv_path = os.path.join(os.getcwd(), csv_filename)
            container_name = "heatmap1vahidi"
            blob_client = blob_service_client.get_blob_client(container=container_name, blob=csv_filename)
            write_output_csv(csv_path, first_row_values, data_rows, repeated_statevalues)
            output_blob_name = generate_unique_filename('output_file', 'csv', timestamp)
            input_csv_file = csv_path
        elif resultread == 2: # SUMMARY
            csv_filename = "output_file.csv"
            csv_path = os.path.join(os.getcwd(), csv_filename)
            container_name = "heatmap1vahidi"
            blob_client = blob_service_client.get_blob_client(container=container_name, blob=csv_filename)
            first_row_values = ['Protein', 'Start', 'End', 'Sequence', 'Modification',
                               'Fragment', 'MaxUptake', 'MHP', 'State', 'Exposure',
                               'File', 'z', 'RT', 'Inten', 'Center']
            grab_expo_sum_col = 7
            grab_start_val_sum_col = 2 # 0 based indexing, therefore 3rd column = 2, grabbing C column "start"
            grab_start_val_sum = grab_columns_summary(input_file_path, grab_start_val_sum_col)
            #print("start val list:")
            #print(grab_start_val_sum)
            num_rows_summary = len(grab_start_val_sum)
            #print("length of start val list:")
            #print(num_rows_summary)
            grab_prot_sum_col = 1
            repeated_protlist = grab_columns_summary(input_file_path, grab_prot_sum_col)
            #print(len(repeated_protlist))
            if len(repeated_protlist) != num_rows_summary:
                string_to_repeat = 'gen_prot'
                repeated_protlist = create_repeated_string_list_sum(num_rows_summary, string_to_repeat)
            grab_end_val_sum_col = 3
            grab_end_val_sum = grab_columns_summary(input_file_path, grab_end_val_sum_col)
            #print(len(grab_end_val_sum))
            grab_seq_sum_col = 4 #grabbing sequences
            grab_seq_sum = grab_columns_summary(input_file_path, grab_seq_sum_col)
            #print(len(grab_seq_sum))
            #Dont Need Mut or Frag Lists
            grab_maxup_sum_col = 8
            grab_maxup_sum = grab_columns_summary(input_file_path, grab_maxup_sum_col)
            #print(len(grab_maxup_sum))
            grab_MHP_sum_col = 5
            grab_MHP_sum =  grab_columns_summary(input_file_path, grab_MHP_sum_col)
            #print(len(grab_MHP_sum))
            grab_state_sum_col = 0
            grab_state_sum = grab_columns_summary(input_file_path, grab_state_sum_col)
            #print(len(grab_state_sum))
            grab_expo_sum = grab_columns_summary(input_file_path, grab_expo_sum_col)
            #print(len(grab_expo_sum))
            string_to_repeat = 'gen_file' #dont need a real file name, therefore gen_file is used to file in file to be read properly
            repeated_filelist = create_repeated_string_list_sum(num_rows_summary, string_to_repeat)
            #print(len(repeated_filelist))
            string_to_repeat = '1'  #HDExaminer outputs processed D vals, therefore want code to "ignore" so treat as 1 charge state
            repeated_charge = create_repeated_string_list_sum(num_rows_summary, string_to_repeat)
            #print(len(repeated_charge))
            grab_RT_sum_col = 6
            grab_RT_sum = grab_columns_summary(input_file_path, grab_RT_sum_col)
            #print(len(grab_RT_sum))
            #Dont need inten, no direct value in HDExaminer files
            grab_center_val_sum_col = 10
            grab_center_val_sum = grab_columns_summary(input_file_path, grab_center_val_sum_col)
            #print(len(grab_center_val_sum))
            write_output_csv_summary(csv_path, first_row_values, num_rows_summary)
            output_blob_name = generate_unique_filename('output_file', 'csv', timestamp)
            input_csv_file = csv_path
        else:
            return "Record not found", 400
    elif file_typeDoH == 2: # HDX workbench
        input_file = request.files['csv_file']
        upload_folder = os.path.join('/home', 'site', 'wwwroot', 'heatmap1vahidi') 
        os.makedirs(upload_folder, exist_ok=True)
        input_file_path = os.path.join(upload_folder, input_file.filename)
        input_file.save(input_file_path)
        #input_csv_file = 'output_file.csv'
        csv_filename = "output_file.csv"
        csv_path = os.path.join(os.getcwd(), csv_filename)
        container_name = "heatmap1vahidi"
        blob_client = blob_service_client.get_blob_client(container=container_name, blob=csv_filename)
        with open(input_file_path, 'r') as f:
            lines = f.readlines()
        start_row = next(i for i, line in enumerate(lines) if line.strip().startswith("peptide"))
        #print('start_row', start_row)
        df = pd.read_csv(input_file_path, skiprows=start_row-1, header=0)
        # If there are too many unnamed headers, adjust to the next row
        if df.columns.str.contains("Unnamed").sum() > len(df.columns) // 2:
            df = pd.read_csv(input_file_path, skiprows=start_row, header=0)
            if df.columns.str.contains("Unnamed").sum() > len(df.columns) // 2:
                df = pd.read_csv(input_file_path, skiprows=start_row+1, header=0)
        #print('df', df)
        columns_to_collect = ['project', 'start', 'end', 'peptide', 'numExHydrogens', 'monoisotopic', 'sample', 'timepoint', 'filename', 'charge', 'rt_start_replicate', 'centroid', 'discarded_replicate']
        df_filtered = df[columns_to_collect]
        #print(f"Filtered DataFrame shape: {df_filtered.shape}")
        #print(df_filtered['discarded_replicate'].dtype)
        #print(df_filtered['discarded_replicate'].unique())
        df_filtered.loc[:, 'discarded_replicate'] = df_filtered['discarded_replicate'].astype(str).str.strip()
        #print("Unique values after strip:")
        #print(df_filtered['discarded_replicate'].unique())
        #df_filtered['discarded_replicate'] = df_filtered['discarded_replicate'].replace({'NA': np.nan, 'TRUE': True, 'FALSE': False})
        #print("Unique values after replacing 'NA', 'TRUE', 'FALSE' with NaN/booleans:")
        #print(df_filtered['discarded_replicate'].unique())
        df_filtered = df_filtered[df_filtered['discarded_replicate'].notna()]
        #print(f"Filtered DataFrame shape after removing NaN: {df_filtered.shape}")
        df_filtered = df_filtered[df_filtered['discarded_replicate'] == 'False']
        #print(f"Filtered DataFrame shape after filtering for False: {df_filtered.shape}")
        df_filtered['timepoint'] = df_filtered['timepoint'].astype(str).str.replace(r'[a-zA-Z]', '', regex=True)
        #print('Filtered Timepoints to drop s')
        column_mapping = {
            'Protein': 'project',
            'Start': 'start',
            'End': 'end',
            'Sequence': 'peptide',
            'Modification': '',
            'Fragment': '',
            'MaxUptake': 'numExHydrogens',
            'MHP': 'monoisotopic',
            'State': 'sample',
            'Exposure': 'timepoint',
            'File': 'filename',
            'z': 'charge',
            'RT': 'rt_start_replicate',
            'Inten': '',
            'Center': 'centroid',
        }
        write_output_csv_WB(csv_path, df_filtered, column_mapping)
        #container_name = "heatmap1vahidi"
        #blob_client = blob_service_client.get_blob_client(container=container_name, blob=csv_filename)
        #with open(csv_path, "rb") as data:
        #    blob_client.upload_blob(data, overwrite=True)
        output_blob_name = generate_unique_filename('output_file', 'csv', timestamp)
        input_csv_file = csv_path
    else: # DynamX - Read as normal
        input_csv_file = request.files['csv_file']

    
    #################################################################################################################################################################
    #################################################################################################################################################################
    ########################################################            Create & Gather Variables            ########################################################
    #################################################################################################################################################################
    #################################################################################################################################################################
    
    
    #Col List
    redc = '#c50f15'
    dredc = '#990000'
    bluec = '#0062cc'
    dbluec = '#08306b'
    whitec = '#FFFFFF'
    blackc = '#000000'
    greyc = '#E0E0E0'
    greenc = '#006600'
    yellowc = '#e69b00'
    purplec = '#643d6e'

    c1 = bluec    #blue 
    c2 = whitec   #white
    c3 = redc     #red
    c_missing = '#bdbdbd' # gray;

    separate_plots_pls = 0
    new_x_axis_title = None

    h_or_v = float(request.form['h_or_v'])
    if h_or_v == 1:
        plot_v = 0
        plot_h = 1
        plot_w = 0
        plot_volc = 0
        plot_l = 0
    elif h_or_v == 2:
        plot_v = 1
        plot_h = 0
        plot_w = 0
        plot_volc = 0
        plot_l = 0
    elif h_or_v == 3:
        plot_v = 0
        plot_h = 0
        plot_w = 1
        plot_volc = 0
        plot_l = 0
    elif h_or_v == 4:
        plot_v = 0
        plot_h = 0
        plot_w = 0
        plot_volc = 1
        plot_l = 0
    elif h_or_v == 6:
        plot_v = 0
        plot_h = 0
        plot_w = 0
        plot_volc = 0
        plot_l = 1
    if plot_volc == 1:
        scatter_plot = 1
    elif plot_volc == 0:
        scatter_plot = 0
    
    expType = float(request.form.get('expType', 1))
    PerformAbsoluteUptake = float(request.form.get('AbsoluteUp',0))
    if PerformAbsoluteUptake == 1:
        expType = 3
        absoluteUptakeValues = 1 #don't want to remake code for RFU type plots, therefore set to RFU and add another variable I can use to turn off unnecessary processing
    if expType == 2:
        PulseLabellingNow = 2
        #print('PulseLabellingNow Value:', PulseLabellingNow)
    # 1 = continuous labelling 2 = pulse labelling 3 = RFU
    
    uploaded_settings = int(request.form.get('use_uploaded_settings', 0))
    if uploaded_settings == 1:
        settings_file = request.files['csv_settings']
        file_content = settings_file.read().decode("utf-8").splitlines()
        variables = {}
        csv_reader = csv.DictReader(file_content)
        #print("Local variables:", locals())  # Prints local scope
        #print("Global variables:", globals())
        for row in csv_reader:
            var_name = row.get("Variable Name", "").strip()
            data_type = row.get("Data Type", "").strip().lower()
            var_value = row.get("Value", "").strip()
            #print(var_name, var_value)
            try:
                if data_type == "int":
                    var_value = int(var_value)
                elif data_type == "float":
                    var_value = float(var_value)
                elif data_type == "bool":
                    var_value = var_value.lower() == "true"
                elif data_type == "list":
                    var_value = eval(var_value) if var_value.startswith("[") else var_value.split(",")
                elif data_type == "dict":
                    var_value = eval(var_value) if var_value.startswith("{") else {}
                elif data_type == "list":
                    # Safely evaluate the list
                    var_value = eval(var_value) if var_value.startswith("[") else var_value.split(",")
                elif data_type == "dict":
                    # Safely evaluate the dictionary
                    var_value = eval(var_value) if var_value.startswith("{") else {}
                else:
                    var_value = str(var_value)
                globals()[var_name] = var_value
                variables[var_name] = var_value
                #print('Variables:', variables)
            except Exception as e:
                print(f"Error processing variable {var_name}: {e}")
        if 'fontchoice' in variables:
            try:
                fontchoice = variables.get('fontchoice')
            except ValueError:
                #print('fontchoice not added/updated')
                pass
            #else:
                #print('fontchoice = ', fontchoice)
        if 'PDFgeneration' in variables:
            try:
                PDFgeneration = variables.get('PDFgeneration')
            except ValueError:
                #print('PDFgeneration not added/updated')
                pass
            #else:
                #print('PDFgeneration = ', PDFgeneration)
        #if 'file_typeDoH' in variables:
        #    try:
        #        file_typeDoH = variables.get('file_typeDoH')
        #    except ValueError:
                #print('file_typeDoH not added/updated')
        #        pass
            #else:
                #print('file_typeDoH = ', file_typeDoH)
        #if 'input_file' in variables:
        #    try:
        #        input_file = variables.get('input_file')
        #    except ValueError:
        #        print('input_file not added/updated')
        #        pass
        #    else:
        #        print('input_file = ', input_file)
        #if 'input_csv_file' in variables:
        #    try:
        #        input_csv_file = variables.get('input_csv_file')
        #    except ValueError:
        #        print('input_csv_file not added/updated')
        #        pass
        #    else:
        #        print('input_csv_file = ', input_csv_file)
        if 'alt_c2' in variables:
            try:
                alt_c2 = variables.get('alt_c2')
            except ValueError:
                #print('alt_c2 not added/updated')
                pass
            #else:
                #print('alt_c2 = ', alt_c2)
        if 'c2' in variables:
            try:
                c2 = variables.get('c2')
            except ValueError:
                #print('c2 not added/updated')
                pass
            #else:
                #print('c2 = ', c2)
        if 'alt_cmissing' in variables:
            try:
                alt_cmissing = variables.get('alt_cmissing')
            except ValueError:
                #print('alt_cmissing not added/updated')
                pass
            #else:
                #print('alt_cmissing = ', alt_cmissing)
        if 'c_missing' in variables:
            try:
                c_missing = variables.get('c_missing')
            except ValueError:
                #print('c_missing not added/updated')
                pass
            #else:
                #print('c_missing = ', c_missing)
        if 'color_by_heatmap' in variables:
            try:
                color_by_heatmap = variables.get('color_by_heatmap')
            except ValueError:
                #print('color_by_heatmap not added/updated')
                pass
            #else:
                #print('color_by_heatmap = ', color_by_heatmap)
        if 'colcutopt' in variables:
            try:
                colcutopt = variables.get('colcutopt')
            except ValueError:
                #print('colcutopt not added/updated')
                pass
            #else:
                #print('colcutopt = ', colcutopt)
        if 'colcutoff' in variables:
            try:
                colcutoff = variables.get('colcutoff')
            except ValueError:
                #print('colcutoff not added/updated')
                pass
            #else:
                #print('colcutoff = ', colcutoff)
        if 'woodsdimen' in variables:
            try:
                woodsdimen = variables.get('woodsdimen')
            except ValueError:
                #print('woodsdimen not added/updated')
                pass
            #else:
                #print('woodsdimen = ', woodsdimen)
        if 'woodsx' in variables:
            try:
                woodsx = variables.get('woodsx')
            except ValueError:
                #print('woodsx not added/updated')
                pass
            #else:
                #print('woodsx = ', woodsx)
        if 'woodsy' in variables:
            try:
                woodsy = variables.get('woodsy')
            except ValueError:
                #print('woodsy not added/updated')
                pass
            #else:
                #print('woodsy = ', woodsy)
        if 'woodscol' in variables:
            try:
                woodscol = variables.get('woodscol')
            except ValueError:
                #print('woodscol not added/updated')
                pass
            #else:
                #print('woodscol = ', woodscol)
        if 'woodscolpos' in variables:
            try:
                woodscolpos = variables.get('woodscolpos')
            except ValueError:
                #print('woodscolpos not added/updated')
                pass
            #else:
                #print('woodscolpos = ', woodscolpos)
        if 'woodscolneu' in variables:
            try:
                woodscolneu = variables.get('woodscolneu')
            except ValueError:
                #print('woodscolneu not added/updated')
                pass
            #else:
                #print('woodscolneu = ', woodscolneu)
        if 'woodscolneg' in variables:
            try:
                woodscolneg = variables.get('woodscolneg')
            except ValueError:
                #print('woodscolneg not added/updated')
                pass
            #else:
                #print('woodscolneg = ', woodscolneg)
        if 'nolines' in variables:
            try:
                nolines = variables.get('nolines')
            except ValueError:
                #print('nolines not added/updated')
                pass
            #else:
                #print('nolines = ', nolines)
        if 'font_size' in variables:
            try:
                font_size = variables.get('font_size')
            except ValueError:
                #print('font_size not added/updated')
                pass
            #else:
                #print('font_size = ', font_size)
        if 'font_size_title' in variables:
            try:
                font_size_title = variables.get('font_size_title')
            except ValueError:
                #print('font_size_title not added/updated')
                pass
            #else:
                #print('font_size_title = ', font_size_title)
        if 'font_size_ticklabel' in variables:
            try:
                font_size_ticklabel = variables.get('font_size_ticklabel')
            except ValueError:
                #print('font_size_ticklabel not added/updated')
                pass
            #else:
                #print('font_size_ticklabel = ', font_size_ticklabel)
        if 'fontcolorc' in variables:
            try:
                fontcolorc = variables.get('fontcolorc')
            except ValueError:
                #print('fontcolorc not added/updated')
                pass
            #else:
                #print('fontcolorc = ', fontcolorc)
        #if 'h_or_v' in variables:
        #    try:
        #        h_or_v = variables.get('h_or_v')
        #    except ValueError:
                #print('h_or_v not added/updated')
        #        pass
            #else:
                #print('h_or_v = ', h_or_v)
        if 'zerobound' in variables:
            try:
                zerobound = variables.get('zerobound')
            except ValueError:
                #print('zerobound not added/updated')
                pass
            #else:
                #print('zerobound = ', zerobound)
        if 'funkybound' in variables:
            try:
                funkybound = variables.get('funkybound')
            except ValueError:
                #print('funkybound not added/updated')
                pass
            #else:
                #print('funkybound = ', funkybound)
        if 'globalmax_delta' in variables:
            try:
                globalmax_delta = variables.get('globalmax_delta')
            except ValueError:
                #print('globalmax_delta not added/updated')
                pass
            #else:
                #print('globalmax_delta = ', globalmax_delta)
        if 'p_threshold' in variables:
            try:
                p_threshold = variables.get('p_threshold')
            except ValueError:
                #print('p_threshold not added/updated')
                pass
            #else:
                #print('p_threshold = ', p_threshold)
        if 'custb' in variables:
            try:
                custb = variables.get('custb')
            except ValueError:
                print('custb not added/updated')
                pass
            else:
                print('custb = ', custb)
        if 'max_range' in variables:
            try:
                max_range = variables.get('max_range')
            except ValueError:
                #print('max_range not added/updated')
                pass
            #else:
                #print('max_range = ', max_range)
        if 'num_shades' in variables:
            try:
                num_shades = variables.get('num_shades')
            except ValueError:
                print('num_shades not added/updated')
                pass
            else:
                print('num_shades = ', num_shades)
        else:
            print('Cant find num_shades')
        if 'alt_col' in variables:
            try:
                alt_col = variables.get('alt_col')
            except ValueError:
                #print('alt_col not added/updated')
                pass
            #else:
                #print('alt_col = ', alt_col)
        if 'c1' in variables:
            try:
                c1 = variables.get('c1')
            except ValueError:
                #print('c1 not added/updated')
                pass
            #else:
                #print('c1 = ', c1)
        if 'c3' in variables:
            try:
                c3 = variables.get('c3')
            except ValueError:
                #print('c3 not added/updated')
                pass
            #else:
                #print('c3 = ', c3)
        if 'nb' in variables:
            try:
                nb = variables.get('nb')
            except ValueError:
                #print('nb not added/updated')
                pass
            #else:
                #print('nb = ', nb)
        if 'pb' in variables:
            try:
                pb = variables.get('pb')
            except ValueError:
                #print('pb not added/updated')
                pass
            #else:
                #print('pb = ', pb)
        if 'neg_col' in variables:
            try:
                neg_col = variables.get('neg_col')
            except ValueError:
                #print('neg_col not added/updated')
                pass
            #else:
                #print('neg_col = ', neg_col)
        if 'pos_col' in variables:
            try:
                pos_col = variables.get('pos_col')
            except ValueError:
                #print('pos_col not added/updated')
                pass
            #else:
                #print('pos_col = ', pos_col)
        if 'nb' in variables:
            try:
                nb = variables.get('nb')
            except ValueError:
                #print('nb not added/updated')
                pass
            #else:
                #print('nb = ', nb)
        if 'pb' in variables:
            try:
                pb = variables.get('pb')
            except ValueError:
                #print('pb not added/updated')
                pass
            #else:
                #print('pb = ', pb)
        if 'colorchoicetype' in variables:
            try:
                colorchoicetype = variables.get('colorchoicetype')
            except ValueError:
                #print('colorchoicetype not added/updated')
                pass
            #else:
                #print('colorchoicetype = ', colorchoicetype)
        if 'negshadein' in variables:
            try:
                negshadein = variables.get('negshadein')
            except ValueError:
                #print('negshadein not added/updated')
                pass
            #else:
                #print('negshadein = ', negshadein)
        if 'posshadein' in variables:
            try:
                posshadein = variables.get('posshadein')
            except ValueError:
                #print('posshadein not added/updated')
                pass
            #else:
                #print('posshadein = ', posshadein)
        if 'nc' in variables:
            try:
                nc = variables.get('nc')
            except ValueError:
                #print('nc not added/updated')
                pass
            #else:
                #print('nc = ', nc)
        if 'pc' in variables:
            try:
                pc = variables.get('pc')
            except ValueError:
                #print('pc not added/updated')
                pass
            #else:
                #print('pc = ', pc)
        if 'dtime' in variables:
            try:
                dtime = variables.get('dtime')
            except ValueError:
                #print('dtime not added/updated')
                pass
            #else:
                #print('dtime = ', dtime)
        if 'dpept' in variables:
            try:
                dpept = variables.get('dpept')
            except ValueError:
                #print('dpept not added/updated')
                pass
            #else:
                #print('dpept = ', dpept)
        if 'dprot' in variables:
            try:
                dprot = variables.get('dprot')
            except ValueError:
                #print('dprot not added/updated')
                pass
            #else:
                #print('dprot = ', dprot)
        if 'renumdict' in variables:
            try:
                renumdict = variables.get('renumdict')
            except ValueError:
                #print('renumdict not added/updated')
                pass
            #else:
                #print('renumdict = ', renumdict)
        if 'muthandling' in variables:
            try:
                muthandling = variables.get('muthandling')
            except ValueError:
                #print('muthandling not added/updated')
                pass
            #else:
                #print('muthandling = ', muthandling)
        if 'slist' in variables:
            try:
                slist = variables.get('slist')
            except ValueError:
                #print('slist not added/updated')
                pass
            #else:
                #print('slist = ', slist)
        if 'annotateHM' in variables:
            try:
                annotateHM = variables.get('annotateHM')
            except ValueError:
                #print('annotateHM not added/updated')
                pass
            #else:
                #print('annotateHM = ', annotateHM)
        if 'chaindictuse' in variables:
            try:
                chaindictuse = variables.get('chaindictuse')
            except ValueError:
                #print('chaindictuse not added/updated')
                pass
            #else:
                #print('chaindictuse = ', chaindictuse)
        if 'numrenum' in variables:
            try:
                numrenum = variables.get('numrenum')
            except ValueError:
                #print('numrenum not added/updated')
                pass
            #else:
                #print('numrenum = ', numrenum)
        if 'key' in variables:
            try:
                key = variables.get('key')
            except ValueError:
                #print('key not added/updated')
                pass
            #else:
                #print('key = ', key)
        if 'value' in variables:
            try:
                value = variables.get('value')
            except ValueError:
                #print('value not added/updated')
                pass
            #else:
                #print('value = ', value)
        if 'numchaindict' in variables:
            try:
                numchaindict = variables.get('numchaindict')
            except ValueError:
                #print('numchaindict not added/updated')
                pass
            #else:
                #print('numchaindict = ', numchaindict)
        if 'chainkey' in variables:
            try:
                chainkey = variables.get('chainkey')
            except ValueError:
                #print('chainkey not added/updated')
                pass
            #else:
                #print('chainkey = ', chainkey)
        if 'chainvalue' in variables:
            try:
                chainvalue = variables.get('chainvalue')
            except ValueError:
                #print('chainvalue not added/updated')
                pass
            #else:
                #print('chainvalue = ', chainvalue)
        if 'numMutProt' in variables:
            try:
                numMutProt = variables.get('numMutProt')
            except ValueError:
                #print('numMutProt not added/updated')
                pass
            #else:
                #print('numMutProt = ', numMutProt)
        if 'keyMutProt' in variables:
            try:
                keyMutProt = variables.get('keyMutProt')
            except ValueError:
                #print('keyMutProt not added/updated')
                pass
            #else:
                #print('keyMutProt = ', keyMutProt)
        if 'valueMutProt' in variables:
            try:
                valueMutProt = variables.get('valueMutProt')
            except ValueError:
                #print('valueMutProt not added/updated')
                pass
            #else:
                #print('valueMutProt = ', valueMutProt)
        if 'numMutRes' in variables:
            try:
                numMutRes = variables.get('numMutRes')
            except ValueError:
                #print('numMutRes not added/updated')
                pass
            #else:
                #print('numMutRes = ', numMutRes)
        if 'resPos_raw' in variables:
            try:
                resPos_raw = variables.get('resPos_raw')
            except ValueError:
                #print('resPos_raw not added/updated')
                pass
            #else:
                #print('resPos_raw = ', resPos_raw)
        if 'WTRes' in variables:
            try:
                WTRes = variables.get('WTRes')
            except ValueError:
                #print('WTRes not added/updated')
                pass
            #else:
                #print('WTRes = ', WTRes)
        if 'num_timedao' in variables:
            try:
                num_timedao = variables.get('num_timedao')
            except ValueError:
                #print('num_timedao not added/updated')
                pass
            #else:
                #print('num_timedao = ', num_timedao)
        if 'value' in variables:
            try:
                value = variables.get('value')
            except ValueError:
                #print('value not added/updated')
                pass
            #else:
                #print('value = ', value)
        if 'state_of_interest' in variables:
            try:
                state_of_interest = variables.get('state_of_interest')
            except ValueError:
                #print('state_of_interest not added/updated')
                pass
            #else:
                #print('state_of_interest = ', state_of_interest)
        if 'ref_time' in variables:
            try:
                ref_time = variables.get('ref_time')
            except ValueError:
                #print('ref_time not added/updated')
                pass
            #else:
                #print('ref_time = ', ref_time)
        if 'numNewTimes' in variables:
            try:
                numNewTimes = variables.get('numNewTimes')
            except ValueError:
                #print('numNewTimes not added/updated')
                pass
            #else:
                #print('numNewTimes = ', numNewTimes)
        if 'valueNT' in variables:
            try:
                valueNT = variables.get('valueNT')
            except ValueError:
                #print('valueNT not added/updated')
                pass
            #else:
                #print('valueNT = ', valueNT)
        if 'num_dpept' in variables:
            try:
                num_dpept = variables.get('num_dpept')
            except ValueError:
                #print('num_dpept not added/updated')
                pass
            #else:
                #print('num_dpept = ', num_dpept)
        if 'dpeppro' in variables:
            try:
                dpeppro = variables.get('dpeppro')
            except ValueError:
                #print('dpeppro not added/updated')
                pass
            #else:
                #print('dpeppro = ', dpeppro)
        if 'dpepst' in variables:
            try:
                dpepst = variables.get('dpepst')
            except ValueError:
                #print('dpepst not added/updated')
                pass
            #else:
                #print('dpepst = ', dpepst)
        if 'dpepend' in variables:
            try:
                dpepend = variables.get('dpepend')
            except ValueError:
                #print('dpepend not added/updated')
                pass
            #else:
                #print('dpepend = ', dpepend)
        if 'dpro1' in variables:
            try:
                dpro1 = variables.get('dpro1')
            except ValueError:
                #print('dpro1 not added/updated')
                pass
            #else:
                #print('dpro1 = ', dpro1)
        if 'dpro2' in variables:
            try:
                dpro2 = variables.get('dpro2')
            except ValueError:
                #print('dpro2 not added/updated')
                pass
            #else:
                #print('dpro2 = ', dpro2)
        if 'dpro3' in variables:
            try:
                dpro3 = variables.get('dpro3')
            except ValueError:
                #print('dpro3 not added/updated')
                pass
            #else:
                #print('dpro3 = ', dpro3)
        if 'dpro3' in variables:
            try:
                dpro3 = variables.get('dpro3')
            except ValueError:
                #print('dpro3 not added/updated')
                pass
            #else:
                #print('dpro3 = ', dpro3)
        if 'dpro4' in variables:
            try:
                dpro4 = variables.get('dpro4')
            except ValueError:
                #print('dpro4 not added/updated')
                pass
            #else:
                #print('dpro4 = ', dpro4)
        if 'dpro5' in variables:
            try:
                dpro5 = variables.get('dpro5')
            except ValueError:
                #print('dpro5 not added/updated')
                pass
            #else:
                #print('dpro5 = ', dpro5)
        if 'dpro6' in variables:
            try:
                dpro6 = variables.get('dpro6')
            except ValueError:
                #print('dpro6 not added/updated')
                pass
            #else:
                #print('dpro6 = ', dpro6)
        if 's1' in variables:
            try:
                s1 = variables.get('s1')
            except ValueError:
                #print('s1 not added/updated')
                pass
            #else:
                #print('s1 = ', s1)
        if 's2' in variables:
            try:
                s2 = variables.get('s2')
            except ValueError:
                #print('s2 not added/updated')
                pass
            #else:
                #print('s2 = ', s2)
        if 'separate_plots_pls' in variables:
            try:
                separate_plots_pls = variables.get('separate_plots_pls')
            except ValueError:
                #print('separate_plots_pls not added/updated')
                pass
            #else:
                #print('separate_plots_pls = ', separate_plots_pls)
        if 'editing_fig_title' in variables:
            try:
                editing_fig_title = variables.get('editing_fig_title')
            except ValueError:
                #print('editing_fig_title not added/updated')
                pass
            #else:
                #print('editing_fig_title = ', editing_fig_title)
        if 'editing_xaxis' in variables:
            try:
                editing_xaxis = variables.get('editing_xaxis')
            except ValueError:
                #print('editing_xaxis not added/updated')
                pass
            #else:
                #print('editing_xaxis = ', editing_xaxis)
        if 'new_x_axis_title' in variables:
            try:
                new_x_axis_title = variables.get('new_x_axis_title')
            except ValueError:
                #print('new_x_axis_title not added/updated')
                pass
            #else:
                #print('new_x_axis_title = ', new_x_axis_title)
        if 'editing_yaxis' in variables:
            try:
                editing_yaxis = variables.get('editing_yaxis')
            except ValueError:
                #print('editing_yaxis not added/updated')
                pass
            #else:
                #print('editing_yaxis = ', editing_yaxis)
        if 'new_y_axis_title' in variables:
            try:
                new_y_axis_title = variables.get('new_y_axis_title')
            except ValueError:
                #print('new_y_axis_title not added/updated')
                pass
            #else:
                #print('new_y_axis_title = ', new_y_axis_title)
        if 'dif_dpi' in variables:
            try:
                dif_dpi = variables.get('dif_dpi')
            except ValueError:
                #print('dif_dpi not added/updated')
                pass
            #else:
                #print('dif_dpi = ', dif_dpi)
        if 'output_bitmap_dpi' in variables:
            try:
                output_bitmap_dpi = variables.get('output_bitmap_dpi')
            except ValueError:
                #print('output_bitmap_dpi not added/updated')
                pass
            #else:
                #print('output_bitmap_dpi = ', output_bitmap_dpi)
        if 'altvolc' in variables:
            try:
                altvolc = variables.get('altvolc')
            except ValueError:
                #print('altvolc not added/updated')
                pass
            #else:
                #print('altvolc = ', altvolc)
        if 'altvolccol' in variables:
            try:
                altvolccol = variables.get('altvolccol')
            except ValueError:
                #print('altvolccol not added/updated')
                pass
            #else:
                #print('altvolccol = ', altvolccol)
        if 'title' in variables:
            try:
                title = variables.get('title')
            except ValueError:
                #print('title not added/updated')
                pass
            #else:
                #print('title = ', title)
        if 'all_in_min' in variables:
            try:
                all_in_min = variables.get('all_in_min')
            except ValueError:
                #print('all_in_min not added/updated')
                pass
            #else:
                #print('all_in_min = ', all_in_min)
        if 'what_t_unit' in variables:
            try:
                what_t_unit = variables.get('what_t_unit')
            except ValueError:
                #print('what_t_unit not added/updated')
                pass
            #else:
                #print('what_t_unit = ', what_t_unit)
        if 'usehmthick' in variables:
            try:
                usehmthick = variables.get('usehmthick')
            except ValueError:
                #print('usehmthick not added/updated')
                pass
            #else:
                #print('usehmthick = ', usehmthick)
        if 'usestatesep' in variables:
            try:
                usestatesep = variables.get('usestatesep')
            except ValueError:
                #print('usestatesep not added/updated')
                pass
            #else:
                #print('usestatesep = ', usestatesep)
        if 'usebordchange' in variables:
            try:
                usebordchange = variables.get('usebordchange')
            except ValueError:
                #print('usebordchange not added/updated')
                pass
            #else:
                #print('usebordchange = ', usebordchange)
        if 'usehmtickchange' in variables:
            try:
                usehmtickchange = variables.get('usehmtickchange')
            except ValueError:
                #print('usehmtickchange not added/updated')
                pass
            #else:
                #print('usehmtickchange = ', usehmtickchange)
        if 'hmspacerthick' in variables:
            try:
                hmspacerthick = variables.get('hmspacerthick')
            except ValueError:
                #print('hmspacerthick not added/updated')
                pass
            #else:
                #print('hmspacerthick = ', hmspacerthick)
        if 'hmcolordivide' in variables:
            try:
                hmcolordivide = variables.get('hmcolordivide')
            except ValueError:
                #print('hmcolordivide not added/updated')
                pass
            #else:
                #print('hmcolordivide = ', hmcolordivide)
        if 'hmsepthick' in variables:
            try:
                hmsepthick = variables.get('hmsepthick')
            except ValueError:
                #print('hmsepthick not added/updated')
                pass
            #else:
                #print('hmsepthick = ', hmsepthick)
        if 'hmsepcolor' in variables:
            try:
                hmsepcolor = variables.get('hmsepcolor')
            except ValueError:
                #print('hmsepcolor not added/updated')
                pass
            #else:
                #print('hmsepcolor = ', hmsepcolor)
        if 'hmbordthick' in variables:
            try:
                hmbordthick = variables.get('hmbordthick')
            except ValueError:
                #print('hmbordthick not added/updated')
                pass
            #else:
                #print('hmbordthick = ', hmbordthick)
        if 'hmbordcolor' in variables:
            try:
                hmbordcolor = variables.get('hmbordcolor')
            except ValueError:
                #print('hmbordcolor not added/updated')
                pass
            #else:
                #print('hmbordcolor = ', hmbordcolor)
        if 'hmtw' in variables:
            try:
                hmtw = variables.get('hmtw')
            except ValueError:
                #print('hmtw not added/updated')
                pass
            #else:
                #print('hmtw = ', hmtw)
        if 'hmtl' in variables:
            try:
                hmtl = variables.get('hmtl')
            except ValueError:
                #print('hmtl not added/updated')
                pass
            #else:
                #print('hmtl = ', hmtl)
        if 'hmtickcolor' in variables:
            try:
                hmtickcolor = variables.get('hmtickcolor')
            except ValueError:
                #print('hmtickcolor not added/updated')
                pass
            #else:
                #print('hmtickcolor = ', hmtickcolor)
        if 'hmtickcolorlabel' in variables:
            try:
                hmtickcolorlabel = variables.get('hmtickcolorlabel')
            except ValueError:
                #print('hmtickcolorlabel not added/updated')
                pass
            #else:
                #print('hmtickcolorlabel = ', hmtickcolorlabel)
        if 'usealtpad' in variables:
            try:
                usealtpad = variables.get('usealtpad')
            except ValueError:
                #print('usealtpad not added/updated')
                pass
            #else:
                #print('usealtpad = ', usealtpad)
        if 'padthick' in variables:
            try:
                padthick = variables.get('padthick')
            except ValueError:
                #print('padthick not added/updated')
                pass
            #else:
                #print('padthick = ', padthick)
        if 'spadthick' in variables:
            try:
                spadthick = variables.get('spadthick')
            except ValueError:
                #print('spadthick not added/updated')
                pass
            #else:
                #print('spadthick = ', spadthick)
        if 'domainlabel' in variables:
            try:
                domainlabel = variables.get('domainlabel')
            except ValueError:
                #print('domainlabel not added/updated')
                pass
            #else:
                #print('domainlabel = ', domainlabel)
        if 'numdomain' in variables:
            try:
                numdomain = variables.get('numdomain')
            except ValueError:
                #print('numdomain not added/updated')
                pass
            #else:
                #print('numdomain = ', numdomain)
        if 'domainName' in variables:
            try:
                domainName = variables.get('domainName')
            except ValueError:
                #print('domainName not added/updated')
                pass
            #else:
                #print('domainName = ', domainName)
        if 'dompepst' in variables:
            try:
                dompepst = variables.get('dompepst')
            except ValueError:
                #print('dompepst not added/updated')
                pass
            #else:
                #print('dompepst = ', dompepst)
        if 'dompepend' in variables:
            try:
                dompepend = variables.get('dompepend')
            except ValueError:
                #print('dompepend not added/updated')
                pass
            #else:
                #print('dompepend = ', dompepend)
        if 'domColour' in variables:
            try:
                domColour = variables.get('domColour')
            except ValueError:
                #print('domColour not added/updated')
                pass
            #else:
                #print('domColour = ', domColour)
        if 'newtitle' in variables:
            try:
                newtitle = variables.get('newtitle')
            except ValueError:
                #print('newtitle not added/updated')
                pass
            #else:
                #print('newtitle = ', newtitle)
        if 'altwoodsthick' in variables:
            try:
                altwoodsthick = variables.get('altwoodsthick')
            except ValueError:
                #print('altwoodsthick not added/updated')
                pass
            #else:
                #print('altwoodsthick = ', altwoodsthick)
        if 'woodscolthick' in variables:
            try:
                woodscolthick = variables.get('woodscolthick')
            except ValueError:
                #print('woodscolthick not added/updated')
                pass
            #else:
                #print('woodscolthick = ', woodscolthick)
        if 'woodsbackthick' in variables:
            try:
                woodsbackthick = variables.get('woodsbackthick')
            except ValueError:
                #print('woodsbackthick not added/updated')
                pass
            #else:
                #print('woodsbackthick = ', woodsbackthick)
        if 'woodsplotYbyglobal' in variables:
            try:
                woodsplotYbyglobal = variables.get('woodsplotYbyglobal')
            except ValueError:
                #print('woodsplotYbyglobal not added/updated')
                pass
            #else:
                #print('woodsplotYbyglobal = ', woodsplotYbyglobal)
        if 'download_pymol' in variables:
            try:
                download_pymol = variables.get('download_pymol')
            except ValueError:
                #print('download_pymol not added/updated')
                pass
            #else:
                #print('download_pymol = ', download_pymol)
        if 'download_chimera' in variables:
            try:
                download_chimera = variables.get('download_chimera')
            except ValueError:
                #print('download_chimera not added/updated')
                pass
            #else:
                #print('download_chimera = ', download_chimera)
        if 'fontcolor' in variables:
            try:
                fontcolor = variables.get('fontcolor')
            except ValueError:
                #print('fontcolor not added/updated')
                pass
            #else:
                #print('fontcolor = ', fontcolor)
        if 'custom_colors' in variables:
            try:
                custom_colors = variables.get('custom_colors')
            except ValueError:
                #print('custom_colors not added/updated')
                pass
            #else:
                #print('custom_colors = ', custom_colors)
        if 'custom_bounds' in variables:
            try:
                custom_bounds = variables.get('custom_bounds')
            except ValueError:
                #print('custom_bounds not added/updated')
                pass
            #else:
                #print('custom_bounds = ', custom_bounds)
        if 'negative_bounds' in variables:
            try:
                negative_bounds = variables.get('negative_bounds')
            except ValueError:
                #print('negative_bounds not added/updated')
                pass
            #else:
                #print('negative_bounds = ', negative_bounds)
        if 'positive_bounds' in variables:
            try:
                positive_bounds = variables.get('positive_bounds')
            except ValueError:
                #print('positive_bounds not added/updated')
                pass
            #else:
                #print('positive_bounds = ', positive_bounds)
        if 'renumbering_dict' in variables:
            try:
                renumbering_dict = variables.get('renumbering_dict')
            except ValueError:
                #print('renumbering_dict not added/updated')
                pass
            #else:
                #print('renumbering_dict = ', renumbering_dict)
        if 'chain_dict' in variables:
            try:
                chain_dict = variables.get('chain_dict')
            except ValueError:
                #print('chain_dict not added/updated')
                pass
            #else:
                #print('chain_dict = ', chain_dict)
        if 'mut_id_dict' in variables:
            try:
                mut_id_dict = variables.get('mut_id_dict')
            except ValueError:
                #print('mut_id_dict not added/updated')
                pass
            #else:
                #print('mut_id_dict = ', mut_id_dict)
        if 'mutation_dict' in variables:
            try:
                mutation_dict = variables.get('mutation_dict')
            except ValueError:
                #print('mutation_dict not added/updated')
                pass
            #else:
                #print('mutation_dict = ', mutation_dict)
        if 'correlation_data' in variables:
            try:
                correlation_data = variables.get('correlation_data')
            except ValueError:
                #print('correlation_data not added/updated')
                pass
            #else:
                #print('correlation_data = ', correlation_data)
        if 'drop_times' in variables:
            try:
                drop_times = variables.get('drop_times')
            except ValueError:
                #print('drop_times not added/updated')
                pass
            #else:
                #print('drop_times = ', drop_times)
        if 'new_times' in variables:
            try:
                new_times = variables.get('new_times')
            except ValueError:
                #print('new_times not added/updated')
                pass
            #else:
                #print('new_times = ', new_times)
        if 'drop_pept' in variables:
            try:
                drop_pept = variables.get('drop_pept')
            except ValueError:
                #print('drop_pept not added/updated')
                pass
            #else:
                #print('drop_pept = ', drop_pept)
        if 'drop_prot' in variables:
            try:
                drop_prot = variables.get('drop_prot')
            except ValueError:
                #print('drop_prot not added/updated')
                pass
            #else:
                #print('drop_prot = ', drop_prot)
        if 'state1_list' in variables:
            try:
                state1_list = variables.get('state1_list')
            except ValueError:
                #print('state1_list not added/updated')
                pass
            #else:
                #print('state1_list = ', state1_list)
        if 'state2_list' in variables:
            try:
                state2_list = variables.get('state2_list')
            except ValueError:
                #print('state2_list not added/updated')
                pass
            #else:
                #print('state2_list = ', state2_list)
        if 'plot_separate' in variables:
            try:
                plot_separate = variables.get('plot_separate')
            except ValueError:
                #print('plot_separate not added/updated')
                pass
            #else:
                #print('plot_separate = ', plot_separate)
        if 'plot_stacked' in variables:
            try:
                plot_stacked = variables.get('plot_stacked')
            except ValueError:
                #print('plot_stacked not added/updated')
                pass
            #else:
                #print('plot_stacked = ', plot_stacked)
        if 'output_bitmap_format' in variables:
            try:
                output_bitmap_format = variables.get('output_bitmap_format')
            except ValueError:
                #print('output_bitmap_format not added/updated')
                pass
            #else:
                #print('output_bitmap_format = ', output_bitmap_format)
        if 'color_sections' in variables:
            try:
                color_sections = variables.get('color_sections')
            except ValueError:
                #print('color_sections not added/updated')
                pass
            #else:
                #print('color_sections = ', color_sections)
        if 'hmlcolor' in variables:
            try:
                hmlcolor = variables.get('hmlcolor')
            except ValueError:
                #print('hmlcolor not added/updated')
                pass
            #else:
                #print('hmlcolor = ', hmlcolor)
        if 'hmsepcolorc' in variables:
            try:
                hmsepcolorc = variables.get('hmsepcolorc')
            except ValueError:
                #print('hmsepcolorc not added/updated')
                pass
            #else:
                #print('hmsepcolorc = ', hmsepcolorc)
        if 'hmbordcolorc' in variables:
            try:
                hmbordcolorc = variables.get('hmbordcolorc')
            except ValueError:
                #print('hmbordcolorc not added/updated')
                pass
            #else:
                #print('hmbordcolorc = ', hmbordcolorc)
        if 'hmtcolor' in variables:
            try:
                hmtcolor = variables.get('hmtcolor')
            except ValueError:
                #print('hmtcolor not added/updated')
                pass
            #else:
                #print('hmtcolor = ', hmtcolor)
        if 'hmtcolor_labels' in variables:
            try:
                hmtcolor_labels = variables.get('hmtcolor_labels')
            except ValueError:
                #print('hmtcolor_labels not added/updated')
                pass
            #else:
                #print('hmtcolor_labels = ', hmtcolor_labels)
        if 'scattercolor' in variables:
            try:
                scattercolor = variables.get('scattercolor')
            except ValueError:
                #print('scattercolor not added/updated')
                pass
            #else:
                #print('scattercolor = ', scattercolor)
        if 'blockposwoods' in variables:
            try:
                blockposwoods = variables.get('blockposwoods')
            except ValueError:
                #print('blockposwoods not added/updated')
                pass
            #else:
                #print('blockposwoods = ', blockposwoods)
        if 'wpadthick' in variables:
            try:
                wpadthick = variables.get('wpadthick')
            except ValueError:
                #print('wpadthick not added/updated')
                pass
            #else:
                #print('wpadthick = ', wpadthick)
        if separate_plots_pls == 1:
            plot_separate = 1
            plot_stacked = 0
        else:    
            plot_stacked  = 1
            plot_separate = 0
        if annotateHM == 1:
            addvalHM = True
        else: 
            addvalHM = False
    
    if uploaded_settings != 1:
        alt_c2 = float(request.form.get('alt_c2',0))
        if alt_c2 == 1:
            c2 = ensure_hex_color(request.form.get('alt_c2_col', '#FFFFFF'))
        else:
            c2 = whitec
        
        alt_cmissing = float(request.form.get('alt_cmissing',0))
        if alt_cmissing == 1:
            c_missing = ensure_hex_color(request.form.get('alt_cmissing_col','#bdbdbd'))
        else:
            c_missing = '#bdbdbd'
            
        #For Woods Plots
        color_by_heatmap = float(request.form.get('colorbyheatmap',0))
        colcutopt = float(request.form.get('colcutopt',0))
        if colcutopt == 1:
            colcutoff = float(request.form.get('colcutoff',0.5))
        else:
            colcutoff = 0.5
        woodsdimen = float(request.form.get('woodsdimen',0))
        if woodsdimen == 1:
            woodsx = float(request.form.get('woodsx',20))
            woodsy = float(request.form.get('woodsy',6))
        else:
            woodsx = 20
            woodsy = 6
        woodscol = float(request.form.get('woodscol',0))
        if woodscol == 1:
            woodscolpos = request.form.get('woodscolpos','#FF0000')
            woodscolneu = request.form.get('woodscolneu',whitec)
            woodscolneg = request.form.get('woodscolneg','#0000FF')
            woodscolpos = ensure_hex_color(woodscolpos)
            woodscolneu = ensure_hex_color(woodscolneu)
            woodscolneg = ensure_hex_color(woodscolneg)
        else:
            woodscolpos = '#FF0000'
            woodscolneu = whitec
            woodscolneg = '#0000FF'
        nolines = float(request.form.get('nolines', 0))
        #Set font sizes
        font_size = float(request.form['font_size'])
        font_size_title = float(request.form['font_size_title'])
        font_size_ticklabel = float(request.form['font_size_ticklabel'])
        fontcolor = blackc
        fontcolorc = float(request.form.get('fontcolor',2))
        if fontcolorc == 1: #White
            fontcolor = whitec 
        elif fontcolorc == 2: #Black
            fontcolor = blackc
        elif fontcolorc == 3: #Blue
            fontcolor = dbluec
        elif fontcolorc == 4: #Red
            fontcolor = dredc
        elif fontcolorc == 5: #Grey
            fontcolor = greyc
        #Set horizontal or vertical plots
        h_or_v = float(request.form['h_or_v'])
        if h_or_v == 1:
            plot_v = 0
            plot_h = 1
            plot_w = 0
            plot_volc = 0
            plot_l = 0
        elif h_or_v == 2:
            plot_v = 1
            plot_h = 0
            plot_w = 0
            plot_volc = 0
            plot_l = 0
        elif h_or_v == 3:
            plot_v = 0
            plot_h = 0
            plot_w = 1
            plot_volc = 0
            plot_l = 0
        elif h_or_v == 4:
            plot_v = 0
            plot_h = 0
            plot_w = 0
            plot_volc = 1
            plot_l = 0
        elif h_or_v == 6:
            plot_v = 0
            plot_h = 0
            plot_w = 0
            plot_volc = 0
            plot_l = 1
        if plot_volc == 1:
            scatter_plot = 1
        elif plot_volc == 0:
            scatter_plot = 0
        zerobound = float(request.form.get('zerobound', 0))
        funkybound = float(request.form.get('change_bounds_abs',3))
        globalmax_delta = float(request.form.get('global_max', 0))
        #relativeUptakeCalc = float(request.form.get('relativeUptakeCalc', 0))
        relativeUptakeCalc = 0
        p_threshold = request.form.get('pthresh', 0.05)
        if p_threshold == '':
            p_threshold = 0.05
        if p_threshold != 0.05 and p_threshold != '':
            p_threshold = float(p_threshold)
        #Choose between setting MaxRange/#Shades and CustomColours/Bounds
        custb = float(request.form['option'])
        if custb == 3: # GET MAX RANGE AND NUM SHADES
            max_range = float(request.form['max_range']) 
            max_range = round(max_range, 1)
            num_shades = float(request.form['num_shades'])
            num_shades = round(num_shades) + 1
            alt_col = float(request.form.get('alt_col', 0))
            if alt_col == 1:
                if max_range <= 2:
                    c1 = request.form.get('negcolalt',bluec).strip()
                    c1 = ensure_hex_color(c1)
                    c3 = request.form.get('poscolalt',redc).strip()
                    c3 = ensure_hex_color(c3)
                else:
                    c1 = request.form.get('negcolalt',dbluec).strip()
                    c1 = ensure_hex_color(c1)
                    c3 = request.form.get('poscolalt',dredc).strip()
                    c3 = ensure_hex_color(c3)
        elif custb == 2: # AUTO SET
            num_shades = 7
            alt_col = float(request.form.get('alt_col', 0))
            if alt_col == 0:
                if funkybound == 1:
                    custom_colors = [c2, '#FC9272', '#FB6A4A', '#EF3B2C', '#CB181D', '#A50F15', '#67000D']
                elif funkybound == 2:
                    custom_colors = ['#023858', '#045A8D', '#0570B0', '#3690C0', '#74A9CF', '#A6BDDB', c2]
                else:
                    if zerobound == 0:
                        custom_colors = ['#023858', '#045A8D', '#0570B0', '#3690C0', '#74A9CF', '#A6BDDB', c2, '#FC9272', '#FB6A4A', '#EF3B2C', '#CB181D', '#A50F15', '#67000D']
                    elif zerobound == 1:
                        custom_colors = ['#023858', '#045A8D', '#0570B0', '#3690C0', '#74A9CF', '#A6BDDB', '#FC9272', '#FB6A4A', '#EF3B2C', '#CB181D', '#A50F15', '#67000D']
            elif alt_col == 1:
                ncolforalt = 7
                if funkybound == 1:
                    c3 = request.form.get('poscolalt', dredc).strip()
                    c3 = ensure_hex_color(c3)
                    #print('c3:', c3)
                    clist = [colorFader(c2,c3,x/ncolforalt) for x in range(ncolforalt)]
                    custom_colors = [c2] + clist[1:]
                elif funkybound == 2:
                    c1 = request.form.get('negcolalt', dbluec).strip()
                    c1 = ensure_hex_color(c1)
                    ncolforaltneg = ncolforalt - 1
                    #print('c1:', c1)
                    clist = [colorFader(c1,c2,x/ncolforaltneg) for x in range(ncolforaltneg)]
                    custom_colors = clist + [c2]
                else:
                    c1 = request.form.get('negcolalt', dbluec).strip()
                    c3 = request.form.get('poscolalt', dredc).strip()
                    c1 = ensure_hex_color(c1)
                    c3 = ensure_hex_color(c3)
                    #print('c1:', c1)
                    #print('c3:', c3)
                    ncolforaltneg = ncolforalt - 1
                    positive_c = [colorFader(c2,c3,x/ncolforalt) for x in range(ncolforalt)]
                    print('positive_c:', positive_c)
                    negative_c = [colorFader(c1,c2,x/ncolforaltneg) for x in range(ncolforaltneg)]
                    print('negative_c:', negative_c)
                    custom_colors = []
                    if zerobound == 0:
                        custom_colors.extend(negative_c)
                        #custom_colors.append(c2)
                        custom_colors.extend(positive_c)
                    elif zerobound == 1:
                        custom_colors.extend(negative_c)
                        custom_colors.extend(positive_c)
                        while '#ffffff' in custom_colors:
                            custom_colors.remove('#ffffff')
                        while '#FFFFFF' in custom_colors:
                            custom_colors.remove('#FFFFFF')
                print('custom colours, custb 2:', custom_colors)
                num_shades = 8
        elif custb == 1: # CHOOSE COL, SET BOUNDS
            max_range = None
            custom_bounds = []
            negative_bounds = []
            positive_bounds = []
            for i in range(1, 9):
                nb = request.form.get(f'inputbn{i}', '').strip()
                pb = request.form.get(f'inputbp{i}', '').strip() 
                if nb:
                    try:
                        nb = float(nb)
                        negative_bounds.append(nb)
                    except ValueError:
                        pass
                if pb:
                    try:
                        pb = float(pb)
                        positive_bounds.append(pb) 
                    except ValueError:
                        pass 
            if zerobound == 1:
                custom_bounds.append(0)
            custom_bounds.extend(negative_bounds)
            custom_bounds.extend(positive_bounds)
            custom_bounds = sorted(custom_bounds)
            #print('Sorted custom bounds:', custom_bounds)
            numinputs = len(custom_bounds) // 2
            if zerobound == 1:
                numinputs == numinputs + 1
            #Gather Colour Choices 
            neg_col = float(request.form['neg_col'])
            pos_col = float(request.form['pos_col'])
            num_col = 2*numinputs - 1
            numbshades_pos = len(positive_bounds) - 1
            numbshades_neg = len(negative_bounds) - 1
            numbshades = numinputs - 1
            #Define negative colours - dark, medium, and light
            if neg_col == 1: #red
                if max(abs(x) for x in custom_bounds) >= 2:
                    dcoln = dredc
                else:
                    dcoln = redc
            elif neg_col == 2: #blue
                if max(abs(x) for x in custom_bounds) >= 2:
                    dcoln = dbluec
                else: 
                    dcoln = bluec
            elif neg_col == 3: #green
                dcoln = greenc
            elif neg_col == 4: #yellow
                dcoln = yellowc
            elif neg_col == 5: #purple
                dcoln = purplec
            elif neg_col == 6: #black
                dcoln = blackc
            #Define positive colours - dark, medium, and light
            if pos_col == 1: #blue
                if max(abs(x) for x in custom_bounds) >= 2:
                    dcolp = dbluec
                else:
                    dcolp = bluec
            elif pos_col == 2: #red
                if max(abs(x) for x in custom_bounds) >= 2:
                    dcolp = dredc
                else:
                    dcolp = redc
            elif pos_col == 3: #green
                dcolp = greenc
            elif pos_col == 4: #yellow
                dcolp = yellowc
            elif pos_col == 5: #purple
                dcolp = purplec
            elif pos_col == 6: #black
                dcolp = blackc
            c1 = dcoln
            c2 = whitec
            c3 = dcolp
            col_mid = [mc.to_hex(c2)]
            cols1 = [colorFader(c1,c2,x/numbshades_neg) for x in range(numbshades_neg+1)]
            cols2 = [colorFader(c2,c3,x/numbshades_pos) for x in range(numbshades_pos+1)]
            if zerobound == 1:
                custom_colors = cols1[:-1] +cols2[1:]
            if zerobound == 0:
                custom_colors = cols1[:-1] + col_mid + cols2[1:]
        elif custb == 4: #INPUT COL AND BOUNDS
            max_range = None
            custom_bounds = []
            negative_bounds = []
            positive_bounds = []
            for i in range(1, 9):
                nb = request.form.get(f'inputbn{i}4', '').strip()
                pb = request.form.get(f'inputbp{i}4', '').strip() 
                if nb:
                    try:
                        nb = float(nb)
                        negative_bounds.append(nb)
                    except ValueError:
                        pass
                if pb:
                    try:
                        pb = float(pb)
                        positive_bounds.append(pb) 
                    except ValueError:
                        pass 
            if zerobound == 1:
                custom_bounds.append(0)
            custom_bounds.extend(negative_bounds)
            custom_bounds.extend(positive_bounds)
            custom_bounds = sorted(custom_bounds)
            #print('Sorted custom bounds:', custom_bounds)
            numinputs = len(custom_bounds) // 2
            numbshades_pos = len(positive_bounds) - 1
            numbshades_neg = len(negative_bounds) - 1
            numbshades = numinputs - 1
            colorchoicetype = float(request.form['optionc'])
            if colorchoicetype == 1:
                negshadein = str(request.form['negcolorid']).strip()
                posshadein = str(request.form['poscolorid']).strip()
                c1 = ensure_hex_color(negshadein)
                c2 = whitec
                c3 = ensure_hex_color(posshadein)
                col_mid = [mc.to_hex(c2)]
                if zerobound == 1:
                    cols1 = [colorFader(c1,c2,x/numbshades_neg) for x in range(numbshades_neg)]
                    cols2 = [colorFader(c2,c3,x/numbshades_pos) for x in range(numbshades_pos)]
                    custom_colors = cols1[:-1] +cols2[1:]
                if zerobound == 0:
                    cols1 = [colorFader(c1,c2,x/numbshades_neg) for x in range(numbshades_neg+1)]
                    cols2 = [colorFader(c2,c3,x/numbshades_pos) for x in range(numbshades_pos+1)]
                    custom_colors = cols1[:-1] + col_mid + cols2[1:]
            elif colorchoicetype ==2: 
                custom_colors = []
                nc_list = []
                pc_list = []
                col_mid = mc.to_rgb(whitec) 
                for i in range(1, 9):
                    nc = request.form.get(f'ninputcol{i}', '').strip()
                    pc = request.form.get(f'pinputcol{i}', '').strip()
                    nc = ensure_hex_color(nc)
                    pc = ensure_hex_color(pc)
                    if nc:
                        try:
                            nc_rgb = mc.to_rgb(nc)
                            nc_list.append(nc_rgb)
                        except ValueError:
                            pass
                    if pc:
                        try:
                            pc_rgb = mc.to_rgb(pc)
                            pc_list.append(pc_rgb)
                        except ValueError:
                            pass
                for nc in nc_list:
                    custom_colors.append(nc)
                if zerobound == 0:
                    custom_colors.append(col_mid)
                for pc in pc_list:
                    custom_colors.append(pc)
        #Determine if Advanced Options were chosen, if yes import those values
        dtime = int(request.form.get('droptime', 0))
        dpept = int(request.form.get('droppept', 0))
        dprot = int(request.form.get('dropprot', 0))
        renumdict = int(request.form.get('renumdict', 0))
        #iddict = int(request.form.get('mutiddict', 0))
        #mutdict = int(request.form.get('mutdict', 0))
        muthandling = int(request.form.get('optionmut',0)) #2 = on
        slist = int(request.form.get('optionst', 1))
        annotateHM = int(request.form.get('annotateHM',0))
        chaindictuse = int(request.form.get('chaindictuse',0))
        if annotateHM == 1:
            addvalHM = True
        else: 
            addvalHM = False
        if renumdict == 1:
            #print('went into renum')
            renumbering_dict = {}
            numrenum = int(request.form.get('numrenum',0))
            for i in range(1, numrenum + 1):
                #print('Checking renum')
                key = str(request.form[f'key{i}']).strip()
                value = request.form[f'value{i}']
                # Only add the key-value pair if the key is not an empty string
                if key:
                    renumbering_dict[key] = value
            #print('renumbering_dict')
            #print(renumbering_dict)
        if chaindictuse == 1:
            print('went into chaindict')
            chain_dict = {}
            numchaindict = int(request.form.get('numchaindict',0))
            for i in range(1, numchaindict + 1):
                #print('Checking renum')
                chainkey = str(request.form[f'keychain{i}']).strip() #Protein Name (in input file)
                chainvalue = str(request.form[f'valuechain{i}']) #Chain Name (in PDB file)
                # Only add the key-value pair if the key is not an empty string
                if chainkey:
                    chain_dict[chainkey] = chainvalue
            #print('Chain Dict:', chain_dict)
        if muthandling == 2:
            mut_id_dict = {}
            numMutProt = int(request.form.get('numMutProt',0))
            for i in range(1, numMutProt + 1):
                keyMutProt = str(request.form[f'keyMutProt{i}']).strip() 
                valueMutProt = str(request.form[f'valueMutProt{i}']) 
                # Only add the key-value pair if the key is not an empty string
                if keyMutProt:
                    mut_id_dict[keyMutProt] = valueMutProt
            #print('mut_id_dict:', mut_id_dict)
            mutation_dict = {}
            numMutRes = int(request.form.get('numMutProt',0))
            for i in range(1, numMutRes + 1):
                resPos_raw = str(request.form[f'resPos{i}']).strip()
                if resPos_raw:
                    resPos = int(resPos_raw)
                    WTRes = str(request.form[f'WTRes{i}'])
                    mutation_dict[resPos] = WTRes
            correlation_data = []  # List to hold dictionaries for each correlated record
            for i in range(1, numMutProt + 1):
                keyMutProt = str(request.form.get(f'keyMutProt{i}', '')).strip()
                valueMutProt = str(request.form.get(f'valueMutProt{i}', '')).strip()
                resPos_raw = str(request.form.get(f'resPos{i}', '')).strip()
                WTRes = str(request.form.get(f'WTRes{i}', '')).strip()
                if keyMutProt and resPos_raw:  # Ensure both keyMutProt and resPos exist
                    correlation_data.append({
                        'keyMutProt': keyMutProt,
                        'valueMutProt': valueMutProt,
                        'resPos': int(resPos_raw),
                        'WTRes': WTRes
                    })
            #print('mutation_dict:', mutation_dict)
        if dtime == 1:
            num_timedao = int(request.form.get('numtimedao', 7))
            drop_times = []
            for i in range(0, num_timedao+7):
                value = request.form.get(f'dt{i}')
                if value:
                    try:
                        drop_times.append(float(value))
                    except ValueError:
                        pass
            #print('drop times list')
            #print(drop_times)
        if PulseLabellingNow == 2:
            state_of_interest = str(request.form.get('state_of_interest', 'NotProvided'))
            ref_time = float(request.form.get('ref_time',2020))
            numNewTimes = int(request.form.get('numTimePulse', 0))
            new_times = []
            for i in range(0, numNewTimes+7):
                valueNT = request.form.get(f'dNT{i}')
                #print(f'ValueNT{i}', '-', i, '-', 'added')
                if valueNT:
                    try:
                        new_times.append(float(valueNT))
                    except ValueError:
                        #print(valueNT, 'Not Added')
                        pass
            #print('new_times when called:', new_times)
        if dpept == 1:
            drop_pept = []
            num_dpept = int(request.form.get('numpeptd', 4))
            for i in range(1, num_dpept+4):
                dpeppro = request.form.get(f'dpeppro{i}', '')
                dpepst = float(request.form.get(f'dpepst{i}', 0))
                dpepend = float(request.form.get(f'dpepend{i}', 0))
                if dpeppro:
                    drop_pept.append((dpeppro, dpepst, dpepend))
        if dprot == 1:
            num_protdao = int(request.form.get('numprotdao', 7))
            print('num prot dao', num_protdao)
            drop_prot = []
            for i in range(0, num_protdao+7):
                value = request.form.get(f'dpro{i}')
                print('value:', i, value)
                if value:
                    try:
                        drop_prot.append(str(value))
                    except ValueError:
                        pass
            print('drop protein list:', drop_prot)
            #print(drop_prot)
        if slist == 2:
            numStatelist = int(request.form.get('numStatelist', 7))
            SpecifyStateProt = float(request.form.get('StateProt', 0)) ##why are you outputting as 0???
            print('Specify State Prot:', SpecifyStateProt)
            #print('starting state list pull')
            state1_list = []
            state2_list = []
            for i in range(0, numStatelist+3):
                s1 = str(request.form.get(f's1{i}', ''))
                s2 = str(request.form.get(f's2{i}', ''))
                if SpecifyStateProt == 1:
                    print('in SpecifyStateProt')
                    prSt1 = str(request.form.get(f'prSt1{i}', ''))
                    prSt2 = str(request.form.get(f'prSt2{i}', ''))
                    #print(f'prSt1{i}')
                    #print(f'prSt2{i}')
                    s1 += prSt1
                    s2 += prSt2
                # Only append non-empty strings to state1_list
                if s1 != '':
                    state1_list.append(s1)
                    # Assign s2 a default value only if it's empty
                    if s2 == '':
                        s2 = str(request.form.get('s21', ''))
                # Only append non-empty strings to state2_list
                if s2 != '':
                    state2_list.append(s2)
            # Ensure that state1_list and state2_list have the same length
            state1_list = state1_list[:len(state2_list)]
            print('state list 1:', state1_list)
            print('state list 2:', state2_list)
        separate_plots_pls = float(request.form.get('separate_plots','0'))
        if separate_plots_pls == 1:
            plot_separate = 1
            plot_stacked = 0
            #print('separate')
        else:    
            plot_stacked  = 1
            plot_separate = 0
            #print('not separate')
        editing_figtitle = float(request.form.get('usealtxaxis', 0))
        editing_xaxis = float(request.form.get('usealtxaxis', 0))
        if editing_xaxis == 1:
            new_x_axis_title = str(request.form.get('altxaxistit','Invalid Text'))
        else:
            if h_or_v == 1 or h_or_v == 2:
                new_x_axis_title = 'Peptides'
            if h_or_v == 3:
                new_x_axis_title = 'Position'
            if h_or_v == 4:
                new_x_axis_title = 'Difference in Deuterium Uptake (Da)'
                if expType == 3:
                    new_x_axis_title = 'Relative Deuterium Uptake'
                    if absoluteUptakeValues == 1:
                        new_x_axis_title = 'Absolute Deuterium Uptake (Da)'
            if h_or_v == 6:
                new_x_axis_title = 'Position'
        editing_yaxis = float(request.form.get('usealtyaxis', 0))
        if editing_yaxis == 1:
            new_y_axis_title = str(request.form.get('altyaxistit','Invalid Text'))
        else:
            if h_or_v == 1 or h_or_v == 2:
                new_y_axis_title = 'H/D exchange time'
            if h_or_v == 3:
                new_y_axis_title = 'Change in Deuterium Uptake (Da)'
                if expType == 3:
                    new_y_axis_title = 'Percent Relative Deuterium Uptake'
                    if absoluteUptakeValues == 1:
                        new_y_axis_title = 'Absolute Deuterium Uptake (Da)'
            if h_or_v == 4:
                new_y_axis_title = '-log$_{10}$(p-value)'
            if h_or_v == 6:
                new_y_axis_title = 'Peptide #'
        buffer = BytesIO()  # saves plot to BytesIO buffer
        output_csv_file  = r'HDX processed data.csv'
        output_pdf_file  = r"uploads/HDX heatmap.pdf"
        output_bitmap = 1
        output_bitmap_name = 'HDX heatmap'
        output_bitmap_dpi = 100
        dif_dpi = float(request.form.get('dif_dpi', 0))
        #print(dif_dpi)
        #print('this is if difdpi is working')
        if dif_dpi != 0:
            output_bitmap_dpi = float(request.form.get('dpi_in'))
        #print(output_bitmap_dpi)
        if PDFgeneration == 1:
            output_bitmap_format = 'pdf'
        else:
            output_bitmap_format = 'png'
    mutation_msg = 0
    split_outp_by_prot = 'all'


    #################################################################################################################################################################
    #################################################################################################################################################################
    #############################################################            Plotting Begins            #############################################################
    #################################################################################################################################################################
    #################################################################################################################################################################
                                                                                                        
    
    if custom_colors != None and custom_bounds != None:
        custom_bounds = list(set(custom_bounds))
        #print(len(custom_colors))
        #print(custom_colors)
        #print(custom_bounds)
        #print(len(custom_bounds))
        if len(custom_colors)>=1 and len(custom_bounds)>=1 and (len(custom_bounds) != len(custom_colors)+1):
            print('It appears that both custom_colors and custom_bounds lists are used together, but there is mismatch between expected numbers of entries.')
            print('custom_bounds list must have exactly one more entry than custom_colors. Aborting processing.')
            raise StopExecution
    if custom_colors != None and custom_bounds != None:
        if len(custom_colors)==0 and len(custom_bounds)==1:
            print('custom_bounds list must contain two or more numbers. Aborting processing.')
            raise StopExecution
    if custom_bounds != None:
        custom_bounds.sort()
    
    #if (not plot_h) and (not plot_v) and (not plot_w):
    #    print('Both variables \'plot_h\' and \'plot_v\' are set to zero. No heatmap plots will be produced.\n')
    #if (not plot_stacked) and (not plot_separate):
    #    print('Both variables \'plot_stacked\' and \'plot_separate\' are set to zero. No heatmap plots will be produced.\n')
    
    # Read the input and start processing:
    data_df = pd.read_csv(input_csv_file)
    
    if 'z' not in data_df.columns:
        #print('No charge data (\'z\' column) in the input file. Assuming z=1 for all peptides.')
        #print('Your D uptake values will be severely underestimated if actually z>1.\n')
        data_df['z'] = 1

    #Handle Renum
    #if renumbering_dict and len(renumbering_dict) >= 1:
    #    data_source_forunique1 = data_df.unstack(level=['Protein'])
    #    unique_proteins = data_source_forunique1.columns.get_level_values('Protein').unique().tolist()
    #    renumbering_dict = {prot: value for prot, value in renumbering_dict.items() if prot in unique_proteins}
    
    #Handle Mut Dict
        
    
    # Renumber residues, if desired:
    #renumbering_dict = None
    if renumbering_dict != None:
        if len(renumbering_dict)>0:
            #print('Renumbering residues of protein IDs specified in renumbering_dict.\n')
            for prot, shift in renumbering_dict.items():
                shift = pd.to_numeric(shift, errors='coerce')
                data_df.loc[data_df['Protein'] == prot, ['Start', 'End']] -= shift
                #data_df.loc[ data_df['Protein']==prot, ['Start','End'] ] = data_df.loc[ data_df['Protein']==prot, ['Start','End'] ] - shift

    def add_missing_peptides(data_df, wt_id, mut_id):
        #print("Into add_missing_peptides")
        wt_peptides = data_df[data_df["Protein"] == wt_id]
        mut_peptides = data_df[data_df["Protein"] == mut_id]
        #print("Just before defining missing peptides"))
        missing_peptides = wt_peptides.merge(
            mut_peptides,
            on=["Start", "End", "State"],
            how="left",
            indicator=True#,
            #suffixes=("", "_mut")
        ).query('_merge == "left_only"').drop(columns=["_merge"])
        if missing_peptides.empty:
            #print("No missing peptides found to add to the mutant protein.")
            return data_df
        new_entries = missing_peptides.copy()
        new_entries["Protein"] = mut_id 
        new_entries["State"] = new_entries["State"].apply(lambda x: x if x.endswith("-m") else f"{x}-m")
        new_entries["Original"] = 0
        columns_to_update = [
            "Protein", "Sequence", "MaxUptake", "MHP", "Exposure", "File", "z", "RT", "Inten", "Center", "Original"
        ]
        for col in columns_to_update:
            if f"{col}_x" in new_entries.columns:
                new_entries[col] = new_entries[f"{col}_x"]
                new_entries = new_entries.drop(columns=[f"{col}_x"])
        # Append the new entries to the original DataFrame
        updated_data_df = pd.concat([data_df, new_entries], ignore_index=True)
        #print(f"Added {len(new_entries)} missing peptides to mutant protein {mut_id}.")
        deduplicated_df = updated_data_df.drop_duplicates(subset=["Protein", "State", "Exposure", "z", "Inten", "Center"], keep="first")
        deduplicated_df = deduplicated_df[["Protein", "Start", "End", "Sequence", "Modification", "Fragment", "MaxUptake", "MHP", "State", "Exposure", "File", "z", "RT", "Inten", "Center", "Original"]]
        #print("Dropped duplicate rows based on 'Protein', 'State', 'z', and 'Center'.")
        return deduplicated_df

    if mutation_dict != None and mut_id_dict != None:
        if len(mutation_dict)>0 and len(mut_id_dict)==0:
            print("A list of residue mutations is provided in mutation_dict, but dictionary mut_id_dict is empty, thus it is impossible to know which protein IDs should be altered. No mutation will be applied.\n")
    
    def validate_peptide_mutant_correlation(data_df, mut_id, keyMutProt):
        # Filter for peptides where the Protein ID matches the mutant ID and correlation key is present
        filtered_df = data_df[data_df["Protein"] == mut_id]
        if filtered_df.empty:
            print(f"No peptides found for mutant protein ID '{mut_id}' with key '{keyMutProt}'.")
        return filtered_df
        
    addedMutState = 0 #always set to 0 first (therefore does not accidentally trigger some mutant handling steps later if mutants are not handled
    # Inside the main mutation processing code
    # mut_id_dict = None
    if mut_id_dict is not None:
        if len(mut_id_dict) > 0:
            data_df['Original'] = 1
            for mut_id, wt_id in mut_id_dict.items():
                # Ensure peptides are correlated to the mutant protein ID
                correlated_peptides = validate_peptide_mutant_correlation(data_df, mut_id, mut_id)
                if correlated_peptides.empty:
                    continue  # Skip processing if no correlation is found
                
                mut_indices = correlated_peptides.index
                print(f"Processing {len(mut_indices)} peptides for mutant protein '{mut_id}' and wild type '{wt_id}'.")
                
                for idx in mut_indices:
                    mut_start = data_df.loc[idx, 'Start']
                    mut_end = data_df.loc[idx, 'End']
                    mut_state = data_df.loc[idx, 'State']
                    mut_seq = data_df.loc[idx, 'Sequence']
                    
                    # Check if the sequence already exists in WT
                    seq_exist = data_df.loc[
                        (data_df['Protein'] == wt_id) &
                        (data_df['Start'] == mut_start) &
                        (data_df['End'] == mut_end) &
                        (data_df['State'] == mut_state) &
                        data_df['Original']
                    ].shape[0]
                    
                    if seq_exist:
                        print(f"Conflict detected for row {idx}: WT protein already contains a peptide with identical 'State' ({mut_state}). Adjusting mutant 'State' to ensure compatibility.")
                        data_df.loc[idx, 'State'] = f"{mut_state}-m"
                        mut_state = data_df.loc[idx, 'State']  # Update mut_state to reflect the new value
                        addedMutState = 1
                    else:
                        addedMutState = 0                
                    # Apply mutations if mutation_dict is provided
                    if mutation_dict:
                        for res, aa in mutation_dict.items():
                            if mut_start <= res <= mut_end:
                                if mut_seq[res - mut_start] == aa:
                                    print(f"Residue #{res} already matches {aa}. Skipping mutation.")
                                else:
                                    wt_seq = mut_seq[:res - mut_start] + aa + mut_seq[res - mut_start + 1:]
                    else:
                        # Fallback to WT sequence
                        wt_seq = data_df.loc[
                            (data_df['Protein'] == wt_id) & 
                            (data_df['Start'] == mut_start) & 
                            (data_df['End'] == mut_end)
                        ].iloc[0]['Sequence']
                    
                    # Update DataFrame with mutated information
                    data_df.loc[idx, 'Protein'] = wt_id
                    data_df.loc[idx, 'Sequence'] = wt_seq
                    data_df.loc[idx, 'Original'] = 0
                
                # Add missing peptides from WT to mutant state
                if addedMutState == 1:
                    data_df = add_missing_peptides(data_df, wt_id, mut_id) 
                
                ##### SAVE DATAFRAME TO A CSV IN AZURE #####
                #data_avg_reset = data_df.reset_index()
                # Save the reset DataFrame to a CSV
                #data_avg_reset.to_csv("data_avg_sample.csv", index=False)
                #save_to_azure_blob(data_avg_reset, CONTAINER_NAME, 'data_avg.csv', AZURE_STORAGE_CONNECTION_STRING)
                #print('')
    
    # Start d_uptake calculations.
    #Decharge:
    if file_typeDoH != 1:  # if file type is dynamx (decharge handled earlier for HDExaminer)
        data_df['Center'] = (data_df['Center'] - mp) * data_df['z']
    
    # Print or log a message about handling multiple charge states
    #if data_df_charge['z'].max() > 1:
    #    data_df = data_df_charge
    #    print("Multiple charge states detected for some peptides. Averaged across charge states.")
    
    #print('Specify State Prot (Before Fixing df States):', SpecifyStateProt)
    if SpecifyStateProt == 1:
        #print('editing states for SpecifyStateProt')
        # Reset index to modify 'State'
        data_df = data_df.reset_index()
        # Modify the 'State' column to include the Protein name
        data_df['State'] = data_df['State'] + data_df['Protein']
        # Print updated State values
        #print("Updated State values:")
        #print(data_df[['Protein', 'State']].head())  # Print first few rows
    
    # Compute average and std.dev. values for replicates at each time point
    data_avg = data_df.groupby(['Protein', 'Start', 'End', 'Sequence', 'State', 'Exposure']).agg({'Center': ['mean', 'std','count']})

    # Unstack the DataFrame to make 'Exposure' a separate level in the columns
    data_avg = data_avg.unstack('Exposure')
    
    #data_avg_reset = data_avg.reset_index()
    # Save the reset DataFrame to a CSV
    #data_avg_reset.to_csv("data_avg_sample.csv", index=False)
    #save_to_azure_blob(data_avg_reset, CONTAINER_NAME, 'data_avg.csv', AZURE_STORAGE_CONNECTION_STRING)
    
    data_avg.loc[:,('Center','mean',slice(None))] = data_avg.loc[:,('Center','mean',slice(None))].sub( data_avg['Center','mean',0], axis=0 )
    data_avg.loc[:,('Center','std',slice(None))]  = data_avg.loc[:,('Center','std',slice(None))].add( data_avg['Center','std',0].fillna(0), axis=0 )
    
    # Ensure that operations align with the multi-index structure
    # Subtract the time=0 column for each state from other time points to calculate D uptake
    #time_zero = data_avg.loc[:, ('Center', 'mean', 0)]
    #for level in ['mean', 'std']:
    #    if level == 'mean':
    #        # Subtract time=0 mean from all means
    #        data_avg.loc[:, ('Center', level, slice(None))] = (
    #            data_avg.loc[:, ('Center', level, slice(None))].sub(time_zero, axis=0)
    #        )
    #    elif level == 'std':
    #        # Add time=0 std to all std values
    #        data_avg.loc[:, ('Center', level, slice(None))] = (
    #            data_avg.loc[:, ('Center', level, slice(None))].add(time_zero.fillna(0), axis=0)
    #        )
    
    # Adjust column names for clarity
    data_avg.columns = data_avg.columns.set_levels(['d_uptake'], level=0)
    data_avg.columns = data_avg.columns.rename('Parameter', level=1)

    # Drop the no longer needed time=0 columns
    data_avg = data_avg.drop(columns=0, level='Exposure')

    percentUseRFU = float(request.form.get('UsePercentRFU', 1))
    if percentUseRFU != 100 or absoluteUptakeValues == 1:
        percentUseRFU = 1
    
    # Perform relative uptake calculation if enabled
    if relativeUptakeCalc == 1 and expType != 3:
        segment_lengths = (data_avg.index.get_level_values('End') - data_avg.index.get_level_values('Start') + 1)
        NumProline = data_avg.index.get_level_values('Sequence').str.count(r'[Pp]')
        MaxUptakesForRFU = (segment_lengths - (NumProline + 1))/percentUseRFU
        MaxUptakesForRFU = MaxUptakesForRFU.to_numpy()[:, None]
        if absoluteUptakeValues == 1:
            MaxUptakesForRFU = 1
        for level in ['mean', 'std']:
            data_avg.loc[:, ('d_uptake', level, slice(None))] /= MaxUptakesForRFU
    if expType == 3:
        segment_lengths = (data_avg.index.get_level_values('End') - data_avg.index.get_level_values('Start') + 1)
        NumProline = data_avg.index.get_level_values('Sequence').str.count(r'[Pp]')
        MaxUptakesForRFU = (segment_lengths - (NumProline + 1))/percentUseRFU
        MaxUptakesForRFU = MaxUptakesForRFU.to_numpy()[:, None]
        if absoluteUptakeValues == 1:
            MaxUptakesForRFU = 1
        for level in ['mean', 'std']:
            data_avg.loc[:, ('d_uptake', level, slice(None))] /= MaxUptakesForRFU
    #data_avg = data_avg.stack('State')

    if PulseLabellingNow == 2:
        #print("Index levels (1):", data_avg.index.names)
        #data_source_mixed = data_avg
        #data_source_mixed = data_source_mixed.stack('Exposure')
        #data_source_mixed = data_source_mixed.unstack(level=['State', 'Exposure'])
        #print("Index levels (2):", data_source_mixed.index.names)
        keep_list = data_avg.index.get_level_values('State') == state_of_interest
        #keep_list = state_of_interest
        #data_avg = data_avg.loc[keep_list,:]

    data_avg=data_avg.stack('Exposure').unstack(level=['State','Exposure'])
    
    def expand_peptides_to_residues(data_avg):
        """
        Expands peptide data in data_avg to residue-level data while maintaining the original structure
        and without modifying the original DataFrame. Replaces the Start and End columns with the residue number.
        """
        # Create a copy of data_avg with reset index for safe processing
        data_avg_reset = data_avg.reset_index()
        #save_to_azure_blob(data_avg_reset, CONTAINER_NAME, 'data_avg.csv', AZURE_STORAGE_CONNECTION_STRING)
        # Extract the column names
        columns = data_avg_reset.columns
        # Check if every sequence is present for every state value
        # Group by 'Sequence' and 'State' and count the occurrences
        state_count = data_avg_reset.groupby(['Sequence', 'State']).size().unstack(fill_value=0)
        # Find sequences that are missing any states by checking for any missing values (0)
        sequences_to_drop = state_count[state_count.min(axis=1) == 0].index
        # Drop the rows corresponding to those sequences
        data_avg_reset = data_avg_reset[~data_avg_reset['Sequence'].isin(sequences_to_drop)]
        # Drop the 'Sequence' column to avoid duplicates when creating residue rows
        #data_avg_reset = data_avg_reset.drop(columns=['Sequence'])
        # Prepare a list to store expanded rows
        expanded_rows = []
        # Iterate over each row in data_avg_reset
        for _, row in data_avg_reset.iterrows():
            start = int(row['Start'].iloc[0])
            end = int(row['End'].iloc[0])   
            # Generate residues for the given range
            for residue_pos in range(start, end + 1):
                # Create a copy of the row and set Start and End to the residue position
                new_row = row.copy()
                new_row['Start'] = residue_pos
                new_row['End'] = residue_pos  # Since Start == End for each residue
                #new_row['Sequence'] = residue_pos #Replace these values
                expanded_rows.append(new_row)
        # Create a new DataFrame with the expanded rows
        averaged_residue_df = pd.DataFrame(expanded_rows)
        # Reorganize columns to maintain the original structure, without the Residue column
        averaged_residue_df = averaged_residue_df[columns]
         # Create a multi-index similar to data_avg for averaging residue-level data
        averaged_residue_df = averaged_residue_df.set_index(['Protein', 'Start', 'End', 'Sequence', 'State', 'Exposure'])
        # Unstack the DataFrame to make 'Exposure' a separate level in the columns
        averaged_residue_df = averaged_residue_df.unstack('Exposure')
        # Adjust column names for clarity
        averaged_residue_df.columns = averaged_residue_df.columns.set_levels(['d_uptake'], level=0)
        averaged_residue_df.columns = averaged_residue_df.columns.rename('Parameter', level=1)
        # Now, average the d_uptake values within each group by 'Start'
        averaged_residue_df = averaged_residue_df.groupby(['Protein', 'Start', 'End', 'State']).mean()
        # Stack the DataFrame to restore the 'Exposure' level
        averaged_residue_df = averaged_residue_df.stack('Exposure', future_stack=True)
        
        return averaged_residue_df

    #PerResidueMap = float(request.form.get('PerResidueMap',0))
    PerResidueMap = 0
    #if uploaded_settings != 1:
        #if PerResidueMap == 1:
        #    averaged_residue_df = expand_peptides_to_residues(data_avg)
            #print(averaged_residue_df)
        #    averaged_residue_df = averaged_residue_df.unstack(level=['State','Exposure'])
        #    output_bitmap_dpi = 50
        #    data_avg = averaged_residue_df
        #else:
    
    #if PulseLabellingNow != 2:
    #    data_avg = data_avg.unstack(level=['State','Exposure'])
    #else:
    #    data_avg = data_avg.unstack(level=['Exposure'])
    
    # Reorganize data for easier column selection
    #print('Data Average:', data_avg)
    #print('Data Average Columns:', data_avg.columns)
    
    # Make a list of states for D uptake difference calculation:
    #data_source = averaged_residue_df if PerResidueMap == 1 else data_avg

    #data_avg_reset = data_avg.reset_index()
    # Save the reset DataFrame to a CSV
    #data_avg_reset.to_csv("data_avg_sample2.csv", index=False)
    #save_to_azure_blob(data_avg_reset, CONTAINER_NAME, 'data_avg2.csv', AZURE_STORAGE_CONNECTION_STRING)
    
    data_source = data_avg

    #data_avg_reset = data_source.reset_index()
    # Save the reset DataFrame to a CSV
    #data_avg_reset.to_csv("data_source_sample3.csv", index=False)
    #save_to_azure_blob(data_avg_reset, CONTAINER_NAME, 'data_source3.csv', AZURE_STORAGE_CONNECTION_STRING)
    
    data_source_forunique = data_source.unstack(level=['Protein', 'Start', 'End'])
    
    all_states = data_source.columns.get_level_values('State').unique()
    #print('all_states:', all_states)
    unique_proteins = data_source_forunique.columns.get_level_values('Protein').unique().tolist()
    #print("Unique Proteins: ", unique_proteins)
    unique_states = all_states.tolist()
    #print("Unique States: ", unique_states)
    unique_exposures = data_source_forunique.columns.get_level_values('Exposure').unique().tolist()
    #print("Unique Exposures: ", unique_exposures)
    protein_values = data_source_forunique.columns.get_level_values('Protein')
    start_values = data_source_forunique.columns.get_level_values('Start')
    end_values = data_source_forunique.columns.get_level_values('End')
    combinations = list(zip(protein_values, start_values, end_values))
    unique_peptides = list(set(combinations))
    #print("Unique Peptides: ", unique_peptides)
    
    #Handle Pulse labelling variables if not provided
    if expType == 2:
        if ref_time == 2020 or ref_time == '':
            ref_time = next((x for x in unique_exposures if x != 0), 0.5)
        if state_of_interest == 'NotProvided':
            state_of_interest = unique_states[0]
        if len(new_times) <= 0:
            new_times = [x for x in unique_exposures if x != 0 and x != ref_time]
        #print('New_times:', new_times)
        valid_exposures = {0, ref_time} | set(new_times)
        exposure_levels = data_source.columns.get_level_values('Exposure') 
        mask = exposure_levels.isin(valid_exposures)
        data_source = data_source.loc[:, mask] 
        unique_exposures = data_source.columns.get_level_values('Exposure').unique().tolist()
        #print("Valid Exposures: ", valid_exposures)
        #print("Unique Exposures (2): ", unique_exposures)
        set_unique_exposures = set(unique_exposures)
        set_new_times = set(valid_exposures)
        difference = set_unique_exposures - set_new_times
        result_list = list(difference)
        #data_source.stack('Exposure')
        #data_avg_reset = data_source.reset_index()
        # Save the reset DataFrame to a CSV
        #data_avg_reset.to_csv("data_source_sample4.csv", index=False)
        #save_to_azure_blob(data_avg_reset, CONTAINER_NAME, 'data_source4.csv', AZURE_STORAGE_CONNECTION_STRING)
        data_source = data_source.drop(columns=result_list, level='Exposure') #updated
        #data_avg_reset = data_source.reset_index()
        # Save the reset DataFrame to a CSV
        #data_avg_reset.to_csv("data_source_sample5.csv", index=False)
        #save_to_azure_blob(data_avg_reset, CONTAINER_NAME, 'data_source5.csv', AZURE_STORAGE_CONNECTION_STRING)
        #data_source.unstack(level = 'Exposure')
        comp_times = new_times
            
    #handle drop_prot
    if drop_prot and len(drop_prot) >= 1:
        drop_prot = [prot for prot in drop_prot if prot in unique_proteins]
        print('Drop Prot list after filtering:', drop_prot)
    #handle drop_times
    if drop_times and len(drop_times) >= 1:
        drop_times = [time for time in drop_times if time in unique_exposures]
    #handle drop_pept
    if drop_pept and len(drop_pept) >= 1:
        drop_pept = [entry for entry in drop_pept if entry in unique_peptides]
    #handle state_list
    if expType != 3:
        if len(all_states) == 1 and expType != 2:
            #print('There is only 1 state in the input csv file: \'%s\'' % all_states[0])
            #print('Need 2 or more states to calculate D uptake difference.\n')
            raise StopExecution
        if state1_list != None and state2_list != None:
            #if len(state1_list) != len(state2_list):
                #print('state1_list and state2_list (the lists of protein states for difference calculation) have unequal lengths.')
                #print('If this is not intentional, there may be some D uptake difference calculations missing from the output.\n')
            filtered_state1 = []
            filtered_state2 = []
            for s1, s2 in zip(state1_list, state2_list):
                if s1 in unique_states and s2 in unique_states:
                    filtered_state1.append(s1)
                    filtered_state2.append(s2)
            #print('Filtered State 1:', filtered_state1)
            #print('Filtered State 2:', filtered_state2)
            if len(filtered_state1) > 0 and len(filtered_state2) > 0:
            #state_list = list(zip(state1_list, state2_list))
                state_list = list(zip(filtered_state1, filtered_state2))
                print('State List:', state_list)
        if not state_list:
            if addedMutState == 0:
                state_list = list(itertools.combinations(all_states, 2))
            else:
                mutant_states = [state for state in all_states if state.endswith('-m')] #artificially created Mutant states
                non_mutant_states = [state for state in all_states if not state.endswith('-m')]
                nm_state_list = list(itertools.combinations(non_mutant_states, 2))
                m_state_list = list(itertools.combinations(mutant_states, 2))
                state_list = nm_state_list + m_state_list
                # print(f"Non-mutant states: {non_mutant_states}")
                # print(f"Mutant states: {mutant_states}")
                # print(f"Combined state list: {combined_state_list}")
            print('A list of protein states for D uptake difference calculation was not provided.')
            print('There are a total of %d states in the input csv file.' % len(all_states))
            print('Will calculate all-against-all D uptake differences, total of %d comparisons.\n' % len(state_list))

    if PulseLabellingNow == 2:
        #data_avg_reset = data_source.reset_index()
        # Save the reset DataFrame to a CSV
        #data_avg_reset.to_csv("data_source_sample6-1.csv", index=False)
        #save_to_azure_blob(data_avg_reset, CONTAINER_NAME, 'data_source6-1.csv', AZURE_STORAGE_CONNECTION_STRING)
        print(keep_list)
        #data_source = data_source.loc[keep_list,:]
        
        #data_avg_reset = data_source.reset_index()
        # Save the reset DataFrame to a CSV
        #data_avg_reset.to_csv("data_source_sample6-2.csv", index=False)
        #save_to_azure_blob(data_avg_reset, CONTAINER_NAME, 'data_source6-2.csv', AZURE_STORAGE_CONNECTION_STRING)
    
    # Calculate D uptake differences between states. Append data to data_avg:
    diffs = pd.DataFrame()
    #print("data source called columns:", data_source.columns)
    #group_levels = ['Protein', 'Start', 'State', 'Exposure'] if PerResidueMap == 1 else data_avg.columns.names

    #print('data_avg df Columns:', data_avg.columns)
    #print('data_avg df Column Names:', data_avg.columns.names)
    #print('data_avg df Column Levels:', data_avg.columns.levels)

    #print('PerResidueMap Value:', PerResidueMap)
    
    #if PerResidueMap == 1:
        #print('Averaged df Columns:', averaged_residue_df.columns)
        #print('Averaged df Column Names:', averaged_residue_df.columns.names)
        #print('Averaged df Column Levels:', averaged_residue_df.columns.levels)
    if expType != 3:  
        if PulseLabellingNow != 2:
            for two_states in state_list:
                state1 = str(two_states[0])
                state2 = str(two_states[1])
                if state1 == state2:
                    print('Found a request to calculate difference between a state \'%s\' and itself.' % state1)
                    print('There may be a typo, or the orders of states might be mixed up, in the \'state1_list\' or \'state2_list\'.\n')
                mean1 = data_source.loc[:, ('d_uptake', 'mean', state1, slice(None))]
                mean2 = data_source.loc[:, ('d_uptake', 'mean', state2, slice(None))]
                #print('mean 1 shape', mean1.shape)
                #print('mean 2 shape', mean2.shape)
                std1 = data_source.loc[:, ('d_uptake', 'std', state1, slice(None))]
                std2 = data_source.loc[:, ('d_uptake', 'std', state2, slice(None))]
                count1 = data_source.loc[:, ('d_uptake', 'count', state1, slice(None))]
                count2 = data_source.loc[:, ('d_uptake', 'count', state2, slice(None))]
                # Calculate differences
                mean = data_source.loc[:,('d_uptake','mean',state1,slice(None))] - data_source.loc[:,('d_uptake','mean',state2,slice(None))].values
                std  = data_source.loc[:,('d_uptake','std',state1,slice(None))]  + data_source.loc[:,('d_uptake','std',state2,slice(None))].values
                pval = ttest_ind_from_stats(data_source.loc[:,('d_uptake','mean',state1,slice(None))], data_source.loc[:,('d_uptake','std',state1,slice(None))], data_source.loc[:,('d_uptake','count',state1,slice(None))], data_source.loc[:,('d_uptake','mean',state2,slice(None))].values, data_source.loc[:,('d_uptake','std',state2,slice(None))].values, data_source.loc[:,('d_uptake','count',state2,slice(None))].values, equal_var=False)[1]
                pval = pval.rename(columns={'count': 'pval'}, level='Parameter')
                new_state = str(state1) + ' - ' + str(state2)
                d = {state1 : new_state}
                mean = mean.rename(columns=d, level=2).round(2)  # Round mean values
                std = std.rename(columns=d, level=2)   
                pval = pval.rename(columns=d, level=2)
                #print("diffs columns:", diffs.columns)
                #print("mean columns:", mean.columns)
                #print("std columns:", std.columns)
                #print("pval columns:", pval.columns)
                #print("diffs index:", diffs.index)
                #print("mean index:", mean.index)
                #print("std index:", std.index)
                #print("pval index:", pval.index)
                #print("diffs column levels:", diffs.columns.names)
                #print("mean column levels:", mean.columns.names)
                #print("std column levels:", std.columns.names)
                #print("pval column levels:", pval.columns.names)
                diffs = pd.concat([diffs, mean, std, pval], axis=1)
            diffs.columns = diffs.columns.set_levels(['Delta_d_uptake'], level=0)
        else:
            names = []
            #data_source.stack('Exposure')
            #data_source.unstack('Exposure')
            #print("Column levels:", data_source.columns.names)
            #print("Example column names:", data_source.columns[:5].tolist())
            #data_source.columns = pd.MultiIndex.from_tuples(
            #    [(a, b, c, float(d)) for a, b, c, d in data_source.columns],
            #    names=['None', 'Parameter', 'State', 'Exposure']  # Set column names explicitly
            #)
            #data_avg_reset = data_source.reset_index()
            #data_avg_reset.to_csv("data_source_sample7.csv", index=False)
            #save_to_azure_blob(data_avg_reset, CONTAINER_NAME, 'data_source7.csv', AZURE_STORAGE_CONNECTION_STRING)
            #print("Column levels (2):", data_source.columns.names)
            #print("Example column names (2):", data_source.columns[:5].tolist())
            for time2 in comp_times:
                #data_source.stack('Exposure')
                available_times = data_source.columns.get_level_values('Exposure').unique()
                #print("Available time points:", available_times)
                #time2 = str(time2)
                #print("Available times types:", [type(t) for t in available_times])
                time2 = min(available_times, key=lambda x: abs(x - float(time2)))
                ref_time = min(available_times, key=lambda x: abs(x - float(ref_time)))
                #expected_type = type(available_times[0])
                #expected_type = available_times.dtype  # Pandas stores the dtype of Index/Series
                time2 = float(time2)  # Convert explicitly
                ref_time = float(ref_time)
                #print("Requested time point:", time2)
                #print("Type of time2 after conversion:", type(time2))
                #print("Type of ref_time after conversion:", type(ref_time))
                #time2 = float(time2)  # Ensure it matches the actual type
                #ref_time = float(ref_time)
                mean = data_source.loc[:,('d_uptake','mean', state_of_interest, time2)] - data_source.loc[:,('d_uptake','mean', state_of_interest, ref_time)].values #UPDATED M5
                pval_temp = ttest_ind_from_stats( data_source.loc[:,('d_uptake','mean', state_of_interest,ref_time)], data_source.loc[:,('d_uptake','std', state_of_interest,ref_time)].map(np.sqrt), data_source.loc[:,('d_uptake','count', state_of_interest,ref_time)], data_source.loc[:,('d_uptake','mean', state_of_interest,time2)].values, data_source.loc[:,('d_uptake','std', state_of_interest,time2)].map(np.sqrt).values, data_source.loc[:,('d_uptake','count', state_of_interest,time2)].values, equal_var=False)[1]
                #pval_temp = ttest_ind_from_stats( data_avg.loc[:,('d_uptake','mean',state1,slice(None))], data_avg.loc[:,('d_uptake','var',state1,slice(None))].map(np.sqrt), data_avg.loc[:,('d_uptake','count',state1,slice(None))], data_avg.loc[:,('d_uptake','mean',state2,slice(None))].values, data_avg.loc[:,('d_uptake','var',state2,slice(None))].map(np.sqrt).values, data_avg.loc[:,('d_uptake','count',state2,slice(None))].values, equal_var=False)[1]
                pval = pd.DataFrame(pval_temp,index=data_source.loc[:,('d_uptake','mean', state_of_interest,ref_time)].index)
                diffs = pd.concat([diffs, mean, pval], axis=1)
                names = names + [mean.name, ('d_uptake', 'pval', state_of_interest, time2)]
            ####
            # Delete if not needed inside for loop
            ####
            #     pval_temp = ttest_ind_from_stats( data_avg.loc[:,('d_uptake','mean',state1,slice(None))], data_avg.loc[:,('d_uptake','var',state1,slice(None))].map(np.sqrt), data_avg.loc[:,('d_uptake','count',state1,slice(None))], data_avg.loc[:,('d_uptake','mean',state2,slice(None))].values, data_avg.loc[:,('d_uptake','var',state2,slice(None))].map(np.sqrt).values, data_avg.loc[:,('d_uptake','count',state2,slice(None))].values, equal_var=False)[1]
            #     pval = pd.DataFrame(pval_temp,columns=data_avg.loc[:,('d_uptake','count',state1,slice(None))].columns,index=data_avg.loc[:,('d_uptake','count',state1,slice(None))].index)
            #     pval = pval.rename(columns={'count': 'pval'}, level='Parameter')
            ####
            diffs.columns = pd.MultiIndex.from_tuples(names, names=["", "Parameter", "State", "Exposure"])
            diffs = diffs.rename(columns={'d_uptake': 'Delta_d_uptake'}, level=0)
            idx = data_source.index
            data_source=pd.concat([data_source, diffs], axis=1)
            data_source.index = idx

    if expType != 3:
        data_source = pd.concat([data_source, diffs], axis=1)

    try:
        data_source.to_csv(output_csv_file)
        print("Processed data was written into a csv file '%s'\n" % output_csv_file)
    except IOError:
        curr_time = time.localtime()
        time_string = "%s.%s.%s_%s.%s.%s" % (curr_time[:6])
        new_output_csv_file = output_csv_file + "_" + time_string + "_.csv"
        data_source.to_csv(new_output_csv_file)
        print("Could not write csv output into file '%s'" % output_csv_file)
        print("On Windows this often happens when file already exists and is open.")
        print("Your csv output was written into a file with a new file name containing an appended unique date and time string.")
        print("Date and time string is in the format of: year.month.day_hour.minute.second")
        print("New file name is '%s'\n" % new_output_csv_file)

    # Drop time points unwanted by the user, individual peptides unwanted by the user, and entire proteins unwanted by the user:
    if drop_times != None:
        data_source = data_source.drop(columns=drop_times, level='Exposure')
    if drop_pept != None:
        data_source = data_source.drop(index=drop_pept)
    if drop_prot != None:    
        data_source = data_source.drop(index=drop_prot)
        
    if expType != 3:
        all_deltas = data_source.loc[:,('Delta_d_uptake','mean',slice(None),slice(None))]
    else:
        all_deltas = data_source.loc[:,('d_uptake','mean',slice(None),slice(None))]
    #print(all_deltas)
    max_delta_pos = all_deltas.max().max()  # axis=None does not work in older versions of pandas
    max_delta_neg = all_deltas.min().min()
    if funkybound == 1 and globalmax_delta != 1:
        #print('max delta positive', max_delta_pos)
        max_delta = max_delta_pos
    elif funkybound == 2 and globalmax_delta != 1:
        #print("max delta negative", max_delta_neg)
        max_delta = abs(max_delta_neg)
    else:
        max_delta = max([abs(max_delta_pos), abs(max_delta_neg)])   
    if custb == 2:
        max_range = round(max_delta, 1)
        if max_range <= 3 and relativeUptakeCalc != 1:
            max_range = 3
    if max_range == 0:
        max_range = 1
    if custb == 2:
        if alt_col == 1:
            c1 = request.form.get('negcolalt',bluec)
            c3 = request.form.get('poscolalt',redc)
    
    # A function to create bounds list for heatmap color bins. Output will be passed into matplotlib.colors.BoundaryNorm()
    def make_bounds(tot_colors, max_delta=max_delta):
        #print('tot_colors being passed:', tot_colors)
        if uploaded_settings != 1:
            custb = float(request.form['option'])
            zerobound = float(request.form.get('zerobound', 0))
            funkybound = float(request.form.get('change_bounds_abs',3))
            if custb ==3:
                max_range = float(request.form['max_range'])
                max_range = round(max_range, 1)
        if uploaded_settings == 1:
            if 'custb' in variables:
                try:
                    custb = variables.get('custb')
                except ValueError:
                    print('custb not added/updated')
                    pass
                else:
                    print('custb = ', custb)
            if 'funkybound' in variables:
                try:
                    funkybound = variables.get('funkybound')
                except ValueError:
                    print('funkybound not added/updated')
                    pass
                else:
                    print('funkybound = ', funkybound)
            if 'zerobound' in variables:
                try:
                    zerobound = variables.get('zerobound')
                except ValueError:
                    print('zerobound not added/updated')
                    pass
                else:
                    print('zerobound = ', zerobound)
        if custb ==2:
            max_range = round(max_delta, 1)
        num_shades = np.ceil(tot_colors / 2).astype(int) # larger by 1 for odd tot_colors  
        if (funkybound == 1 or funkybound == 2) and (custb == 2 or custb == 3): 
            num_shades = np.ceil(tot_colors).astype(int)
        num_shades = round(num_shades)
        if max_range:
            incr = (max_range - 0.5) / num_shades
        else:
            incr = np.ceil(10 * max_delta / num_shades)/10
            max_range = incr * num_shades
        if tot_colors % 2 == 0 and (((custb != 3 and custb != 3.0) and (custb != 2 and custb != 2.0)) and ((funkybound != 1 and funkybound != 1.0) and (funkybound != 2 and funkybound != 2.0))):
            #print("In tot_colors % 2")
            #print('custb value:', custb)
            #print('funkybound value:', funkybound)
            bounds = np.round(np.linspace(-max_range, max_range, num = tot_colors+1), 2)
            if zerobound == 1:
                bounds = np.round(np.linspace(-max_range, max_range, num = tot_colors), 2)
                bounds = np.append(bounds, 0)
        else:
            if custb == 2 and alt_col == 0: 
                #print('num shades:',num_shades)
                if max_range <= 10:
                    half = np.round(np.linspace(0.5, max_range, num=num_shades), 2)
                else:
                    half = np.round(np.linspace(0, max_range, num=num_shades+1), 2)
                    half = half[1:]
                #print('half:', half)
            elif custb == 3 and (funkybound == 1 or funkybound == 2):
                #print('num shades:',num_shades)
                if max_range <= 10:
                    half = np.round(np.linspace(0.5, max_range, num=num_shades), 2)
                else:
                    half = np.round(np.linspace(0, max_range, num=num_shades+1), 2)
                    half = half[1:]
                #print('half:', half)
            elif custb == 2 and alt_col == 1:
                #if zerobound == 1:
                #    num_shades = num_shades - 1
                if max_range <= 10:
                    half = np.round(np.linspace(0.5, max_range, num=num_shades), 2)
                else:
                    half = np.round(np.linspace(0, max_range, num=num_shades+1), 2)
                    half = half[1:]
                #print('num shades:', num_shades, 'half:', half)
            else:    
                half = np.round(np.linspace(0.5, max_range, num = num_shades-1),2)
                #print('num shades - 1, half:', half)
            if zerobound == 0:
                if funkybound == 1:
                    bounds = np.copy(half) 
                elif funkybound == 2:
                    bounds = np.copy(-half)  
                else:
                    bounds = np.append(-half, half) 
            elif zerobound == 1:
                if funkybound == 1:
                    bounds = np.copy(half)  
                elif funkybound == 2:
                    bounds = np.copy(-half) 
                else:
                    bounds = np.append(-half, half)
                bounds = np.append(bounds, 0) 
            bounds = np.sort(bounds)
        if zerobound == 1 or funkybound == 1 or funkybound == 2:
            bounds = np.append(bounds, 0)
            bounds = list(set(bounds))
            bounds = np.sort(bounds)
        #print('Bound List:', bounds)
        #print('Funky Bound Num:', funkybound)
        return bounds

    print("Maximal difference in D uptake between protein states in the data is %.2f." % max_delta)
    if max_range != None:
        if max_range:
            print("The user has set max_range=%.2f; using this value for the color scale instead of above maximal difference." % max_range)
    print('')
    
    # Drop columns and rows that are completely empty in all_deltas:
    old_shape = all_deltas.shape
    all_deltas = all_deltas.dropna(axis=0, how='all') # rows
    all_deltas = all_deltas.dropna(axis=1, how='all') # columns
    if all_deltas.shape[0] < old_shape[0]:
        print('Before plotting, removed %d rows from D uptake difference table that are completely empty. This happens when\nsome peptides are present in only one of the protein states.\nYou may wish to inspect output csv file to make sure processed data looks reasonable.\n' % (old_shape[0] - all_deltas.shape[0]))
    if all_deltas.shape[1] < old_shape[1]:
        print('Before plotting, removed %d columns from D uptake difference table that are completely empty. This happens when\nsome time points are present in only one of the protein states.\nYou may wish to inspect output csv file to make sure processed data looks reasonable.\n' % (old_shape[1] - all_deltas.shape[1]))
    data_source.to_csv("data_source.csv",index=False)
    #save_to_azure_blob(data_source, CONTAINER_NAME, 'data_source.csv', AZURE_STORAGE_CONNECTION_STRING)
    if expType != 3:
        all_pvals           = data_source.loc[:,('Delta_d_uptake','pval',slice(None),slice(None))]        # after per forming .dropna() on all_deltas, potentially there are more columns in all_pvals than in all_deltas
    #else: 
    #    all_pvals           = data_source.loc[:,('d_uptake','pval',slice(None),slice(None))]
    #all_pvals.to_csv("all_pvals1.csv",index=False)
    #save_to_azure_blob(all_pvals, CONTAINER_NAME, 'all_pvals1.csv', AZURE_STORAGE_CONNECTION_STRING)
    if expType != 3:
        columns_of_interest = all_deltas.rename(columns={'mean': 'pval'}, level='Parameter').columns   # columns from all_deltas, renamed to look like from all_pvals
        all_pvals           = all_pvals.loc[:, columns_of_interest]                                    # now same columns are in all_pvals as in all_deltas
        #all_pvals.to_csv("all_pvals2.csv",index=False)
        #save_to_azure_blob(all_pvals, CONTAINER_NAME, 'all_pvals2.csv', AZURE_STORAGE_CONNECTION_STRING)
        all_pvals           = all_pvals.loc[all_deltas.index, :]                                       # same for rows
    
    # Make p-value vs. D uptake difference volcano plots: create a directory; then loop through all states and times to plot data
    all_states = all_deltas.columns.get_level_values('State').unique()

    #print('prot list')
    if PulseLabellingNow != 2:
        all_prot = all_deltas.groupby('Protein')
    else: 
        all_prot = all_deltas
    #print(all_prot)
    numprotinfile = len(all_prot)
    #print(numprotinfile)

    def is_valid_hex_color(value):
        hex_pattern = re.compile(r'^#[0-9A-Fa-f]{6}$')
        return bool(hex_pattern.match(value))

    #print("Shape of all_deltas after cleaning:", all_deltas.shape)
    #print("Shape of all_pvals after alignment:", all_pvals.shape)
    #print("NaN count in all_deltas:", all_deltas.isna().sum().sum())
    #print("NaN count in all_pvals:", all_pvals.isna().sum().sum())
    #print("Columns in all_deltas:", all_deltas.columns)
    #print("Columns in all_pvals:", all_pvals.columns)


    ##################################################################################
    ###########################      SCATTER PLOTTING      ###########################
    ##################################################################################


    if uploaded_settings != 1:
        p_threshold = request.form.get('pthresh', 0.05)
        altvolc = request.form.get('altvolc', 0)
        scattercolor = blackc
    if altvolc == 1:
        if uploaded_settings != 1:
            altvolccol = request.form.get('altvolccol', blackc)
        if is_valid_hex_color(altvolccol) == True:
            scattercolor = altvolccol
        else:
            scattercolor = blackc
    else:
        scattercolor = blackc
    if p_threshold == '':
        p_threshold = 0.05
    if p_threshold != 0.05 and p_threshold != '':
        p_threshold = float(p_threshold)
    if scatter_plot == 1 and expType != 3: #RIGHT NOW DOESN'T WORK W RFU BC USES Delta_d_uptake
        all_deltas_grouped = all_deltas.groupby('Protein')
        if scatter_dir:
            Path(scatter_dir).mkdir(parents=True, exist_ok=True)
        scatter_buffer = io.BytesIO()  # Create an in-memory zip file buffer
        with zipfile.ZipFile(scatter_buffer, 'a', zipfile.ZIP_DEFLATED, False) as zip_file:
            for prot, data_subset in all_deltas_grouped:
                for state in all_states:
                    all_times = all_deltas.loc[:, ('Delta_d_uptake', 'mean', state, slice(None))].columns.get_level_values('Exposure').unique()
                    for i in all_times:
                        #print("Checking data for scatter plot:")
                        #print(all_deltas.loc[:, ('Delta_d_uptake', 'mean', state, i)].head())
                        #print(all_pvals.loc[:, ('Delta_d_uptake', 'pval', state, i)].head())
                        if uploaded_settings != 1:
                            if editing_figtitle == 1:
                                title = str(request.form.get('alttitleword','Invalid Text Entry'))
                            else:
                                title = f'{state}, exposure: {i}'
                        fig = plt.figure(figsize=(12, 10))
                        #print(all_deltas.loc[:, ('Delta_d_uptake', 'mean', state, i)], -all_pvals.loc[:, ('Delta_d_uptake', 'pval', state, i)].map(np.log10))
                        plt.plot(all_deltas.loc[:, ('Delta_d_uptake', 'mean', state, i)], -all_pvals.loc[:, ('Delta_d_uptake', 'pval', state, i)].map(np.log10), markerfacecolor=scattercolor, markeredgecolor=scattercolor, marker='o', linestyle='None')
                        plt.axhline(y=-np.log10(p_threshold), color='k', linewidth=1)
                        plt.axvline(x=-Dd_threshold, color='k', linewidth=1)
                        plt.axvline(x=Dd_threshold, color='k', linewidth=1)
                        plt.title(title)
                        plt.ylabel(new_y_axis_title)
                        plt.xlabel(new_x_axis_title)
                        plt.xlim([-max_delta * 1.05, max_delta * 1.05])
                        plt.ylim(bottom=0)
                        plt.tick_params(axis='both', labelsize=font_size_ticklabel)
                        file_name = f"volcano_plot_{prot}_{state}_{str(i)}.{output_bitmap_format}"
                        plt.savefig(file_name, dpi=output_bitmap_dpi)
                        plt.close()
                        if scatter_dir:
                            file_path = f"{scatter_dir}/{file_name}"
                        else:
                            file_path = file_name
                        zip_file.write(file_name)  # Add the file to the zip archive
                        Path(file_name).unlink()  # Remove the file after adding it to the zip archive
        # Save the zip buffer to a file and then upload to Azure Blob Storage
        scatter_buffer.seek(0)
        blob_client = blob_service_client.get_blob_client(container=CONTAINER_NAME, blob='scatter_plots.zip')
        # Upload the in-memory zip file buffer to Azure Blob Storage
        blob_client.upload_blob(scatter_buffer.getvalue(), overwrite=True)


    ##################################################################################
    ##################################################################################
    ##################################################################################


    # Data filtration based on p-values:
    if significant_only and expType != 3:
        count_pvals = all_pvals.notna().sum()
        no_reps_bool = (count_pvals==0)
        no_reps_tot  = no_reps_bool.sum()
        all_pvals_bool = (all_pvals>p_threshold)
        #print("all_deltas shape:", all_deltas.shape)
        #print("all_pvals_bool shape:", all_pvals_bool.shape)
        #print("all_deltas columns unique:", all_deltas.columns.is_unique)
        #print("all_pvals_bool columns unique:", all_pvals_bool.columns.is_unique)
        #print("all_deltas index unique:", all_deltas.index.is_unique)
        #print("all_pvals_bool index unique:", all_pvals_bool.index.is_unique)
        #print("all_deltas columns:", all_deltas.columns.tolist())
        #print("all_pvals_bool columns:", all_pvals_bool.columns.tolist())
        #print("Duplicate columns in all_pvals_bool:", all_pvals_bool.columns[all_pvals_bool.columns.duplicated()])
        if PulseLabellingNow == 2:
            all_deltas = all_deltas.loc[:, ~all_deltas.columns.duplicated()]
            all_pvals_bool = all_pvals_bool.loc[:, ~all_pvals_bool.columns.duplicated()]
            all_pvals_bool = all_pvals_bool.reindex_like(all_deltas)#, fill_value=False)
            all_deltas = all_deltas.where(all_pvals_bool.astype(bool), 0)              
        else:
            all_deltas = all_deltas.mask(all_pvals_bool.values, 0)
        if no_reps_tot == len(count_pvals):   # all data is missing replicates
            print('None of %d protein states and exposure time combinations appear to contain replicated data. No P-value filtering will be applied.\n' % no_reps_tot)
        else:
            if no_reps_tot > 0:   # looks like some, but not all data has missing replicates
                print('States and expo time combinations do not have replicated data')
                #print('%d of %d protein states and exposure time combinations do not have replicated data. P-value filtering is not applied to this data.\n' % no_reps_tot, len(count_pvals))
            print('Set %d D uptake differences with larger than %.2e p-values to 0 before plotting the data.\n' % (all_pvals_bool.sum().sum(), p_threshold))

            
    if custom_colors != None: 
        if custom_colors:
            colors = custom_colors
            if custom_bounds != None:
                if custom_bounds:
                    bounds = custom_bounds
                    print('bounds and ccol')
            else:
                bounds = make_bounds(len(colors))
                print('no bounds found')
    elif custom_bounds != None: 
        if custom_bounds:
            print("The user has provided custom_bounds list and no custom_colors list. Using custom_bounds to calculate colors, and ingnoring max_range and num_shades variables.\n")
            b = np.array(custom_bounds)
            neg_shades = len(b[b<0])
            pos_shades = len(b[b>0])
            col_mid = []
            if len(b[b==0])==0:
                if neg_shades>0 and pos_shades>0:
                    col_mid = [ mc.to_hex(c2) ]
                    if zerobound == 0:
                        neg_shades = len(b[b<0])-1
                        pos_shades = len(b[b>0])-1
                    if zerobound == 1:
                        neg_shades = len(b[b<0])
                        pos_shades = len(b[b>0])
            if neg_shades > 0:
                cols1 = [colorFader(c1,c2,x/neg_shades) for x in range(neg_shades+1)]
            else:
                cols1 = [ mc.to_hex(c2) ]
            if pos_shades>0:
                cols2 = [colorFader(c2,c3,x/pos_shades) for x in range(pos_shades+1)]
            else:
                cols2 = [ mc.to_hex(c2) ]
            if zerobound == 0:
                colors = cols1[:-1] + col_mid + cols2[1:]
            if zerobound == 1:
                colors = cols1[:-1] + cols2[1:]
            bounds = custom_bounds
            #print(bounds)
            #print('no ccol')
    else: #no custom colours or bounds (auto or set range option)
        #c_missing = '#bdbdbd' # gray;
        num_shades = round(num_shades) # Convert num_shades to integer if it's a float
        #print('num_shades:', num_shades)
        #print('custb:', custb)
        #print('zerobound:', zerobound)
        if zerobound == 0 and custb != 3:
            cols1 = [colorFader(c1,c2,x/num_shades) for x in range(num_shades+1)]
            cols2 = [colorFader(c2,c3,x/num_shades) for x in range(num_shades+1)]
            if funkybound == 1: #only positive
                colors = cols2[1:] 
            elif funkybound == 2: #only negative
                colors = cols1
                print('not custb 3 and fb == 2, colors:', colors)
            else: #positive and negative
                colors = cols1 + cols2[1:]
        elif zerobound == 0 and custb == 3:
            if funkybound == 1 or funkybound == 2:
                cols1 = [colorFader(c1,c2,x/num_shades) for x in range(num_shades-1)]
                cols2 = [colorFader(c2,c3,x/num_shades) for x in range(num_shades)]
                if funkybound == 1:
                    colors = [c2] + cols2[1:]
                    print('custb 3 and fb 1, colors:', colors)
                if funkybound == 2:
                    colors = cols1 + [c2]
                    print('custb 3 and fb 2, colors:', colors)
            else:
                cols1 = [colorFader(c1,c2,x/num_shades) for x in range(num_shades+1)]
                cols2 = [colorFader(c2,c3,x/num_shades) for x in range(num_shades+1)]
                colors = cols1 + cols2[1:]
        if zerobound == 1:
            cols1 = [colorFader(c1, c2, (x/num_shades)) for x in range(num_shades+1)]
            cols2 = [colorFader(c2, c3, (x/num_shades)) for x in range(num_shades+1)]
            if funkybound == 1:
                colors = cols2[1:]
            elif funkybound == 2:
                colors = cols1[:-1]
                print('in 0 bound')
            else:
                colors = cols1[:-1] + cols2[1:]
        bounds = make_bounds(len(colors))
        if zerobound == 1:
            bounds = np.append(bounds, 0)
            bounds = list(set(bounds))
            bounds = np.sort(bounds)
        #print(bounds)
        #print('no ccol or bounds found')

    colormap = mc.ListedColormap(colors)
    #print(colors)
    #print(bounds)
    #print('right before calling')
    my_norm = mc.BoundaryNorm(bounds, ncolors=len(colors))

    if PulseLabellingNow == 2:
        all_deltas.index.set_names(['Protein','Start','End','Sequence'], inplace=True)
        #all_deltas_reset = all_deltas.reset_index()
        #all_deltas_reset.to_csv("all_deltas_full3750.csv", index=False)
        #save_to_azure_blob(all_deltas_reset, CONTAINER_NAME, 'all_deltas_full2.csv', AZURE_STORAGE_CONNECTION_STRING)
    
    all_deltas.reset_index(level=['Protein','Start','End'],inplace=True)
        
    # Open pdf output file:
    pdf=PdfPages(output_pdf_file) 
    try:
        pdf.attach_note('')
    except IOError:
        curr_time = time.localtime()
        time_string = "%s.%s.%s_%s.%s.%s" % (curr_time[:6])
        new_output_pdf_file = output_pdf_file + "_" + time_string + "_.pdf"
        pdf=PdfPages(new_output_pdf_file)
        print("Could not write pdf output into file '%s'" % output_pdf_file)
        print("On Windows this often happens when file already exists and is open.")
        print("Your pdf output was written into a file with a new file name containing an appended unique date and time string.")
        print("Date and time string is in the format of: year.month.day_hour.minute.second")
        print("New file name is '%s'\n" % new_output_pdf_file)
    plt.rcParams['font.size'] = font_size
    plt.rcParams['text.color'] = fontcolor
    plt.rcParams['axes.labelcolor'] = fontcolor
    
    #what_t_unit = 0
    if uploaded_settings != 1:
        all_in_min = float(request.form.get('all_in_min',0))
        if all_in_min == 1:
            what_t_unit = float(request.form.get('what_t_unit',0))

    if PulseLabellingNow == 2 and float(ref_time) >= 1:
        ref_time = round(float(ref_time), 1)
        ref_time = f'{ref_time}'.rstrip('0').rstrip('.')
    
    def make_time_tick_labels(time_list):
        labels = []
        if all_in_min == 0 or what_t_unit == 0:
            for hdx_time in time_list:
                if hdx_time < 1:
                    string = '%s s' % np.round(hdx_time*60).astype(int)
                elif hdx_time < 60:
                    mins = np.floor(hdx_time).astype(int)
                    sec = np.round(np.remainder(hdx_time, 1)*60).astype(int)
                    if sec > 0:
                        string = '%s m %s s' % (mins, sec)
                    else:
                        string = '%s m' % mins
                else:
                    hour = np.round(hdx_time/60).astype(int)
                    mins  = np.round(np.remainder(hdx_time, 60)).astype(int)
                    if mins > 0:
                        string = '%s h %s m' % (hour, mins)
                    else:
                        string = '%s h' % hour
                if PulseLabellingNow == 2:
                    string += f' - {ref_time} m'
                labels.append(string)
        else:
            if what_t_unit == 1:  # s
                for hdx_time in time_list:
                    secs = (np.round(hdx_time).astype(int)) * 60
                    string = f'{secs} s'
                    if PulseLabellingNow == 2:
                        string += f' - {ref_time} m'
                    labels.append(string)
            if what_t_unit == 2:  # m
                for hdx_time in time_list:
                    if hdx_time >= 1:
                        mins = np.round(hdx_time).astype(int)
                        string = f'{mins} m'
                    else:
                        string = f'{hdx_time} m'.rstrip('0').rstrip('.')
                    if PulseLabellingNow == 2:
                        string += f' - {ref_time} m'
                    labels.append(string)
            if what_t_unit == 3:  # h
                for hdx_time in time_list:
                    if hdx_time >= 60:
                        hours = (np.round(hdx_time).astype(int)) / 60
                        string = f'{hours} h'.rstrip('0').rstrip('.')
                    else:
                        hours = (np.round(hdx_time).astype(int)) / 60
                        string = f'{hours:.2f} h'.rstrip('0').rstrip('.')
                    if PulseLabellingNow == 2:
                        string += f' - {ref_time} m'
                    labels.append(string)
        return(labels)

    if uploaded_settings != 1:
        hmspacerthick = 4
        hmlcolor = greyc
        hmsepthick = 4
        hmsepcolorc = blackc
        hmbordthick = 10
        hmbordcolorc = blackc
        hmtl = 25 
        hmtw = 4 
        hmtcolor = blackc
        hmtcolor_labels = blackc
        usehmthick = float(request.form.get('usehmthick', 0))
        usestatesep = float(request.form.get('stackedplotdivl',0))
        usebordchange = float(request.form.get('hmplotbord',0))
        usehmtickchange = float(request.form.get('changehmtick',0))
        if usehmthick == 1:
            hmspacerthick = float(request.form.get('hmthickness',4)) 
            hmcolordivide = float(request.form.get('hmcolor',1))
            if hmcolordivide == 1: #White
                hmlcolor = whitec 
            elif hmcolordivide == 2: #Black
                hmlcolor = blackc
            elif hmcolordivide == 3: #Blue
                hmlcolor = dbluec
            elif hmcolordivide == 4: #Red
                hmlcolor = dredc
            elif hmcolordivide == 5: #Grey
                hmlcolor = greyc
        if usestatesep == 1:
            hmsepthick = float(request.form.get('hmsepthickness',4)) 
            hmsepcolor = float(request.form.get('hmsepcolor',2))
            if hmsepcolor == 1: #White
                hmsepcolorc = whitec 
            elif hmsepcolor == 2: #Black
                hmsepcolorc = blackc
            elif hmsepcolor == 3: #Blue
                hmsepcolorc = dbluec
            elif hmsepcolor == 4: #Red
                hmsepcolorc = dredc
            elif hmsepcolor == 5: #Grey
                hmsepcolorc = greyc    
        if usebordchange == 1:
            hmbordthick = float(request.form.get('hmbordthickness',10))
            hmbordcolor = float(request.form.get('hmbordcolor',2))
            if hmbordcolor == 1: #White
                hmbordcolorc = whitec 
            elif hmbordcolor == 2: #Black
                hmbordcolorc = blackc
            elif hmbordcolor == 3: #Blue
                hmbordcolorc = dbluec
            elif hmbordcolor == 4: #Red
                hmbordcolorc = dredc
            elif hmbordcolor == 5: #Grey
                hmbordcolorc = greyc 
        if usehmtickchange == 1:
            hmtw = float(request.form.get('hmtickwidth',4))
            hmtl = float(request.form.get('hmticklength',25))
            hmtickcolor = float(request.form.get('hmtickcolor',2))
            hmtickcolorlabel = float(request.form.get('hmtickcolorlabel',2))
            #Color Ticks
            if hmtickcolor == 1: #White
                hmtcolor = whitec
            elif hmtickcolor == 2: #Black
                hmtcolor = blackc
            elif hmtickcolor == 3: #Blue
                hmtcolor = dbluec
            elif hmtickcolor == 4: #Red
                hmtcolor = dredc
            elif hmtickcolor == 5: #Grey
                hmtcolor = greyc
            #Color Labels
            if hmtickcolorlabel == 1: #White
                hmtcolor_labels = whitec 
            elif hmtickcolorlabel == 2: #Black
                hmtcolor_labels = blackc
            elif hmtickcolorlabel == 3: #Blue
                hmtcolor_labels = dbluec
            elif hmtickcolorlabel == 4: #Red
                hmtcolor_labels = dredc
            elif hmtickcolorlabel == 5: #Grey
                hmtcolor_labels = greyc
        padthick = 20 #padding for titles
        spadthick = 10 #padding for axis labels
        usealtpad = float(request.form.get('usealtpad',20))
        if usealtpad == 1:
            padthick = float(request.form.get('altpad',20))
            spadthick = float(request.form.get('altpads',20))
    #color_sections = [                    AN EXAMPLE OF THE LIST USED HERE
    #    (5, 10, 'red', 'Section 1'),
    #    (15, 20, 'blue', 'Section 2'),
    #    (25, 30, 'green', 'Section 3')
    #]
    if uploaded_settings != 1:
        domainlabel = float(request.form.get('domainlabel',0))
        if domainlabel == 1:
            numdomain = int(request.form.get('numdomain', 4))
            padthick = padthick + 40 + (font_size_title-32)*2
            for i in range(1, numdomain+4):
                domainName = request.form.get(f'domainName{i}', '')
                dompepst = float(request.form.get(f'dompepst{i}', 0)) - 1
                dompepend = float(request.form.get(f'dompepend{i}', 0))
                domColour = request.form.get(f'domColour{i}', '#000000')
                if domainName:
                    color_sections.append((dompepst, dompepend, domColour, domainName))
    blockpositioning = -2.55 + (-font_size+28)/15
    blocktextpos = -0.4 - 0.5*(font_size_title-32)/font_size_title 
    blockcolthick = 0.25
    heatmap_buffer = io.BytesIO()

    
    ##################################################################################
    ###########################      HEATMAP PLOTTING      ###########################
    ##################################################################################

    
    def h_plot(pdf, data, title, pept, color_sections):
        global output_bitmap_h_count
        if isinstance(data, np.ma.MaskedArray):
            data = data.filled(0)
        all_prot = all_deltas.groupby('Protein')
        if expType != 3:
            ###############################  CONTINUOUS LABELLING PLOTTING  ###############################
            if len(all_prot) > 1 and PDFgeneration != 1 and plot_stacked == 1:
                with zipfile.ZipFile(heatmap_buffer, 'a', zipfile.ZIP_DEFLATED, False) as zip_file:
                    # Generate the heatmap for each protein
                    fig_margin_x_l = 3
                    fig_margin_x_r = 1
                    fig_margin_y_t = 4
                    fig_margin_y_b = 3
                    fig_length_x = data.shape[0]
                    fig_length_y = data.shape[1]
                    fig_x = fig_margin_x_l + fig_length_x + fig_margin_x_r
                    fig_y = fig_margin_y_b + fig_length_y + fig_margin_y_t
                    cbar_spacer = 0.5  # Gap between heatmap plot and colorbar
                    cbar_length = min(len(colors), fig_length_x)  # 10 or fig_length_x, whichever is smaller
                    # Make tick labels
                    pept['Start'] = pept['Start'].round(0).astype(int)
                    pept['End'] = pept['End'].round(0).astype(int)
                    pept_tick_labels = pept['Start'].map(str) + '-' + pept['End'].map(str)
                    time_tick_labels = make_time_tick_labels(data.columns.get_level_values('Exposure'))
                    fig = plt.figure(figsize=(fig_x, fig_y), dpi=output_bitmap_dpi)
                    ax = plt.axes((fig_margin_x_l/fig_x, fig_margin_y_b/fig_y, fig_length_x/fig_x, fig_length_y/fig_y))
                    ax_4_cbar = plt.axes(((fig_margin_x_l+(fig_length_x-cbar_length)/2)/fig_x, (fig_margin_y_b - cbar_spacer - cbar_length/len(colors))/fig_y, cbar_length/fig_x, cbar_length/len(colors)/fig_y))
                    if relativeUptakeCalc == 1:
                        hmlabelforbar = u'Relative Δ Deuterium Uptake'
                    else:
                        hmlabelforbar = u'Δ Deuterium Uptake'
                    sns.heatmap(data.transpose(), ax=ax, cmap=colormap, norm=my_norm, xticklabels=pept_tick_labels, yticklabels=time_tick_labels, linewidths=hmspacerthick, linecolor=hmlcolor, 
                                square=True, cbar_ax=ax_4_cbar, cbar_kws={"orientation": "horizontal", 'label': hmlabelforbar}, annot = addvalHM)
                    plt.sca(ax)
                    ax.set_facecolor(c_missing)
                    if editing_figtitle == 1:
                        if uploaded_settings != 1:
                            newtitle = str(request.form.get('alttitleword',title))
                        plt.title(newtitle, fontsize=font_size_title, color=fontcolor, pad=padthick)
                    else:
                        plt.title(title, fontsize=font_size_title, color=fontcolor, pad=padthick)
                    plt.setp(ax.spines.values(), linewidth=hmsepthick, color=hmsepcolorc)
                    plt.ylabel(new_y_axis_title, labelpad=spadthick)
                    plt.xlabel(new_x_axis_title, labelpad=spadthick)
                    ax.xaxis.tick_top()
                    ax.xaxis.set_label_position('top')
                    ax.tick_params(axis='x', labelrotation=90, length = hmtl, width = hmtw, colors=hmtcolor)
                    ax.tick_params(axis='y', labelrotation=0, length = hmtl, width = hmtw, colors=hmtcolor)
                    ax.xaxis.set_tick_params(labelcolor=hmtcolor_labels)  # X-axis tick label color
                    ax.yaxis.set_tick_params(labelcolor=hmtcolor_labels)  # Y-axis tick label color
                    ax.set_xticklabels(pept_tick_labels, fontsize = font_size_ticklabel)
                    ax.set_yticklabels(time_tick_labels, fontsize = font_size_ticklabel)
                    # Placing black line around the plot:
                    ax.axhline(y=0, color=hmbordcolorc, linewidth=hmbordthick)
                    ax.axhline(y=fig_length_y, color=hmbordcolorc, linewidth=hmbordthick)
                    ax.axvline(x=0, color=hmbordcolorc, linewidth=hmbordthick)
                    ax.axvline(x=fig_length_x, color=hmbordcolorc, linewidth=hmbordthick)
                    # Add separator lines between states:
                    states = data.columns.get_level_values('State').unique()
                    sep_line_y = 0
                    for i in range(len(states)-1):
                        exposures = data.loc[:,('Delta_d_uptake','mean',states[i],slice(None))].columns.get_level_values('Exposure')
                        sep_line_y += len(exposures)
                        ax.axhline(y=sep_line_y, color=hmsepcolorc, linewidth=hmsepthick)
                    # Add colored blocks above the peptide axis
                    for section in color_sections:
                        start, end, color, text = section
                        # Adjust y position of the rectangle and annotation to be above the top ticks
                        rect_y = blockpositioning  # Position above the top of the heatmap
                        rect = plt.Rectangle((start, rect_y), end-start, blockcolthick, color=color, clip_on=False)
                        ax.add_patch(rect)
                        text_y = blocktextpos + rect_y + blockcolthick / 2
                        ax.annotate(text, xy=((start+end)/2, text_y), xytext=(0,0), textcoords='offset points',
                                    ha='center', va='center', fontsize=font_size_title, color=fontcolor, clip_on=False, annotation_clip=False)
                        #print('adding section')
                    #pdf.savefig(bbox_inches='tight')  # Saves the current figure into a pdf page
                    if output_bitmap:
                        if plot_separate == 1:
                            output_bitmap_file = output_bitmap_name + '_h_' + str(output_bitmap_h_count) + '_' + prot + '_' + state + '.' + output_bitmap_format
                        else: 
                            output_bitmap_file = output_bitmap_name + '_h_' + str(output_bitmap_h_count) + '_' + prot + '.' + output_bitmap_format
                        plt.savefig(output_bitmap_file, dpi=output_bitmap_dpi, bbox_inches='tight')
                        output_bitmap_h_count += 1
                    if plot_separate == 1:
                        output_bitmap_file = output_bitmap_name + '_h_' + str(output_bitmap_h_count) + '_' + state + '.' + output_bitmap_format
                    else: 
                        output_bitmap_file = output_bitmap_name + '_h_' + str(output_bitmap_h_count) + '.' + output_bitmap_format
                    plt.savefig(output_bitmap_file, format='png', dpi=output_bitmap_dpi, bbox_inches='tight')
                    plt.close()
                    zip_file.write(output_bitmap_file)  # Add the file to the zip archive
                    # Remove the file after adding it to the zip archive
                    Path(output_bitmap_file).unlink()
                with open('heatmaps.zip', 'wb') as f:
                    f.write(heatmap_buffer.getvalue())
            else:
                with zipfile.ZipFile(heatmap_buffer, 'a', zipfile.ZIP_DEFLATED, False) as zip_file:
                    fig_margin_x_l = 3
                    fig_margin_x_r = 1
                    fig_margin_y_t = 4
                    fig_margin_y_b = 3
                    fig_length_x = data.shape[0]
                    fig_length_y = data.shape[1]
                    fig_x = fig_margin_x_l + fig_length_x + fig_margin_x_r
                    fig_y = fig_margin_y_b + fig_length_y + fig_margin_y_t
                    cbar_spacer = 0.5  # Gap between heatmap plot and colorbar
                    cbar_length = min(len(colors), fig_length_x)  # 10 or fig_length_x, whichever is smaller
                    # Make tick labels
                    pept['Start'] = pept['Start'].round(0).astype(int)
                    pept['End'] = pept['End'].round(0).astype(int)
                    pept_tick_labels = pept['Start'].map(str) + '-' + pept['End'].map(str)
                    time_tick_labels = make_time_tick_labels(data.columns.get_level_values('Exposure'))
                    fig = plt.figure(figsize=(fig_x, fig_y), dpi=output_bitmap_dpi)
                    ax = plt.axes((fig_margin_x_l/fig_x, fig_margin_y_b/fig_y, fig_length_x/fig_x, fig_length_y/fig_y))
                    ax_4_cbar = plt.axes(((fig_margin_x_l+(fig_length_x-cbar_length)/2)/fig_x, (fig_margin_y_b - cbar_spacer - cbar_length/len(colors))/fig_y, cbar_length/fig_x, cbar_length/len(colors)/fig_y))
                    if relativeUptakeCalc == 1:
                        hmlabelforbar = u'Relative Δ Deuterium Uptake'
                    else:
                        hmlabelforbar = u'Δ Deuterium Uptake'
                    sns.heatmap(data.transpose(), ax=ax, cmap=colormap, norm=my_norm, xticklabels=pept_tick_labels, yticklabels=time_tick_labels, linewidths=hmspacerthick, linecolor=hmlcolor, 
                                square=True, cbar_ax=ax_4_cbar, cbar_kws={"orientation": "horizontal", 'label': hmlabelforbar}, annot = addvalHM)
                    plt.sca(ax)
                    ax.set_facecolor(c_missing)
                    if editing_figtitle == 1:
                        if uploaded_settings != 1:
                            newtitle = str(request.form.get('alttitleword','Invalid Text'))
                        plt.title(newtitle, fontsize=font_size_title, color=fontcolor, pad=padthick)
                    else:
                        plt.title(title, fontsize=font_size_title, color=fontcolor, pad=padthick)
                    plt.setp(ax.spines.values(), linewidth=hmsepthick, color=hmsepcolorc)
                    plt.ylabel(new_y_axis_title, labelpad=spadthick)
                    plt.xlabel(new_x_axis_title, labelpad=spadthick)
                    ax.xaxis.tick_top()
                    ax.xaxis.set_label_position('top')
                    ax.tick_params(axis='x', labelrotation=90, length = hmtl, width = hmtw, colors=hmtcolor)
                    ax.tick_params(axis='y', labelrotation=0, length = hmtl, width = hmtw, colors=hmtcolor)
                    ax.xaxis.set_tick_params(labelcolor=hmtcolor_labels)  # X-axis tick label color
                    ax.yaxis.set_tick_params(labelcolor=hmtcolor_labels)  # Y-axis tick label color
                    ax.set_xticklabels(pept_tick_labels, fontsize = font_size_ticklabel)
                    ax.set_yticklabels(time_tick_labels, fontsize = font_size_ticklabel)
                    # Placing black line around the plot:
                    ax.axhline(y=0, color=hmbordcolorc, linewidth=hmbordthick)
                    ax.axhline(y=fig_length_y, color=hmbordcolorc, linewidth=hmbordthick)
                    ax.axvline(x=0, color=hmbordcolorc, linewidth=hmbordthick)
                    ax.axvline(x=fig_length_x, color=hmbordcolorc, linewidth=hmbordthick)
                    # Add separator lines between states:
                    states = data.columns.get_level_values('State').unique()
                    sep_line_y = 0
                    for i in range(len(states)-1):
                        exposures = data.loc[:,('Delta_d_uptake','mean',states[i],slice(None))].columns.get_level_values('Exposure')
                        sep_line_y += len(exposures)
                        ax.axhline(y=sep_line_y, color=hmsepcolorc, linewidth=hmsepthick)
                    # Add colored blocks above the peptide axis
                    for section in color_sections:
                        start, end, color, text = section
                        # Adjust y position of the rectangle and annotation to be above the top ticks
                        rect_y = blockpositioning  # Position above the top of the heatmap
                        rect = plt.Rectangle((start, rect_y), end-start, blockcolthick, color=color, clip_on=False)
                        ax.add_patch(rect)
                        text_y = blocktextpos + rect_y + blockcolthick / 2 
                        ax.annotate(text, xy=((start+end)/2, text_y), xytext=(0,0), textcoords='offset points',
                                    ha='center', va='center', fontsize=font_size_title, color=fontcolor, clip_on=False, annotation_clip=False)
                        print('adding section')
                    #pdf.savefig(bbox_inches='tight')  # Saves the current figure into a pdf page
                    if output_bitmap:
                        if plot_separate == 1:
                            output_bitmap_file = output_bitmap_name + '_h_' + str(output_bitmap_h_count) + '_' + prot + '_' + state + '.' + output_bitmap_format
                        else: 
                            output_bitmap_file = output_bitmap_name + '_h_' + str(output_bitmap_h_count) + '_' + prot + '.' + output_bitmap_format
                        plt.savefig(output_bitmap_file, dpi=output_bitmap_dpi, bbox_inches='tight')
                        output_bitmap_h_count += 1
                    if plot_separate == 1:
                        output_bitmap_file = output_bitmap_name + '_h_' + str(output_bitmap_h_count) + '_' + prot + '_' + state + '.' + output_bitmap_format
                    else: 
                        output_bitmap_file = output_bitmap_name + '_h_' + str(output_bitmap_h_count) + '_' + prot + '.' + output_bitmap_format
                    if PDFgeneration == 1:
                        plt.savefig(output_bitmap_file, format='pdf', dpi=output_bitmap_dpi, bbox_inches='tight')
                    else:
                        plt.savefig(output_bitmap_file, format='png', dpi=output_bitmap_dpi, bbox_inches='tight')
                    plt.close()
                    zip_file.write(output_bitmap_file)  # Add the file to the zip archive
                with open('heatmaps.zip', 'wb') as f:
                    f.write(heatmap_buffer.getvalue())
        ###############################  RELATIVE FRACTIONAL UPTAKE PLOTTING  ###############################
        elif expType == 3: #RFU
            if len(all_prot) > 1 and PDFgeneration != 1 and plot_stacked == 1:
                with zipfile.ZipFile(heatmap_buffer, 'a', zipfile.ZIP_DEFLATED, False) as zip_file:
                    # Generate the heatmap for each protein
                    fig_margin_x_l = 3
                    fig_margin_x_r = 1
                    fig_margin_y_t = 4
                    fig_margin_y_b = 3
                    fig_length_x = data.shape[0]
                    fig_length_y = data.shape[1]
                    fig_x = fig_margin_x_l + fig_length_x + fig_margin_x_r
                    fig_y = fig_margin_y_b + fig_length_y + fig_margin_y_t
                    cbar_spacer = 0.5  # Gap between heatmap plot and colorbar
                    cbar_length = min(len(colors), fig_length_x)  # 10 or fig_length_x, whichever is smaller
                    # Make tick labels
                    pept['Start'] = pept['Start'].round(0).astype(int)
                    pept['End'] = pept['End'].round(0).astype(int)
                    pept_tick_labels = pept['Start'].map(str) + '-' + pept['End'].map(str)
                    time_tick_labels = make_time_tick_labels(data.columns.get_level_values('Exposure'))
                    fig = plt.figure(figsize=(fig_x, fig_y), dpi=output_bitmap_dpi)
                    ax = plt.axes((fig_margin_x_l/fig_x, fig_margin_y_b/fig_y, fig_length_x/fig_x, fig_length_y/fig_y))
                    ax_4_cbar = plt.axes(((fig_margin_x_l+(fig_length_x-cbar_length)/2)/fig_x, (fig_margin_y_b - cbar_spacer - cbar_length/len(colors))/fig_y, cbar_length/fig_x, cbar_length/len(colors)/fig_y))
                    if percentUseRFU != 100:
                        hmlabelforbar = u'Relative Deuterium Uptake'
                    elif absoluteUptakeValues == 1:
                        hmlabelforbar = u'Absolute Deuterium Uptake'
                    else:
                        hmlabelforbar = u'Percent Relative Deuterium Uptake'
                    sns.heatmap(data.transpose(), ax=ax, cmap=colormap, norm=my_norm, xticklabels=pept_tick_labels, yticklabels=time_tick_labels, linewidths=hmspacerthick, linecolor=hmlcolor, 
                                square=True, cbar_ax=ax_4_cbar, cbar_kws={"orientation": "horizontal", 'label': hmlabelforbar}, annot = addvalHM)
                    plt.sca(ax)
                    ax.set_facecolor(c_missing)
                    if editing_figtitle == 1:
                        if uploaded_settings != 1:
                            newtitle = str(request.form.get('alttitleword',title))
                        plt.title(newtitle, fontsize=font_size_title, color=fontcolor, pad=padthick)
                    else:
                        plt.title(title, fontsize=font_size_title, color=fontcolor, pad=padthick)
                    plt.setp(ax.spines.values(), linewidth=hmsepthick, color=hmsepcolorc)
                    plt.ylabel(new_y_axis_title, labelpad=spadthick)
                    plt.xlabel(new_x_axis_title, labelpad=spadthick)
                    ax.xaxis.tick_top()
                    ax.xaxis.set_label_position('top')
                    ax.tick_params(axis='x', labelrotation=90, length = hmtl, width = hmtw, colors=hmtcolor)
                    ax.tick_params(axis='y', labelrotation=0, length = hmtl, width = hmtw, colors=hmtcolor)
                    ax.xaxis.set_tick_params(labelcolor=hmtcolor_labels)  # X-axis tick label color
                    ax.yaxis.set_tick_params(labelcolor=hmtcolor_labels)  # Y-axis tick label color
                    ax.set_xticklabels(pept_tick_labels, fontsize = font_size_ticklabel)
                    ax.set_yticklabels(time_tick_labels, fontsize = font_size_ticklabel)
                    # Placing black line around the plot:
                    ax.axhline(y=0, color=hmbordcolorc, linewidth=hmbordthick)
                    ax.axhline(y=fig_length_y, color=hmbordcolorc, linewidth=hmbordthick)
                    ax.axvline(x=0, color=hmbordcolorc, linewidth=hmbordthick)
                    ax.axvline(x=fig_length_x, color=hmbordcolorc, linewidth=hmbordthick)
                    # Add separator lines between states:
                    states = data.columns.get_level_values('State').unique()
                    sep_line_y = 0
                    for i in range(len(states)-1):
                        exposures = data.loc[:,('d_uptake','mean',states[i],slice(None))].columns.get_level_values('Exposure')
                        sep_line_y += len(exposures)
                        ax.axhline(y=sep_line_y, color=hmsepcolorc, linewidth=hmsepthick)
                    # Add colored blocks above the peptide axis
                    for section in color_sections:
                        start, end, color, text = section
                        # Adjust y position of the rectangle and annotation to be above the top ticks
                        rect_y = blockpositioning  # Position above the top of the heatmap
                        rect = plt.Rectangle((start, rect_y), end-start, blockcolthick, color=color, clip_on=False)
                        ax.add_patch(rect)
                        text_y = blocktextpos + rect_y + blockcolthick / 2
                        ax.annotate(text, xy=((start+end)/2, text_y), xytext=(0,0), textcoords='offset points',
                                    ha='center', va='center', fontsize=font_size_title, color=fontcolor, clip_on=False, annotation_clip=False)
                        #print('adding section')
                    #pdf.savefig(bbox_inches='tight')  # Saves the current figure into a pdf page
                    if output_bitmap:
                        if plot_separate == 1:
                            output_bitmap_file = output_bitmap_name + '_h_' + str(output_bitmap_h_count) + '_' + prot + '_' + state + '.' + output_bitmap_format
                        else: 
                            output_bitmap_file = output_bitmap_name + '_h_' + str(output_bitmap_h_count) + '_' + prot + '.' + output_bitmap_format
                        plt.savefig(output_bitmap_file, dpi=output_bitmap_dpi, bbox_inches='tight')
                        output_bitmap_h_count += 1
                    if plot_separate == 1:
                        output_bitmap_file = output_bitmap_name + '_h_' + str(output_bitmap_h_count) + '_' + state + '.' + output_bitmap_format
                    else: 
                        output_bitmap_file = output_bitmap_name + '_h_' + str(output_bitmap_h_count) + '.' + output_bitmap_format
                    plt.savefig(output_bitmap_file, format='png', dpi=output_bitmap_dpi, bbox_inches='tight')
                    plt.close()
                    zip_file.write(output_bitmap_file)  # Add the file to the zip archive
                    # Remove the file after adding it to the zip archive
                    Path(output_bitmap_file).unlink()
                with open('heatmaps.zip', 'wb') as f:
                    f.write(heatmap_buffer.getvalue())
            else:
                with zipfile.ZipFile(heatmap_buffer, 'a', zipfile.ZIP_DEFLATED, False) as zip_file:
                    fig_margin_x_l = 3
                    fig_margin_x_r = 1
                    fig_margin_y_t = 4
                    fig_margin_y_b = 3
                    fig_length_x = data.shape[0]
                    fig_length_y = data.shape[1]
                    fig_x = fig_margin_x_l + fig_length_x + fig_margin_x_r
                    fig_y = fig_margin_y_b + fig_length_y + fig_margin_y_t
                    cbar_spacer = 0.5  # Gap between heatmap plot and colorbar
                    cbar_length = min(len(colors), fig_length_x)  # 10 or fig_length_x, whichever is smaller
                    # Make tick labels
                    pept['Start'] = pept['Start'].round(0).astype(int)
                    pept['End'] = pept['End'].round(0).astype(int)
                    pept_tick_labels = pept['Start'].map(str) + '-' + pept['End'].map(str)
                    time_tick_labels = make_time_tick_labels(data.columns.get_level_values('Exposure'))
                    fig = plt.figure(figsize=(fig_x, fig_y), dpi=output_bitmap_dpi)
                    ax = plt.axes((fig_margin_x_l/fig_x, fig_margin_y_b/fig_y, fig_length_x/fig_x, fig_length_y/fig_y))
                    ax_4_cbar = plt.axes(((fig_margin_x_l+(fig_length_x-cbar_length)/2)/fig_x, (fig_margin_y_b - cbar_spacer - cbar_length/len(colors))/fig_y, cbar_length/fig_x, cbar_length/len(colors)/fig_y))
                    if percentUseRFU != 100:
                        hmlabelforbar = u'Relative Deuterium Uptake'
                    elif absoluteUptakeValues == 1:
                        hmlabelforbar = u'Absolute Deuterium Uptake'
                    else:
                        hmlabelforbar = u'Percent Relative Deuterium Uptake'
                    sns.heatmap(data.transpose(), ax=ax, cmap=colormap, norm=my_norm, xticklabels=pept_tick_labels, yticklabels=time_tick_labels, linewidths=hmspacerthick, linecolor=hmlcolor, 
                                square=True, cbar_ax=ax_4_cbar, cbar_kws={"orientation": "horizontal", 'label': hmlabelforbar}, annot = addvalHM)
                    plt.sca(ax)
                    ax.set_facecolor(c_missing)
                    if editing_figtitle == 1:
                        if uploaded_settings != 1:
                            newtitle = str(request.form.get('alttitleword','Invalid Text'))
                        plt.title(newtitle, fontsize=font_size_title, color=fontcolor, pad=padthick)
                    else:
                        plt.title(title, fontsize=font_size_title, color=fontcolor, pad=padthick)
                    plt.setp(ax.spines.values(), linewidth=hmsepthick, color=hmsepcolorc)
                    plt.ylabel(new_y_axis_title, labelpad=spadthick)
                    plt.xlabel(new_x_axis_title, labelpad=spadthick)
                    ax.xaxis.tick_top()
                    ax.xaxis.set_label_position('top')
                    ax.tick_params(axis='x', labelrotation=90, length = hmtl, width = hmtw, colors=hmtcolor)
                    ax.tick_params(axis='y', labelrotation=0, length = hmtl, width = hmtw, colors=hmtcolor)
                    ax.xaxis.set_tick_params(labelcolor=hmtcolor_labels)  # X-axis tick label color
                    ax.yaxis.set_tick_params(labelcolor=hmtcolor_labels)  # Y-axis tick label color
                    ax.set_xticklabels(pept_tick_labels, fontsize = font_size_ticklabel)
                    ax.set_yticklabels(time_tick_labels, fontsize = font_size_ticklabel)
                    # Placing black line around the plot:
                    ax.axhline(y=0, color=hmbordcolorc, linewidth=hmbordthick)
                    ax.axhline(y=fig_length_y, color=hmbordcolorc, linewidth=hmbordthick)
                    ax.axvline(x=0, color=hmbordcolorc, linewidth=hmbordthick)
                    ax.axvline(x=fig_length_x, color=hmbordcolorc, linewidth=hmbordthick)
                    # Add separator lines between states:
                    states = data.columns.get_level_values('State').unique()
                    sep_line_y = 0
                    for i in range(len(states)-1):
                        exposures = data.loc[:,('d_uptake','mean',states[i],slice(None))].columns.get_level_values('Exposure')
                        sep_line_y += len(exposures)
                        ax.axhline(y=sep_line_y, color=hmsepcolorc, linewidth=hmsepthick)
                    # Add colored blocks above the peptide axis
                    for section in color_sections:
                        start, end, color, text = section
                        # Adjust y position of the rectangle and annotation to be above the top ticks
                        rect_y = blockpositioning  # Position above the top of the heatmap
                        rect = plt.Rectangle((start, rect_y), end-start, blockcolthick, color=color, clip_on=False)
                        ax.add_patch(rect)
                        text_y = blocktextpos + rect_y + blockcolthick / 2 
                        ax.annotate(text, xy=((start+end)/2, text_y), xytext=(0,0), textcoords='offset points',
                                    ha='center', va='center', fontsize=font_size_title, color=fontcolor, clip_on=False, annotation_clip=False)
                        print('adding section')
                    #pdf.savefig(bbox_inches='tight')  # Saves the current figure into a pdf page
                    if output_bitmap:
                        if plot_separate == 1:
                            output_bitmap_file = output_bitmap_name + '_h_' + str(output_bitmap_h_count) + '_' + prot + '_' + state + '.' + output_bitmap_format
                        else: 
                            output_bitmap_file = output_bitmap_name + '_h_' + str(output_bitmap_h_count) + '_' + prot + '.' + output_bitmap_format
                        plt.savefig(output_bitmap_file, dpi=output_bitmap_dpi, bbox_inches='tight')
                        output_bitmap_h_count += 1
                    if plot_separate == 1:
                        output_bitmap_file = output_bitmap_name + '_h_' + str(output_bitmap_h_count) + '_' + prot + '_' + state + '.' + output_bitmap_format
                    else: 
                        output_bitmap_file = output_bitmap_name + '_h_' + str(output_bitmap_h_count) + '_' + prot + '.' + output_bitmap_format
                    if PDFgeneration == 1:
                        plt.savefig(output_bitmap_file, format='pdf', dpi=output_bitmap_dpi, bbox_inches='tight')
                    else:
                        plt.savefig(output_bitmap_file, format='png', dpi=output_bitmap_dpi, bbox_inches='tight')
                    plt.close()
                    zip_file.write(output_bitmap_file)  # Add the file to the zip archive
                with open('heatmaps.zip', 'wb') as f:
                    f.write(heatmap_buffer.getvalue())

    
    ##################################################################################
    ###########################      VERTICAL PLOTTING      ##########################
    ##################################################################################

    
    def v_plot(pdf, data, title, pept):
        global output_bitmap_v_count
        all_prot = all_deltas.groupby('Protein')
        if isinstance(data, np.ma.MaskedArray):
            data = data.filled(0)
        if len(all_prot) > 1 and PDFgeneration != 1 and plot_stacked == 1:
            with zipfile.ZipFile(heatmap_buffer, 'a', zipfile.ZIP_DEFLATED, False) as zip_file:
                fig_margin_x_l = 3
                fig_margin_x_r = 3
                fig_margin_y_t = 3
                fig_margin_y_b = 1
                fig_length_x = data.shape[1]
                fig_length_y = data.shape[0]
                fig_x = fig_margin_x_l + fig_length_x + fig_margin_x_r
                fig_y = fig_margin_y_b + fig_length_y + fig_margin_y_t
                cbar_spacer = 0.5  # Gap between heatmap plot and colorbar
                cbar_length  = min(len(colors), fig_length_y)  # 10 or fig_length_x, whichever is smaller
                # make tick labels:
                pept['Start'] = pept['Start'].round(0).astype(int)
                pept['End'] = pept['End'].round(0).astype(int)
                pept_tick_labels = pept['Start'].map(str) + '-' + pept['End'].map(str)
                time_tick_labels = make_time_tick_labels( data.columns.get_level_values('Exposure') )
                fig = plt.figure(figsize=(fig_x , fig_y), dpi=output_bitmap_dpi)
                ax = plt.axes((fig_margin_x_l/fig_x, fig_margin_y_b/fig_y, fig_length_x/fig_x, fig_length_y/fig_y ))
                ax_4_cbar = plt.axes(((fig_margin_x_l+fig_length_x+cbar_spacer)/fig_x, (fig_margin_y_b + (fig_length_y-cbar_length)/2)/fig_y, cbar_length/len(colors)/fig_x, cbar_length/fig_y ))
                if relativeUptakeCalc == 1:
                    hmlabelforbar = u'Relative Δ Deuterium Uptake'
                else:
                    hmlabelforbar = u'Δ Deuterium Uptake'
                sns.heatmap(data, ax=ax, cmap=colormap, norm=my_norm, xticklabels=time_tick_labels, yticklabels=pept_tick_labels, linewidths=hmspacerthick, linecolor=hmlcolor,
                            square=True, cbar_ax=ax_4_cbar, cbar_kws={"orientation": "vertical", 'label': hmlabelforbar}, annot = addvalHM)
                plt.sca(ax)
                ax.set_facecolor(c_missing)
                if editing_figtitle == 1:
                    if uploaded_settings != 1:
                        newtitle = str(request.form.get('alttitleword','Invalid Text'))
                    plt.title(newtitle, fontsize=font_size_title, color=fontcolor, pad=padthick)
                else:
                    plt.title(title, fontsize=font_size_title, color=fontcolor, pad=padthick)
                plt.setp(ax.spines.values(), linewidth=hmsepthick, color=hmsepcolorc)
                plt.xlabel(new_y_axis_title, labelpad=spadthick)
                plt.ylabel(new_x_axis_title, labelpad=spadthick)
                ax.xaxis.tick_top()
                ax.xaxis.set_label_position('top')
                ax.tick_params(axis = 'x', labelrotation = 30, length = hmtl, width = hmtw, colors=hmtcolor)
                ax.tick_params(axis = 'y', labelrotation = 0, length = hmtl, width = hmtw, colors=hmtcolor)
                ax.xaxis.set_tick_params(labelcolor=hmtcolor_labels)  # X-axis tick label color
                ax.yaxis.set_tick_params(labelcolor=hmtcolor_labels)  # Y-axis tick label color
                ax.set_yticklabels(pept_tick_labels, fontsize = font_size_ticklabel)
                ax.set_xticklabels(time_tick_labels, fontsize = font_size_ticklabel)
                #Placing black line around the plot:
                ax.axhline(y=0, color=hmbordcolorc,linewidth=hmbordthick)
                ax.axhline(y=fig_length_y, color=hmbordcolorc,linewidth=hmbordthick)
                ax.axvline(x=0, color=hmbordcolorc,linewidth=hmbordthick)
                ax.axvline(x=fig_length_x, color=hmbordcolorc,linewidth=hmbordthick)
                # Add separator lines between states:
                states = data.columns.get_level_values('State').unique()
                sep_line_x = 0
                for i in range(len(states)-1):
                    exposures = data.loc[:,('Delta_d_uptake','mean',states[i],slice(None))].columns.get_level_values('Exposure')
                    sep_line_x += len(exposures)
                    ax.axvline(x=sep_line_x, color=hmsepcolorc, linewidth=hmsepthick)
                #print(pept_tick_labels)
                #print(time_tick_labels)
                #pdf.savefig()  # saves the current figure into a pdf page
                if output_bitmap:
                    if plot_separate == 1:
                        output_bitmap_file = output_bitmap_name + '_v_' + str(output_bitmap_v_count) + '_' + prot + '_' + state + '.' + output_bitmap_format
                    else: 
                        output_bitmap_file = output_bitmap_name + '_v_' + str(output_bitmap_v_count) + '_' + prot + '.' + output_bitmap_format
                    plt.savefig(output_bitmap_file, dpi=output_bitmap_dpi)
                    output_bitmap_v_count += 1
                if plot_separate == 1:
                    output_bitmap_file = output_bitmap_name + '_v_' + str(output_bitmap_v_count) + '_' + prot + '_' + state + '.' + output_bitmap_format
                else: 
                    output_bitmap_file = output_bitmap_name + '_v_' + str(output_bitmap_v_count) + '_' + prot + '.' + output_bitmap_format
                plt.savefig(output_bitmap_file, format='png', dpi=output_bitmap_dpi)
                buffer.seek(0)    
                plt.close()
                zip_file.write(output_bitmap_file)  # Add the file to the zip archive
                # Remove the file after adding it to the zip archive
                Path(output_bitmap_file).unlink()
                #print('s4')
            # Save the zip buffer to a file or return it as a response in Flask, depending on your use case
            with open('heatmaps.zip', 'wb') as f:
                f.write(heatmap_buffer.getvalue())
        else:
            with zipfile.ZipFile(heatmap_buffer, 'a', zipfile.ZIP_DEFLATED, False) as zip_file:
                fig_margin_x_l = 3
                fig_margin_x_r = 3
                fig_margin_y_t = 3
                fig_margin_y_b = 1
                fig_length_x = data.shape[1]
                fig_length_y = data.shape[0]
                fig_x = fig_margin_x_l + fig_length_x + fig_margin_x_r
                fig_y = fig_margin_y_b + fig_length_y + fig_margin_y_t
                cbar_spacer = 0.5  # Gap between heatmap plot and colorbar
                cbar_length  = min(len(colors), fig_length_y)  # 10 or fig_length_x, whichever is smaller
                # make tick labels:
                pept['Start'] = pept['Start'].round(0).astype(int)
                pept['End'] = pept['End'].round(0).astype(int)
                pept_tick_labels = pept['Start'].map(str) + '-' + pept['End'].map(str)
                time_tick_labels = make_time_tick_labels( data.columns.get_level_values('Exposure') )
                fig = plt.figure(figsize=(fig_x , fig_y), dpi=output_bitmap_dpi)
                ax = plt.axes((fig_margin_x_l/fig_x, fig_margin_y_b/fig_y, fig_length_x/fig_x, fig_length_y/fig_y ))
                ax_4_cbar = plt.axes(((fig_margin_x_l+fig_length_x+cbar_spacer)/fig_x, (fig_margin_y_b + (fig_length_y-cbar_length)/2)/fig_y, cbar_length/len(colors)/fig_x, cbar_length/fig_y ))
                if relativeUptakeCalc == 1:
                    hmlabelforbar = u'Relative Δ Deuterium Uptake'
                else:
                    hmlabelforbar = u'Δ Deuterium Uptake'
                sns.heatmap(data, ax=ax, cmap=colormap, norm=my_norm, xticklabels=time_tick_labels, yticklabels=pept_tick_labels, linewidths=hmspacerthick, linecolor=hmlcolor, 
                            square=True, cbar_ax=ax_4_cbar, cbar_kws={"orientation": "vertical", 'label': hmlabelforbar}, annot = addvalHM)
                plt.sca(ax)
                ax.set_facecolor(c_missing)
                if editing_figtitle == 1:
                    if uploaded_settings != 1:
                        newtitle = str(request.form.get('alttitleword','Invalid Text'))
                    plt.title(newtitle, fontsize=font_size_title, color=fontcolor, pad=padthick)
                else:
                    plt.title(title, fontsize=font_size_title, color=fontcolor, pad=padthick)
                plt.setp(ax.spines.values(), linewidth=hmsepthick, color=hmsepcolorc)
                plt.xlabel(new_y_axis_title, labelpad=spadthick)
                plt.ylabel(new_x_axis_title, labelpad=spadthick)
                ax.xaxis.tick_top()
                ax.xaxis.set_label_position('top')
                ax.tick_params(axis = 'x', labelrotation = 30, length = hmtl, width = hmtw, colors=hmtcolor)
                ax.tick_params(axis = 'y', labelrotation = 0, length = hmtl, width = hmtw, colors=hmtcolor)
                ax.xaxis.set_tick_params(labelcolor=hmtcolor_labels)  # X-axis tick label color
                ax.yaxis.set_tick_params(labelcolor=hmtcolor_labels)  # Y-axis tick label color
                ax.set_yticklabels(pept_tick_labels, fontsize = font_size_ticklabel)
                ax.set_xticklabels(time_tick_labels, fontsize = font_size_ticklabel)
                #Placing black line around the plot:
                ax.axhline(y=0, color=hmbordcolorc,linewidth=hmbordthick)
                ax.axhline(y=fig_length_y, color=hmbordcolorc,linewidth=hmbordthick)
                ax.axvline(x=0, color=hmbordcolorc,linewidth=hmbordthick)
                ax.axvline(x=fig_length_x, color=hmbordcolorc,linewidth=hmbordthick)
                # Add separator lines between states:
                states = data.columns.get_level_values('State').unique()
                sep_line_x = 0
                for i in range(len(states)-1):
                    exposures = data.loc[:,('Delta_d_uptake','mean',states[i],slice(None))].columns.get_level_values('Exposure')
                    sep_line_x += len(exposures)
                    ax.axvline(x=sep_line_x, color=hmsepcolorc, linewidth=hmsepthick)
                #print(pept_tick_labels)
                #print(time_tick_labels)
                #pdf.savefig()  # saves the current figure into a pdf page
                if output_bitmap:
                    if plot_separate == 1:
                        output_bitmap_file = output_bitmap_name + '_v_' + str(output_bitmap_v_count) + '_' + prot + '_' + state + '.' + output_bitmap_format
                    else: 
                        output_bitmap_file = output_bitmap_name + '_v_' + str(output_bitmap_v_count) + '_' + prot + '.' + output_bitmap_format
                    plt.savefig(output_bitmap_file, dpi=output_bitmap_dpi)
                    output_bitmap_v_count += 1
                if plot_separate == 1:
                    output_bitmap_file = output_bitmap_name + '_v_' + str(output_bitmap_v_count) + '_' + prot + '_' + state + '.' + output_bitmap_format
                else: 
                    output_bitmap_file = output_bitmap_name + '_v_' + str(output_bitmap_v_count) + '_' + prot + '.' + output_bitmap_format
                if PDFgeneration == 1:
                    plt.savefig(output_bitmap_file, format='pdf', dpi=output_bitmap_dpi)
                else:
                    plt.savefig(output_bitmap_file, format='png', dpi=output_bitmap_dpi)
                plt.close()
                zip_file.write(output_bitmap_file)  # Add the file to the zip archive
            with open('heatmaps.zip', 'wb') as f:
                f.write(heatmap_buffer.getvalue())


    ##################################################################################
    ############################      WOODS PLOTTING      ############################
    ##################################################################################

    
    if uploaded_settings != 1:
        altwoodsthick = float(request.form.get('altwoodsthick',0))
        if altwoodsthick == 1:
            woodscolthick = float(request.form.get('woodscolthick',3))
            woodsbackthick = float(request.form.get('woodsbackthick',0.5))
        woodsplotYbyglobal = float(request.form.get('woodsplotYbyglobal',1))
        if domainlabel == 1:
            wpadthick = wpadthick + 35
        blockposwoods = max_delta + max_delta/10
    elif 'blockposwoods' not in variables:
        blockposwoods = max_delta + max_delta/10 #for w_plot
        
    
    def w_plot(pdf, data, title, pept, state, data_source, max_delta, colcutoff):
        global output_bitmap_h_count
        if isinstance(data, np.ma.MaskedArray):
            data = data.filled(0)
        unique_exposures = data.columns.get_level_values('Exposure').unique()
        numplots = len(unique_exposures)
        data_avgrot = data_source.transpose()
        first_start = None
        last_end = None
        for exposure in unique_exposures:
            plt.figure(figsize=(woodsx, woodsy))
            if expType != 3:
                subset = data.xs(('Delta_d_uptake', 'mean', state, exposure), level=[None, 'Parameter', 'State', 'Exposure'], axis=1) 
            else:
                subset = data.xs(('d_uptake', 'mean', state, exposure), level=[None, 'Parameter', 'State', 'Exposure'], axis=1) 
            if colcutopt == 0:
                colcutoff = 0.5
            if nolines == 0:
                plt.axhline(y = 0, color = blackc) #0 line
                plt.axhline(y = colcutoff*percentUseRFU, color = blackc, linestyle = 'dashed')
                plt.axhline(y = -colcutoff*percentUseRFU, color = blackc, linestyle = 'dashed')
                plt.axhline(y = 0.75*colcutoff*percentUseRFU, color = blackc, linestyle = 'dotted') 
                plt.axhline(y = -0.75*colcutoff*percentUseRFU, color = blackc, linestyle = 'dotted')
            #print('Values for pept:', pept)
            #print('Index for pept:', pept.index)
            if PerResidueMap != 1: 
                for peptide in pept.index.get_level_values('Sequence').unique():
                    if peptide in subset.index:
                        mean_values = subset.loc[peptide].values
                        start_end_values = pept.loc[peptide, ['Start', 'End']].values
                        starts = start_end_values[0]
                        ends = start_end_values[1]
                        startsl = []
                        endsl = []
                        startsl.append(start_end_values[0])
                        endsl.append(start_end_values[1])
                        if first_start is None:
                            first_start = startsl[0]
                        if last_end is None or last_end < endsl[-1]:
                            last_end = endsl[-1]
                        # Plot the line for current peptide
                        if not pd.isna(mean_values).all():
                            mean_value = mean_values[0]
                            if color_by_heatmap == 1:
                                cmap = colormap
                                norm = my_norm
                                color = cmap(norm(mean_value))
                            elif mean_values > colcutoff:
                                color = woodscolpos
                            elif mean_values < -colcutoff:
                                color = woodscolneg
                            else: 
                                color = woodscolneu
                            plt.plot([starts, ends], [mean_value, mean_value], label=f"{peptide}", color=color, lw=woodscolthick, path_effects=[pe.Stroke(linewidth=woodscolthick+woodsbackthick, foreground='black'), pe.Normal()])
                        else:
                            print(f"Skipping plotting for peptide {peptide}: data length mismatch.")
            else:
                for peptide in pept['Start'].unique():
                    if peptide in subset.index:
                        mean_values = subset.loc[peptide].values
                        start_end_values = pept.loc[peptide, ['Start', 'End']].values
                        starts = start_end_values[0]
                        ends = start_end_values[1]+1 
                        startsl = []
                        endsl = []
                        startsl.append(start_end_values[0])
                        endsl.append(start_end_values[1]+1) 
                        if first_start is None:
                            first_start = startsl[0]
                        if last_end is None or last_end < endsl[-1]:
                            last_end = endsl[-1]
                        if not pd.isna(mean_values).all():
                            mean_value = mean_values[0]
                            if color_by_heatmap == 1:
                                cmap = colormap
                                norm = my_norm
                                color = cmap(norm(mean_value))
                            elif mean_values > colcutoff:
                                color = woodscolpos
                            elif mean_values < -colcutoff:
                                color = woodscolneg
                            else: 
                                color = woodscolneu
                            plt.plot([starts, ends], [mean_value, mean_value], label=f"{peptide}", color=color, lw=woodscolthick, path_effects=[pe.Stroke(linewidth=woodscolthick+woodsbackthick, foreground='black'), pe.Normal()])
                        else:
                            print(f"Skipping plotting for peptide {peptide}: data length mismatch.")
            if editing_figtitle == 1:
                if uploaded_settings != 1:
                    newtitle = str(request.form.get('alttitleword','Invalid Text'))
                plt.title(newtitle, fontsize=font_size_title, color=fontcolor, pad=wpadthick)
            else:
                if expType != 3:
                    plt.title(f"{title} - State Comparison: {state}, Exposure: {exposure}", fontsize=font_size_title, color=fontcolor, pad=wpadthick)
                else:
                    plt.title(f"{title} - State: {state}, Exposure: {exposure}", fontsize=font_size_title, color=fontcolor, pad=wpadthick)
            plt.xlabel(new_x_axis_title)
            plt.ylabel(new_y_axis_title)
            ax = plt.gca()
            ax.tick_params(axis='both', labelsize=font_size_ticklabel)
            if first_start-(last_end*1/20) <= 0:
                ax.set_xlim([0, last_end+(last_end*1/20)])
            else:
                ax.set_xlim([first_start-(last_end*1/20), last_end+(last_end*1/20)])
            if woodsplotYbyglobal == 1:
                ax.set_ylim([-max_delta-(max_delta/10), max_delta+(max_delta/10)])
            else:
                if subset.size > 0: 
                    max_y_woods = np.nanmax(np.abs(subset.values))  
                    if np.isnan(max_y_woods):  
                        max_y_woods = 2.5
                else:
                    max_y_woods = 2.5
                ax.set_ylim([-max_y_woods-(max_y_woods/10), max_y_woods+(max_y_woods/10)])
            for section in color_sections:
                    start, end, color, text = section
                    rect_y = blockposwoods
                    rect = plt.Rectangle((start+1, rect_y), end-start, blockthickwoods, color=color, clip_on=False)
                    ax.add_patch(rect)
                    text_y = blocktextposwoods + rect_y + blockthickwoods / 2 
                    ax.annotate(text, xy=((start+end)/2, text_y), xytext=(0,0), textcoords='offset points',
                                ha='center', va='center', fontsize=20, color=fontcolor, clip_on=False, annotation_clip=False)
                    print('adding section')
            plt.grid(True)
            #pdf.savefig(bbox_inches='tight')
            img_buffer = BytesIO()
            if PDFgeneration == 1:
                plt.savefig(img_buffer, format='pdf', dpi=output_bitmap_dpi, bbox_inches='tight')
            else:
                plt.savefig(img_buffer, format='png', dpi=output_bitmap_dpi, bbox_inches='tight')
            img_buffer.seek(0)
            if PDFgeneration == 1:
                zip_file.writestr(f"{title}_State_{state}_Exposure_{exposure}.pdf", img_buffer.read())
            else:
                zip_file.writestr(f"{title}_State_{state}_Exposure_{exposure}.png", img_buffer.read())
            plt.close()


    ##################################################################################
    ###########################      LADDER PLOTTING      ############################
    ##################################################################################


    if uploaded_settings != 1:
        altwoodsthick = float(request.form.get('altwoodsthick',0))
        if altwoodsthick == 1:
            woodscolthick = float(request.form.get('woodscolthick',3))
            woodsbackthick = float(request.form.get('woodsbackthick',0.5))
        woodsplotYbyglobal = float(request.form.get('woodsplotYbyglobal',1))
        if domainlabel == 1:
            wpadthick = wpadthick + 35
        blockposwoods = max_delta + max_delta/10
    elif 'blockposwoods' not in variables:
        blockposwoods = max_delta + max_delta/10 #for w_plot
    ladderplotcol = whitec
    flipYLadder = float(request.form.get('flipYLadder', 0))
    
    def l_plot(pdf, data, title, pept, state, data_source, max_delta, colcutoff):
        global output_bitmap_h_count
        if isinstance(data, np.ma.MaskedArray):
            data = data.filled(0)
        unique_exposures = data.columns.get_level_values('Exposure').unique()
        numplots = len(unique_exposures)
        data_avgrot = data_source.transpose()
        first_start = None
        last_end = None
        for exposure in unique_exposures:
            plt.figure(figsize=(woodsx, woodsy))
            if expType != 3:
                subset = data.xs(('Delta_d_uptake', 'mean', state, exposure), level=[None, 'Parameter', 'State', 'Exposure'], axis=1) 
            else:
                subset = data.xs(('d_uptake', 'mean', state, exposure), level=[None, 'Parameter', 'State', 'Exposure'], axis=1) 
            if colcutopt == 0:
                colcutoff = 0.5
            #print('Values for pept:', pept)
            #print('Index for pept:', pept.index)
            if PerResidueMap != 1: 
                peptide_count = 1
                for peptide in pept.index.get_level_values('Sequence').unique():
                    if peptide in subset.index:
                        mean_values = subset.loc[peptide].values
                        start_end_values = pept.loc[peptide, ['Start', 'End']].values
                        starts = start_end_values[0]
                        ends = start_end_values[1]
                        startsl = []
                        endsl = []
                        startsl.append(start_end_values[0])
                        endsl.append(start_end_values[1])
                        if first_start is None:
                            first_start = startsl[0]
                        if last_end is None or last_end < endsl[-1]:
                            last_end = endsl[-1]
                        # Plot the line for current peptide
                        if color_by_heatmap == 1:
                            if not pd.isna(mean_values).all():
                                mean_value = mean_values[0]
                                cmap = colormap
                                norm = my_norm
                                color = cmap(norm(mean_value))
                            else:
                                color = greyc
                                #print(f"Missing Mean data for peptide {peptide}")
                        else: 
                            color = ladderplotcol
                        plt.plot([starts, ends], [peptide_count, peptide_count], label=f"{peptide}", color=color, lw=woodscolthick, path_effects=[pe.Stroke(linewidth=woodscolthick+woodsbackthick, foreground='black'), pe.Normal()])
                        peptide_count += 1
            if editing_figtitle == 1:
                if uploaded_settings != 1:
                    newtitle = str(request.form.get('alttitleword','Invalid Text'))
                plt.title(newtitle, fontsize=font_size_title, color=fontcolor, pad=wpadthick)
            else:
                plt.title(f"{title} - State Comparison: {state}, Exposure: {exposure}", fontsize=font_size_title, color=fontcolor, pad=wpadthick)
            plt.xlabel(new_x_axis_title)
            plt.ylabel(new_y_axis_title)
            ax = plt.gca()
            ax.tick_params(axis='both', labelsize=font_size_ticklabel)
            if first_start-(last_end*1/20) <= 0:
                ax.set_xlim([0, last_end+(last_end*1/20)])
            else:
                ax.set_xlim([first_start-(last_end*1/20), last_end+(last_end*1/20)])
            ax.set_ylim([0, peptide_count+peptide_count/10])
            if flipYLadder == 1:
                ax.yaxis.set_inverted(True)
            for section in color_sections:
                    start, end, color, text = section
                    rect_y = blockposwoods
                    rect = plt.Rectangle((start+1, rect_y), end-start, blockthickwoods, color=color, clip_on=False)
                    ax.add_patch(rect)
                    text_y = blocktextposwoods + rect_y + blockthickwoods / 2 
                    ax.annotate(text, xy=((start+end)/2, text_y), xytext=(0,0), textcoords='offset points',
                                ha='center', va='center', fontsize=20, color=fontcolor, clip_on=False, annotation_clip=False)
                    print('adding section')
            plt.grid(True)
            #pdf.savefig(bbox_inches='tight')
            img_buffer = BytesIO()
            if PDFgeneration == 1:
                plt.savefig(img_buffer, format='pdf', dpi=output_bitmap_dpi, bbox_inches='tight')
            else:
                plt.savefig(img_buffer, format='png', dpi=output_bitmap_dpi, bbox_inches='tight')
            img_buffer.seek(0)
            if PDFgeneration == 1:
                zip_file.writestr(f"{title}_State_{state}_Exposure_{exposure}.pdf", img_buffer.read())
            else:
                zip_file.writestr(f"{title}_State_{state}_Exposure_{exposure}.png", img_buffer.read())
            plt.close()
        
    
    #################################################################################################################################################################
    #################################################################################################################################################################
    ###############################################################            WHAT TO MAKE            ##############################################################
    #################################################################################################################################################################
    #################################################################################################################################################################

    
    if h_or_v == 1:
        plot_v = 0
        plot_h = 1
        plot_w = 0
        plot_volc = 0
        plot_l = 0
    elif h_or_v == 2:
        plot_v = 1
        plot_h = 0
        plot_w = 0
        plot_volc = 0
        plot_l = 0
    elif h_or_v == 3:
        plot_v = 0
        plot_h = 0
        plot_w = 1
        plot_volc = 0
        plot_l = 0
    elif h_or_v == 4:
        plot_v = 0
        plot_h = 0
        plot_w = 0
        plot_volc = 1
        plot_l = 0
    elif h_or_v == 6:
        plot_v = 0
        plot_h = 0
        plot_w = 0
        plot_volc = 0
        plot_l = 1
    if plot_volc == 1:
        scatter_plot = 1
    elif plot_volc == 0:
        scatter_plot = 0

    heatmapPeptCutoffUse = float(request.form.get('heatmapPeptCutoffUse', 0))
    if heatmapPeptCutoffUse == 1:
        heatmapPeptCutoff = float(request.form.get('heatmapPeptCutoff', 0))
    else:
        heatmapPeptCutoff = 0
    heatmapPeptCutoff = int(heatmapPeptCutoff)
    
    download_store_values = int(request.form.get('store_values', 0))
    if download_store_values != 1:
        if plot_w == 1:
            all_deltas_grouped = all_deltas.groupby('Protein')
            output_buffer = BytesIO()
            with zipfile.ZipFile(output_buffer, 'w', zipfile.ZIP_DEFLATED, False) as zip_file:
                for prot, data_subset in all_deltas_grouped:
                    for state in all_states:
                        if expType != 3:
                            data = data_subset.loc[:,('Delta_d_uptake','mean',state,slice(None))]
                        else:
                            data = data_subset.loc[:,('d_uptake','mean',state,slice(None))]
                        pept = data_subset.loc[:,['Start','End']]
                        title = prot
                        w_plot(pdf,data,title,pept,state,data_source,max_delta,colcutoff)
            output_buffer.seek(0)
            #zip_filename = generate_unique_filename('Plots', 'zip', timestamp)
            #upload_to_blob_storage(container_client, heatmap_buffer.getvalue(), zip_filename)
        elif plot_l == 1:
            all_deltas_grouped = all_deltas.groupby('Protein')
            output_buffer = BytesIO()
            with zipfile.ZipFile(output_buffer, 'w', zipfile.ZIP_DEFLATED, False) as zip_file:
                for prot, data_subset in all_deltas_grouped:
                    for state in all_states:
                        if expType != 3:
                            data = data_subset.loc[:,('Delta_d_uptake','mean',state,slice(None))]
                        else:
                            data = data_subset.loc[:,('d_uptake','mean',state,slice(None))]
                        pept = data_subset.loc[:,['Start','End']]
                        title = prot
                        l_plot(pdf,data,title,pept,state,data_source,max_delta,colcutoff)
            output_buffer.seek(0)
        elif split_outp_by_prot != None:    
            if split_outp_by_prot == 'all':
                all_deltas_grouped = all_deltas.groupby('Protein')
                for prot, data_subset in all_deltas_grouped:
                    if plot_stacked != None:
                        if plot_stacked:
                            if expType != 3:
                                data = data_subset.loc[:,('Delta_d_uptake','mean',slice(None),slice(None))]
                            else:
                                data = data_subset.loc[:,('d_uptake','mean',slice(None),slice(None))]
                            pept = data_subset.loc[:,['Start','End']]
                            if plot_h == 1:
                                if expType != 3:
                                    title = 'D uptake difference; Protein: ' + prot + '\nTop to Bottom: ' + ', '.join(all_states.to_list())
                                elif absoluteUptakeValues != 1:
                                    if percentUseRFU == 100:
                                        title = 'Percent Relative D uptake; Protein: ' + prot + '\nTop to Bottom: ' + ', '.join(all_states.to_list())
                                    else:
                                        title = 'Relative D uptake; Protein: ' + prot + '\nTop to Bottom: ' + ', '.join(all_states.to_list())
                                else:
                                    title = 'Absolute D uptake; Protein: ' + prot + '\nTop to Bottom: ' + ', '.join(all_states.to_list())
                                if heatmapPeptCutoff != 0 and len(pept) > heatmapPeptCutoff:
                                    chunked_data = np.array_split(data, np.ceil(len(data) / heatmapPeptCutoff))
                                    chunked_pept = np.array_split(pept, np.ceil(len(pept) / heatmapPeptCutoff))
                                    for i, (chunk, pept_chunk) in enumerate(zip(chunked_data, chunked_pept)):
                                        h_plot(pdf, chunk, f"{title} (Part {i+1})", pept_chunk, color_sections)
                                else:
                                    h_plot(pdf, data, title, pept, color_sections)
                            if plot_v == 1:
                                if expType != 3:
                                    title = 'D uptake difference; Protein: ' + prot + '\nLeft to Right: ' + ', '.join(all_states.to_list())
                                elif absoluteUptakeValues != 1:
                                    if percentUseRFU == 100:
                                        title = 'Percent Relative D uptake; Protein: ' + prot + '\nLeft to Right: ' + ', '.join(all_states.to_list())
                                    else:
                                        title = 'Relative D uptake; Protein: ' + prot + '\nLeft to Right: ' + ', '.join(all_states.to_list())
                                else:
                                    title = 'Absolute D uptake; Protein: ' + prot + '\nLeft to Right: ' + ', '.join(all_states.to_list())
                                if heatmapPeptCutoff != 0 and len(pept) > heatmapPeptCutoff:
                                    chunked_data = np.array_split(data, np.ceil(len(data) / heatmapPeptCutoff))
                                    chunked_pept = np.array_split(pept, np.ceil(len(pept) / heatmapPeptCutoff))
                                    for i, (chunk, pept_chunk) in enumerate(zip(chunked_data, chunked_pept)):
                                        v_plot(pdf, chunk, f"{title} (Part {i+1})", pept_chunk)
                                else:
                                    v_plot(pdf, data, title, pept)
                    if plot_separate != None:
                        if plot_separate:
                            for state in all_states:
                                if expType != 3:
                                    data = data_subset.loc[:,('Delta_d_uptake','mean',state,slice(None))]
                                else: 
                                    data = data_subset.loc[:,('d_uptake','mean',state,slice(None))]
                                pept = data_subset.loc[:,['Start','End']]
                                if expType != 3:
                                    title = 'D uptake difference: ' + state + '\nProtein: ' + prot
                                elif absoluteUptakeValues != 1:
                                    if percentUseRFU == 100:
                                        title = 'Percent Relative D uptake: ' + state + '\nProtein: ' + prot
                                    else:
                                        title = 'Relative D uptake: ' + state + '\nProtein: ' + prot
                                else:
                                    title = 'Absolute D uptake: ' + state + '\nProtein: ' + prot
                                if plot_h == 1:
                                    if heatmapPeptCutoff != 0 and len(pept) > heatmapPeptCutoff:
                                        chunked_data = np.array_split(data, np.ceil(len(data) / heatmapPeptCutoff))
                                        chunked_pept = np.array_split(pept, np.ceil(len(pept) / heatmapPeptCutoff))
                                        for i, (chunk, pept_chunk) in enumerate(zip(chunked_data, chunked_pept)):
                                            h_plot(pdf, chunk, f"{title} (Part {i+1})", pept_chunk, color_sections)
                                    else:
                                        h_plot(pdf, data, title, pept, color_sections)
                                if plot_v == 1:
                                    if heatmapPeptCutoff != 0 and len(pept) > heatmapPeptCutoff:
                                        chunked_data = np.array_split(data, np.ceil(len(data) / heatmapPeptCutoff))
                                        chunked_pept = np.array_split(pept, np.ceil(len(pept) / heatmapPeptCutoff))
                                        for i, (chunk, pept_chunk) in enumerate(zip(chunked_data, chunked_pept)):
                                            v_plot(pdf, chunk, f"{title} (Part {i+1})", pept_chunk)
                                    else:
                                        v_plot(pdf, data, title, pept)
        elif type(split_outp_by_prot) == list and len(split_outp_by_prot) > 0:
            all_deltas_grouped = all_deltas.groupby('Protein')
            for prot in split_outp_by_prot:
                if type(prot) == tuple:
                    data_subset = pd.concat([all_deltas_grouped.get_group(key) for key in prot])
                    prot_string = '; Proteins: ' + ' '.join(prot)
                else:
                    data_subset = all_deltas_grouped.get_group(prot)
                    prot_string = '; Protein: ' + prot
                if plot_stacke == 1:
                    if expType != 3:
                        data = data_subset.loc[:,('Delta_d_uptake','mean',slice(None),slice(None))]
                    else:
                        data = data_subset.loc[:,('d_uptake','mean',slice(None),slice(None))]
                    pept = data_subset.loc[:,['Start','End']]
                    if plot_h == 1:
                        if expType != 3:
                            title = 'D uptake difference; Protein: ' + prot_string + '\nTop to bottom: ' + ', '.join(all_states.to_list())
                        elif absoluteUptakeValues != 1:
                            if percentUseRFU == 100:
                                title = 'Percent Relative D uptake; Protein: ' + prot_string + '\nTop to bottom: ' + ', '.join(all_states.to_list())
                            else:
                                title = 'Relative D uptake; Protein: ' + prot_string + '\nTop to bottom: ' + ', '.join(all_states.to_list())
                        else:
                            title = 'Absolute D uptake; Protein: ' + prot_string + '\nTop to bottom: ' + ', '.join(all_states.to_list())
                        if heatmapPeptCutoff != 0 and len(pept) > heatmapPeptCutoff:
                            chunked_data = np.array_split(data, np.ceil(len(data) / heatmapPeptCutoff))
                            chunked_pept = np.array_split(pept, np.ceil(len(pept) / heatmapPeptCutoff))
                            for i, (chunk, pept_chunk) in enumerate(zip(chunked_data, chunked_pept)):
                                h_plot(pdf, chunk, f"{title} (Part {i+1})", pept_chunk, color_sections)
                        else:
                            h_plot(pdf, data, title, pept, color_sections)
                    if plot_v == 1:
                        if expType != 3:
                            title = 'D uptake difference; Protein: ' + prot_string + '\nLeft to Right: ' + ', '.join(all_states.to_list())
                        elif absoluteUptakeValues != 1:
                            if percentUseRFU == 100:
                                title = 'Percent Relative D uptake; Protein: ' + prot_string + '\nLeft to Right: ' + ', '.join(all_states.to_list())
                            else:
                                title = 'Relative D uptake; Protein: ' + prot_string + '\nLeft to Right: ' + ', '.join(all_states.to_list())
                        else:
                            title = 'Absolute D uptake; Protein: ' + prot_string + '\nLeft to Right: ' + ', '.join(all_states.to_list())
                        if heatmapPeptCutoff != 0 and len(pept) > heatmapPeptCutoff:
                            chunked_data = np.array_split(data, np.ceil(len(data) / heatmapPeptCutoff))
                            chunked_pept = np.array_split(pept, np.ceil(len(pept) / heatmapPeptCutoff))
                            for i, (chunk, pept_chunk) in enumerate(zip(chunked_data, chunked_pept)):
                                v_plot(pdf, chunk, f"{title} (Part {i+1})", pept_chunk)
                        else:
                            v_plot(pdf, data, title, pept)
                elif plot_separate == 1:
                    for state in all_states:
                        if expType != 3:
                            data = data_subset.loc[:,('Delta_d_uptake','mean',state,slice(None))]
                        else:
                            data = data_subset.loc[:,('d_uptake','mean',state,slice(None))]
                        pept = data_subset.loc[:,['Start','End']]
                        if expType != 3:
                            title = 'D uptake difference: ' + state + prot_string
                        elif absoluteUptakeValues != 1:
                            if percentUseRFU == 100:
                                title = 'Percent Relative D uptake: ' + state + prot_string
                            else:
                                title = 'Relative D uptake: ' + state + prot_string
                        else:
                            title = 'Absolute D uptake: ' + state + prot_string
                        if plot_h == 1:
                            if heatmapPeptCutoff != 0 and len(pept) > heatmapPeptCutoff:
                                chunked_data = np.array_split(data, np.ceil(len(data) / heatmapPeptCutoff))
                                chunked_pept = np.array_split(pept, np.ceil(len(pept) / heatmapPeptCutoff))
                                for i, (chunk, pept_chunk) in enumerate(zip(chunked_data, chunked_pept)):
                                    h_plot(pdf, chunk, f"{title} (Part {i+1})", pept_chunk, color_sections)
                            else:
                                h_plot(pdf, data, title, pept, color_sections)
                        if plot_v == 1:
                            if heatmapPeptCutoff != 0 and len(pept) > heatmapPeptCutoff:
                                chunked_data = np.array_split(data, np.ceil(len(data) / heatmapPeptCutoff))
                                chunked_pept = np.array_split(pept, np.ceil(len(pept) / heatmapPeptCutoff))
                                for i, (chunk, pept_chunk) in enumerate(zip(chunked_data, chunked_pept)):
                                    v_plot(pdf, chunk, f"{title} (Part {i+1})", pept_chunk)
                            else:
                                v_plot(pdf, data, title, pept)
        else:
            if plot_stacked == 1:
                if expType != 3:
                    data = all_deltas.loc[:,('Delta_d_uptake','mean',slice(None),slice(None))]
                else:
                    data = all_deltas.loc[:,('d_uptake','mean',slice(None),slice(None))]
                pept = all_deltas.loc[:,['Start','End']]
                if plot_h == 1:
                    if expType != 3:
                        title = 'D uptake difference\nTop to bottom: ' + ', '.join(all_states.to_list())
                    elif absoluteUptakeValues != 1:
                        if percentUseRFU == 100:
                            title = 'Percent Relative D uptake \nTop to bottom: ' + ', '.join(all_states.to_list())
                        else:
                            title = 'Relative D uptake \nTop to bottom: ' + ', '.join(all_states.to_list())
                    else:
                        title = 'Absolute D uptake \nTop to bottom: ' + ', '.join(all_states.to_list())
                    if heatmapPeptCutoff != 0 and len(pept) > heatmapPeptCutoff:
                        chunked_data = np.array_split(data, np.ceil(len(data) / heatmapPeptCutoff))
                        chunked_pept = np.array_split(pept, np.ceil(len(pept) / heatmapPeptCutoff))
                        for i, (chunk, pept_chunk) in enumerate(zip(chunked_data, chunked_pept)):
                            h_plot(pdf, chunk, f"{title} (Part {i+1})", pept_chunk, color_sections)
                    else:
                        h_plot(pdf, data, title, pept, color_sections)
                if plot_v == 1:
                    if expType != 3:
                        title = 'D uptake difference\nLeft to right: ' + ', '.join(all_states.to_list())
                    elif absoluteUptakeValues != 1:
                        if percentUseRFU == 100:
                            title = 'Percent Relative D uptake \nLeft to right: ' + ', '.join(all_states.to_list())
                        else:
                            title = 'Relative D uptake \nLeft to right: ' + ', '.join(all_states.to_list())
                    else:
                        title = 'Absolute D uptake \nLeft to right: ' + ', '.join(all_states.to_list())
                    if heatmapPeptCutoff != 0 and len(pept) > heatmapPeptCutoff:
                        chunked_data = np.array_split(data, np.ceil(len(data) / heatmapPeptCutoff))
                        chunked_pept = np.array_split(pept, np.ceil(len(pept) / heatmapPeptCutoff))
                        for i, (chunk, pept_chunk) in enumerate(zip(chunked_data, chunked_pept)):
                            v_plot(pdf, chunk, f"{title} (Part {i+1})", pept_chunk)
                    else:
                        v_plot(pdf, data, title, pept)
            elif plot_separate == 1:
                for state in all_states:
                    if expType != 3:
                        data = all_deltas.loc[:,('Delta_d_uptake','mean',state,slice(None))]
                    else: 
                        data = all_deltas.loc[:,('d_uptake','mean',state,slice(None))]
                    if expType != 3:
                        title = 'D uptake difference: ' + state
                    elif absoluteUptakeValues != 1:
                        if percentUseRFU == 100:
                            title = 'Percent Relative D uptake: ' + state
                        else:
                            title = 'Relative D uptake: ' + state
                    else:
                        title = 'Absolute D uptake: ' + state
                    pept = all_deltas.loc[:,['Start','End']]
                    if plot_h == 1:
                        if heatmapPeptCutoff != 0 and len(pept) > heatmapPeptCutoff:
                            chunked_data = np.array_split(data, np.ceil(len(data) / heatmapPeptCutoff))
                            chunked_pept = np.array_split(pept, np.ceil(len(pept) / heatmapPeptCutoff))
                            for i, (chunk, pept_chunk) in enumerate(zip(chunked_data, chunked_pept)):
                                h_plot(pdf, chunk, f"{title} (Part {i+1})", pept_chunk, color_sections)
                        else:
                            h_plot(pdf, data, title, pept, color_sections)
                    if plot_v == 1:
                        if heatmapPeptCutoff != 0 and len(pept) > heatmapPeptCutoff:
                            chunked_data = np.array_split(data, np.ceil(len(data) / heatmapPeptCutoff))
                            chunked_pept = np.array_split(pept, np.ceil(len(pept) / heatmapPeptCutoff))
                            for i, (chunk, pept_chunk) in enumerate(zip(chunked_data, chunked_pept)):
                                v_plot(pdf, chunk, f"{title} (Part {i+1})", pept_chunk)
                        else:
                            v_plot(pdf, data, title, pept)
        #pdf.close()


    #################################################################################################################################################################
    #################################################################################################################################################################
    ############################################################0##            WHAT TO SEND            ##############################################################
    #################################################################################################################################################################
    #################################################################################################################################################################

    
    download_pymol = int(request.form.get('download_pymol', 0))
    download_chimera = int(request.form.get('download_chimera', 0))
    output_files = []
    if uploaded_settings != 1:
        h_or_v = float(request.form['h_or_v'])
    if h_or_v == 1:
        plot_v = 0
        plot_h = 1
        plot_w = 0
        plot_volc = 0
        plot_l = 0
    elif h_or_v == 2:
        plot_v = 1
        plot_h = 0
        plot_w = 0
        plot_volc = 0
        plot_l = 0
    elif h_or_v == 3:
        plot_v = 0
        plot_h = 0
        plot_w = 1
        plot_volc = 0
        plot_l = 0
    elif h_or_v == 4:
        plot_v = 0
        plot_h = 0
        plot_w = 0
        plot_volc = 1
        plot_l = 0
    elif h_or_v == 6:
        plot_v = 0
        plot_h = 0
        plot_w = 0
        plot_volc = 0
        plot_l = 1
    if plot_volc == 1:
        scatter_plot = 1
    elif plot_volc == 0:
        scatter_plot = 0
    
    def download_collect(chain_dict, chain_id, colors, bounds, all_deltas, all_states, buffer, heatmap_buffer, output_buffer, mk_pymol):
        h_or_v = float(request.form['h_or_v'])
        if h_or_v == 1:
            plot_v = 0
            plot_h = 1
            plot_w = 0
            plot_volc = 0
            plot_l = 0
        elif h_or_v == 2:
            plot_v = 1
            plot_h = 0
            plot_w = 0
            plot_volc = 0
            plot_l = 0
        elif h_or_v == 3:
            plot_v = 0
            plot_h = 0
            plot_w = 1
            plot_volc = 0
            plot_l = 0
        elif h_or_v == 4:
            plot_v = 0
            plot_h = 0
            plot_w = 0
            plot_volc = 1
            plot_l = 0
        elif h_or_v == 6:
            plot_v = 0
            plot_h = 0
            plot_w = 0
            plot_volc = 0
            plot_l = 1
        if plot_volc == 1:
            scatter_plot = 1
        elif plot_volc == 0:
            scatter_plot = 0
        pymol_dir = 'pymol_macros'
        Path(pymol_dir).mkdir(parents=True, exist_ok=True)
        # Populate chain_dict if not provided
        if not chain_dict:
            for prot in all_deltas['Protein'].unique():
                chain_dict[prot] = chain_id
        # Initialize PyMOL header - THIS DOESN'T LOOK GOOD FOR COLOURING MULTIPLE CHAINS, JUST AVOID IT LET PEOPLE COLOUR BASE WHATEVER
        #pymol_header = "alter all, b=1000\ncolor %s, all\n\n" % ('0x' + mc.to_hex(c_missing)[1:])
        #pymol_header = ''
        for prot, data_subset in all_deltas_grouped:
            protein_data = data_subset[data_subset['Protein'] == prot]
            chains = chain_dict.get(prot, '').split(',')
            chain_selection = '+'.join(chains)
            pymol_header = (
                f"alter chain {chain_selection}, b=1000\n"
                f"color %s, chain {chain_selection}\n\n" % ('0x' + mc.to_hex(c_missing)[1:])
            )
            pymol_footer = f"\n\ncolor %s, chain {chain_selection} and b<%.3f\n" % ('0x' + mc.to_hex(colors[0])[1:], bounds[0])
            for i in range(len(bounds) - 1):
                curr_color = '0x' + mc.to_hex(colors[i])[1:]
                pymol_footer += f"color %s, chain {chain_selection} and (b>%.3f or b=%.3f) and b<%.3f\n" % (
                    curr_color, bounds[i], bounds[i], bounds[i + 1]
                )
            curr_color = '0x' + mc.to_hex(colors[-1])[1:]
            pymol_footer += f"color %s, chain {chain_selection} and (b>%.3f or b=%.3f) and b<999\n" % (
                curr_color, bounds[-1], bounds[-1]
            )
            for state in all_states:
                if expType != 3:
                    data = protein_data.loc[:, ('Delta_d_uptake', 'mean', state, slice(None))]
                else:
                    data = protein_data.loc[:, ('d_uptake', 'mean', state, slice(None))]
                time_list = data.columns.get_level_values('Exposure').unique()
                for i in time_list:
                    if expType != 3:
                        mask = data.loc[:, ('Delta_d_uptake', 'mean', state, i)].notna().to_list()
                        pymol_str = pd.Series(
                            list(map(mk_pymol, protein_data['Protein'], protein_data['Start'], protein_data['End'], data.loc[:, ('Delta_d_uptake', 'mean', state, i)]))
                        )[mask]
                    else:
                        mask = data.loc[:, ('d_uptake', 'mean', state, i)].notna().to_list()
                        pymol_str = pd.Series(
                            list(map(mk_pymol, protein_data['Protein'], protein_data['Start'], protein_data['End'], data.loc[:, ('d_uptake', 'mean', state, i)]))
                        )[mask]
                    print('Pymol Str:', pymol_str)
                    pymol_script = pymol_str.to_string(header=False, index=False)
                    print('Pymol Script:', pymol_script)
                    file_name = "%s_%s_%s.pml" % (prot, state, str(i))
                    if pymol_dir:
                        file_name = os.path.join(pymol_dir, file_name)
                    with open(file_name, 'w') as f_pymol:
                        f_pymol.write(pymol_header + '\n')
                        for command in pymol_str:
                            f_pymol.write(command + '\n')
                        f_pymol.write(pymol_footer + '\n')
                    output_files.append(file_name)
        zip_file_name = 'output_files.zip'
        with zipfile.ZipFile(zip_file_name, 'w') as zipf:
            for file in output_files:
                zipf.write(file, os.path.basename(file))
        combined_buffer = io.BytesIO()
        with zipfile.ZipFile(combined_buffer, 'w') as zipf:
            zipf.write(zip_file_name, os.path.basename(zip_file_name))
            if PerResidueMap != 1:
                if h_or_v == 3:
                    zipf.writestr('WoodsPlots.zip', output_buffer.getbuffer())
                elif len(all_prot) > 1:
                    zipf.writestr('heatmaps.zip', heatmap_buffer.getbuffer())
                elif PDFgeneration == 1:
                    zipf.writestr('heatmap.zip', heatmap_buffer.getbuffer())
                else:
                    zipf.writestr('heatmap.zip', heatmap_buffer.getbuffer())
        combined_buffer.seek(0)
        plot_and_pymol_name = generate_unique_filename('Plot_and_Pymol', 'zip', timestamp)
        upload_pymol = upload_to_blob_storage(container_client, combined_buffer.getvalue(), plot_and_pymol_name)
        zip_blob_stream = download_blob_as_bytes(container_client, plot_and_pymol_name)
        response = make_response(send_file(zip_blob_stream, as_attachment=True, download_name=plot_and_pymol_name))
        response.headers['Content-Disposition'] = f'attachment; filename={plot_and_pymol_name}'
        response.headers['Content-Type'] = 'application/zip'
        return response
        #return send_file(zip_blob_stream, as_attachment=True, download_name=plot_and_pymol_name, mimetype='application/zip')

    def download_collect_chimerax(chain_dict, chain_id, colors, bounds, all_deltas, all_states, buffer, heatmap_buffer, output_buffer, mk_chimerax):
        h_or_v = float(request.form['h_or_v'])
        if h_or_v == 1:
            plot_v = 0
            plot_h = 1
            plot_w = 0
            plot_volc = 0
            plot_l = 0
        elif h_or_v == 2:
            plot_v = 1
            plot_h = 0
            plot_w = 0
            plot_volc = 0
            plot_l = 0
        elif h_or_v == 3:
            plot_v = 0
            plot_h = 0
            plot_w = 1
            plot_volc = 0
            plot_l = 0
        elif h_or_v == 4:
            plot_v = 0
            plot_h = 0
            plot_w = 0
            plot_volc = 1
            plot_l = 0
        elif h_or_v == 6:
            plot_v = 0
            plot_h = 0
            plot_w = 0
            plot_volc = 0
            plot_l = 1
        if plot_volc == 1:
            scatter_plot = 1
        elif plot_volc == 0:
            scatter_plot = 0
        chimerax_dir = 'chimerax_macros'
        Path(chimerax_dir).mkdir(parents=True, exist_ok=True)
        if chain_dict:
            if chain_id:
                print('Both variables chain_dict and chain_id are not empty; Using values from chain_dict, which overrides chain_id.\n')
        else:
            for prot in all_deltas['Protein'].unique():
                chain_dict[prot] = chain_id
        chimerax_header = "# ChimeraX coloring script\n\n"
        usable_bounds = bounds[1:-1]
        usable_colors = colors[:-1] 
        palette_definition = '"'
        for i in range(len(usable_bounds)):
            palette_definition += f"{usable_bounds[i]:.3f},{to_hex(usable_colors[i])}:"
            if i < len(usable_bounds) - 1:
                palette_definition += f"{usable_bounds[i]:.3f},{to_hex(usable_colors[i + 1])}:"
        palette_definition += f"{usable_bounds[-1]:.3f},{to_hex(colors[-1])}"
        palette_definition = palette_definition.rstrip(':')
        palette_definition += '"'
        chimerax_footer = f"\n\ncolor byattribute colorval palette {palette_definition}\n"
        if chaindictuse == 1:
            chain_to_col = {key: '/' + ','.join(value if isinstance(value, list) else [value])
                    for key, value in chain_dict.items()}
        else: 
            chain_to_col = {}
        output_files = []
        for prot, data_subset in all_deltas.groupby('Protein'):  
            chain_tocol = chain_to_col.get(prot, '')
            print('Chain to col:', chain_tocol)
            for state in all_states:
                if expType != 3:
                    data = all_deltas.loc[:, ('Delta_d_uptake', 'mean', state, slice(None))]
                else:
                    data = all_deltas.loc[:, ('d_uptake', 'mean', state, slice(None))]
                time_list = data.columns.get_level_values('Exposure').unique()
                for i, time in enumerate(time_list):
                    chimerax_commands = []
                    if expType != 3:
                        mask = data.loc[:, ('Delta_d_uptake', 'mean', state, time)].notna()
                        chimerax_str2 = pd.Series(
                            list(map(
                                lambda start, end, bfactor: mk_chimerax(chain_tocol, start, end, bfactor),
                                data_subset.loc[mask, 'Start'],
                                data_subset.loc[mask, 'End'],
                                data.loc[mask, ('Delta_d_uptake', 'mean', state, time)]
                            ))
                        )
                    else:
                        mask = data.loc[:, ('d_uptake', 'mean', state, time)].notna()
                        chimerax_str2 = pd.Series(
                            list(map(
                                lambda start, end, bfactor: mk_chimerax(chain_tocol, start, end, bfactor),
                                data_subset.loc[mask, 'Start'],
                                data_subset.loc[mask, 'End'],
                                data.loc[mask, ('d_uptake', 'mean', state, time)]
                            ))
                        )
                    chimerax_commands.extend(chimerax_str2.tolist())
                    chimerax_script = "\n".join(chimerax_commands)
                    file_name = f"{prot.replace(' ', '_')}_{state.replace(' ', '_')}_{str(time).replace(' ', '_')}.cxc"
                    if chimerax_dir:
                        file_name = os.path.join(chimerax_dir, file_name)
                    with open(file_name, 'w') as f_chimerax:
                        print(chimerax_header, chimerax_script, chimerax_footer, file=f_chimerax)
                    output_files.append(file_name)
        zip_file_name = 'output_files.zip'
        with zipfile.ZipFile(zip_file_name, 'w') as zipf:
            for file in output_files:
                zipf.write(file, os.path.basename(file))
        combined_buffer = io.BytesIO()
        with zipfile.ZipFile(combined_buffer, 'w') as zipf:
            zipf.write(zip_file_name, os.path.basename(zip_file_name))
            if PerResidueMap != 1:
                if h_or_v == 3:
                    zipf.writestr('WoodsPlots.zip', output_buffer.getbuffer())
                elif len(all_deltas['Protein'].unique()) > 1:
                    zipf.writestr('heatmaps.zip', heatmap_buffer.getbuffer())
                else:
                    zipf.writestr('heatmap.zip', heatmap_buffer.getbuffer())
        combined_buffer.seek(0)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        zip_filename =  generate_unique_filename('Plot_and_ChimeraX', 'zip', timestamp)
        response = make_response(send_file(combined_buffer, as_attachment=True, download_name=zip_filename))
        response.headers['Content-Disposition'] = f'attachment; filename={zip_filename}'
        response.headers['Content-Type'] = 'application/zip'
        return response

    if download_store_values == 1:
        AZURE_STORAGE_CONNECTION_STRINGsv = "DefaultEndpointsProtocol=https;AccountName=heatmap1vahidi;AccountKey=4R083EcQh6qmzDyTxZ0MChBPJ3hWYFz3f7Y/Zf3SyHEyach77fV/Sl0YSMZ/6B8bbdXHsUzNjfKq+AStG05iCQ==;EndpointSuffix=core.windows.net"
        CONTAINER_NAMEsv = "heatmap1vahidi"
        timestamp = datetime.now().strftime("%Y.%m.%d_%H.%M")
        BLOB_NAMEsv = f'HDgraphiX_settings_{timestamp}.csv'
        blob_service_clientsv = BlobServiceClient.from_connection_string(AZURE_STORAGE_CONNECTION_STRINGsv)
        container_clientsv = blob_service_client.get_container_client(CONTAINER_NAMEsv)
        redc = '#c50f15'
        dredc = '#990000'
        bluec = '#0062cc'
        dbluec = '#08306b'
        whitec = '#FFFFFF'
        blackc = '#000000'
        greyc = '#E0E0E0'
        greenc = '#006600'
        yellowc = '#e69b00'
        purplec = '#643d6e'
        c1 = bluec    #blue 
        c2 = whitec   #white
        c3 = redc     #red
        c_missing = '#bdbdbd' # gray;
        fontchoice = float(request.form.get('fontchoiceu', 1))
        PDFgeneration = float(request.form.get('gen_pdf',0))
        download_pymol = int(request.form.get('download_pymol', 0))
        download_chimera = int(request.form.get('download_chimera', 0))
        alt_c2 = float(request.form.get('alt_c2',0))
        if alt_c2 == 1:
            c2 = ensure_hex_color(request.form.get('alt_c2_col', '#FFFFFF'))
        else:
            c2 = whitec
        alt_cmissing = float(request.form.get('alt_cmissing',0))
        if alt_cmissing == 1:
            c_missing = ensure_hex_color(request.form.get('alt_cmissing_col','#bdbdbd'))
        else:
            c_missing = '#bdbdbd'
        color_by_heatmap = float(request.form.get('colorbyheatmap',0))
        colcutopt = float(request.form.get('colcutopt',0))
        if colcutopt == 1:
            colcutoff = float(request.form.get('colcutoff',0.5))
        else:
            colcutoff = 0.5
        woodsdimen = float(request.form.get('woodsdimen',0))
        if woodsdimen == 1:
            woodsx = float(request.form.get('woodsx',20))
            woodsy = float(request.form.get('woodsy',6))
        else:
            woodsx = 20
            woodsy = 6
        woodscol = float(request.form.get('woodscol',0))
        if woodscol == 1:
            woodscolpos = request.form.get('woodscolpos','#FF0000')
            woodscolneu = request.form.get('woodscolneu',whitec)
            woodscolneg = request.form.get('woodscolneg','#0000FF')
            woodscolpos = ensure_hex_color(woodscolpos)
            woodscolneu = ensure_hex_color(woodscolneu)
            woodscolneg = ensure_hex_color(woodscolneg)
        else:
            woodscolpos = '#FF0000'
            woodscolneu = whitec
            woodscolneg = '#0000FF'
        nolines = float(request.form.get('nolines', 0))
        font_size = float(request.form['font_size'])
        font_size_title = float(request.form['font_size_title'])
        font_size_ticklabel = float(request.form['font_size_ticklabel'])
        fontcolor = blackc
        fontcolorc = float(request.form.get('fontcolor',2))
        if fontcolorc == 1: #White
            fontcolor = whitec 
        elif fontcolorc == 2: #Black
            fontcolor = blackc
        elif fontcolorc == 3: #Blue
            fontcolor = dbluec
        elif fontcolorc == 4: #Red
            fontcolor = dredc
        elif fontcolorc == 5: #Grey
            fontcolor = greyc
        zerobound = float(request.form.get('zerobound', 0))
        funkybound = float(request.form.get('change_bounds_abs',3))
        globalmax_delta = float(request.form.get('global_max', 0))
        #relativeUptakeCalc = float(request.form.get('relativeUptakeCalc', 0))
        relativeUptakeCalc = 0
        p_threshold = request.form.get('pthresh', 0.05)
        if p_threshold == '':
            p_threshold = 0.05
        if p_threshold != 0.05 and p_threshold != '':
            p_threshold = float(p_threshold)
        custb = float(request.form['option'])
        if custb == 3:
            max_range = float(request.form['max_range']) 
            max_range = round(max_range, 1)
            num_shades = float(request.form['num_shades'])
            num_shades = round(num_shades) + 1
            alt_col = float(request.form.get('alt_col', 0))
            if alt_col == 1:
                if max_range <= 2:
                    c1 = request.form.get('negcolalt',bluec).strip()
                    c1 = ensure_hex_color(c1)
                    c3 = request.form.get('poscolalt',redc).strip()
                    c3 = ensure_hex_color(c3)
                else:
                    c1 = request.form.get('negcolalt',dbluec).strip()
                    c1 = ensure_hex_color(c1)
                    c3 = request.form.get('poscolalt',dredc).strip()
                    c3 = ensure_hex_color(c3)
        elif custb == 2: 
            num_shades = 7
            alt_col = float(request.form.get('alt_col', 0))
            if alt_col == 0:
                if funkybound == 1:
                    custom_colors = [c2, '#FC9272', '#FB6A4A', '#EF3B2C', '#CB181D', '#A50F15', '#67000D']
                elif funkybound == 2:
                    custom_colors = ['#023858', '#045A8D', '#0570B0', '#3690C0', '#74A9CF', '#A6BDDB', c2]
                else:
                    if zerobound == 0:
                        custom_colors = ['#023858', '#045A8D', '#0570B0', '#3690C0', '#74A9CF', '#A6BDDB', c2, '#FC9272', '#FB6A4A', '#EF3B2C', '#CB181D', '#A50F15', '#67000D']
                    elif zerobound == 1:
                        custom_colors = ['#023858', '#045A8D', '#0570B0', '#3690C0', '#74A9CF', '#A6BDDB', '#FC9272', '#FB6A4A', '#EF3B2C', '#CB181D', '#A50F15', '#67000D']
            elif alt_col == 1:
                ncolforalt = 7
                if funkybound == 1:
                    c3 = request.form.get('poscolalt', dredc).strip()
                    c3 = ensure_hex_color(c3)
                    clist = [colorFader(c2,c3,x/ncolforalt) for x in range(ncolforalt)]
                    custom_colors = [c2] + clist[1:]
                elif funkybound == 2:
                    c1 = request.form.get('negcolalt', dbluec).strip()
                    c1 = ensure_hex_color(c1)
                    ncolforaltneg = ncolforalt - 1
                    clist = [colorFader(c1,c2,x/ncolforaltneg) for x in range(ncolforaltneg)]
                    custom_colors = clist + [c2]
                else:
                    c1 = request.form.get('negcolalt', dbluec).strip()
                    c3 = request.form.get('poscolalt', dredc).strip()
                    c1 = ensure_hex_color(c1)
                    c3 = ensure_hex_color(c3)
                    ncolforaltneg = ncolforalt - 1
                    positive_c = [colorFader(c2,c3,x/ncolforalt) for x in range(ncolforalt)]
                    print('positive_c:', positive_c)
                    negative_c = [colorFader(c1,c2,x/ncolforaltneg) for x in range(ncolforaltneg)]
                    print('negative_c:', negative_c)
                    custom_colors = []
                    if zerobound == 0:
                        custom_colors.extend(negative_c)
                        custom_colors.extend(positive_c)
                    elif zerobound == 1:
                        custom_colors.extend(negative_c)
                        custom_colors.extend(positive_c)
                        while '#ffffff' in custom_colors:
                            custom_colors.remove('#ffffff')
                        while '#FFFFFF' in custom_colors:
                            custom_colors.remove('#FFFFFF')
                print('custom colours, custb 2:', custom_colors)
                num_shades = 8
        elif custb == 1: 
            max_range = None
            custom_bounds = []
            negative_bounds = []
            positive_bounds = []
            for i in range(1, 9):
                nb = request.form.get(f'inputbn{i}', '').strip()
                pb = request.form.get(f'inputbp{i}', '').strip() 
                if nb:
                    try:
                        nb = float(nb)
                        negative_bounds.append(nb)
                    except ValueError:
                        pass
                if pb:
                    try:
                        pb = float(pb)
                        positive_bounds.append(pb) 
                    except ValueError:
                        pass 
            if zerobound == 1:
                custom_bounds.append(0)
            custom_bounds.extend(negative_bounds)
            custom_bounds.extend(positive_bounds)
            custom_bounds = sorted(custom_bounds)
            numinputs = len(custom_bounds) // 2
            if zerobound == 1:
                numinputs == numinputs + 1
            neg_col = float(request.form['neg_col'])
            pos_col = float(request.form['pos_col'])
            num_col = 2*numinputs - 1
            numbshades_pos = len(positive_bounds) - 1
            numbshades_neg = len(negative_bounds) - 1
            numbshades = numinputs - 1
            if neg_col == 1: #red
                if max(abs(x) for x in custom_bounds) >= 2:
                    dcoln = dredc
                else:
                    dcoln = redc
            elif neg_col == 2: #blue
                if max(abs(x) for x in custom_bounds) >= 2:
                    dcoln = dbluec
                else: 
                    dcoln = bluec
            elif neg_col == 3: #green
                dcoln = greenc
            elif neg_col == 4: #yellow
                dcoln = yellowc
            elif neg_col == 5: #purple
                dcoln = purplec
            elif neg_col == 6: #black
                dcoln = blackc
            if pos_col == 1: #blue
                if max(abs(x) for x in custom_bounds) >= 2:
                    dcolp = dbluec
                else:
                    dcolp = bluec
            elif pos_col == 2: #red
                if max(abs(x) for x in custom_bounds) >= 2:
                    dcolp = dredc
                else:
                    dcolp = redc
            elif pos_col == 3: #green
                dcolp = greenc
            elif pos_col == 4: #yellow
                dcolp = yellowc
            elif pos_col == 5: #purple
                dcolp = purplec
            elif pos_col == 6: #black
                dcolp = blackc
            c1 = dcoln
            c2 = whitec
            c3 = dcolp
            col_mid = [mc.to_hex(c2)]
            cols1 = [colorFader(c1,c2,x/numbshades_neg) for x in range(numbshades_neg+1)]
            cols2 = [colorFader(c2,c3,x/numbshades_pos) for x in range(numbshades_pos+1)]
            if zerobound == 1:
                custom_colors = cols1[:-1] +cols2[1:]
            if zerobound == 0:
                custom_colors = cols1[:-1] + col_mid + cols2[1:]
        elif custb == 4: #INPUT COL AND BOUNDS
            max_range = None
            custom_bounds = []
            negative_bounds = []
            positive_bounds = []
            for i in range(1, 9):
                nb = request.form.get(f'inputbn{i}4', '').strip()
                pb = request.form.get(f'inputbp{i}4', '').strip() 
                if nb:
                    try:
                        nb = float(nb)
                        negative_bounds.append(nb)
                    except ValueError:
                        pass
                if pb:
                    try:
                        pb = float(pb)
                        positive_bounds.append(pb) 
                    except ValueError:
                        pass 
            if zerobound == 1:
                custom_bounds.append(0)
            custom_bounds.extend(negative_bounds)
            custom_bounds.extend(positive_bounds)
            custom_bounds = sorted(custom_bounds)
            numinputs = len(custom_bounds) // 2
            numbshades_pos = len(positive_bounds) - 1
            numbshades_neg = len(negative_bounds) - 1
            numbshades = numinputs - 1
            colorchoicetype = float(request.form['optionc'])
            if colorchoicetype == 1:
                negshadein = str(request.form['negcolorid']).strip()
                posshadein = str(request.form['poscolorid']).strip()
                c1 = ensure_hex_color(negshadein)
                c2 = whitec
                c3 = ensure_hex_color(posshadein)
                col_mid = [mc.to_hex(c2)]
                if zerobound == 1:
                    cols1 = [colorFader(c1,c2,x/numbshades_neg) for x in range(numbshades_neg)]
                    cols2 = [colorFader(c2,c3,x/numbshades_pos) for x in range(numbshades_pos)]
                    custom_colors = cols1[:-1] +cols2[1:]
                if zerobound == 0:
                    cols1 = [colorFader(c1,c2,x/numbshades_neg) for x in range(numbshades_neg+1)]
                    cols2 = [colorFader(c2,c3,x/numbshades_pos) for x in range(numbshades_pos+1)]
                    custom_colors = cols1[:-1] + col_mid + cols2[1:]
            elif colorchoicetype ==2: 
                custom_colors = []
                nc_list = []
                pc_list = []
                col_mid = mc.to_rgb(whitec) 
                for i in range(1, 9):
                    nc = request.form.get(f'ninputcol{i}', '').strip()
                    pc = request.form.get(f'pinputcol{i}', '').strip()
                    nc = ensure_hex_color(nc)
                    pc = ensure_hex_color(pc)
                    if nc:
                        try:
                            nc_rgb = mc.to_rgb(nc)
                            nc_list.append(nc_rgb)
                        except ValueError:
                            pass
                    if pc:
                        try:
                            pc_rgb = mc.to_rgb(pc)
                            pc_list.append(pc_rgb)
                        except ValueError:
                            pass
                for nc in nc_list:
                    custom_colors.append(nc)
                if zerobound == 0:
                    custom_colors.append(col_mid)
                for pc in pc_list:
                    custom_colors.append(pc)
        dtime = int(request.form.get('droptime', 0))
        dpept = int(request.form.get('droppept', 0))
        dprot = int(request.form.get('dropprot', 0))
        renumdict = int(request.form.get('renumdict', 0))
        #iddict = int(request.form.get('mutiddict', 0))
        #mutdict = int(request.form.get('mutdict', 0))
        muthandling = int(request.form.get('optionmut',0))
        slist = int(request.form.get('optionst', 1))
        annotateHM = int(request.form.get('annotateHM',0))
        chaindictuse = int(request.form.get('chaindictuse',0))
        if annotateHM == 1:
            addvalHM = True
        else: 
            addvalHM = False
        if renumdict == 1:
            renumbering_dict = {}
            numrenum = int(request.form.get('numrenum',0))
            for i in range(1, numrenum + 1):
                key = str(request.form[f'key{i}']).strip()
                value = request.form[f'value{i}']
                if key:
                    renumbering_dict[key] = value
        if chaindictuse == 1:
            print('went into chaindict')
            chain_dict = {}
            numchaindict = int(request.form.get('numchaindict',0))
            for i in range(1, numchaindict + 1):
                chainkey = str(request.form[f'keychain{i}']).strip() 
                chainvalue = str(request.form[f'valuechain{i}']) 
                if chainkey:
                    chain_dict[chainkey] = chainvalue
        if muthandling == 2:
            mut_id_dict = {}
            numMutProt = int(request.form.get('numMutProt',0))
            for i in range(1, numMutProt + 1):
                keyMutProt = str(request.form[f'keyMutProt{i}']).strip() 
                valueMutProt = str(request.form[f'valueMutProt{i}']) 
                if keyMutProt:
                    mut_id_dict[keyMutProt] = valueMutProt
            mutation_dict = {}
            numMutRes = int(request.form.get('numMutProt',0))
            for i in range(1, numMutRes + 1):
                resPos_raw = str(request.form[f'resPos{i}']).strip()
                if resPos_raw:
                    resPos = int(resPos_raw)
                    WTRes = str(request.form[f'WTRes{i}'])
                    mutation_dict[resPos] = WTRes
            correlation_data = [] 
            for i in range(1, numMutProt + 1):
                keyMutProt = str(request.form.get(f'keyMutProt{i}', '')).strip()
                valueMutProt = str(request.form.get(f'valueMutProt{i}', '')).strip()
                resPos_raw = str(request.form.get(f'resPos{i}', '')).strip()
                WTRes = str(request.form.get(f'WTRes{i}', '')).strip()
                if keyMutProt and resPos_raw: 
                    correlation_data.append({
                        'keyMutProt': keyMutProt,
                        'valueMutProt': valueMutProt,
                        'resPos': int(resPos_raw),
                        'WTRes': WTRes
                    })
        if dtime == 1:
            num_timedao = int(request.form.get('numtimedao', 7))
            drop_times = []
            for i in range(0, num_timedao+7):
                value = request.form.get(f'dt{i}')
                if value:
                    try:
                        drop_times.append(float(value))
                    except ValueError:
                        pass
        #if PulseLabellingNow == 2:                                 # NOT REMAKING BC POTENTIALLY GRABBED FROM DATA
        #    state_of_interest = str(request.form.get('state_of_interest', 'NotProvided'))
        #    ref_time = float(request.form.get('ref_time',2020))
        #    numNewTimes = int(request.form.get('numTimePulse', 0))
        #    new_times = []
        #    for i in range(0, numNewTimes+3):
        #        valueNT = request.form.get(f'dNT{i}')
        #        if valueNT:
        #            try:
        #                new_times.append(float(valueNT))
        #            except ValueError:
        #                pass
        if dpept == 1:
            drop_pept = []
            num_dpept = int(request.form.get('numpeptd', 4))
            for i in range(1, num_dpept+4):
                dpeppro = request.form.get(f'dpeppro{i}', '')
                dpepst = float(request.form.get(f'dpepst{i}', 0))
                dpepend = float(request.form.get(f'dpepend{i}', 0))
                if dpeppro:
                    drop_pept.append((dpeppro, dpepst, dpepend))
        if dprot == 1:
            num_protdao = int(request.form.get('numprotdao', 7))
            drop_prot = []
            for i in range(0, num_protdao+7):
                value = request.form.get(f'dpro{i}')
                if value:
                    try:
                        drop_prot.append(float(value))
                    except ValueError:
                        pass
        if slist == 2:
            numStatelist = int(request.form.get('numStatelist', 7))
            SpecifyStateProt = int(request.form.get('StateProt', 0))
            state1_list = []
            state2_list = []
            for i in range(0, numStatelist+3):
                s1 = str(request.form.get(f's1{i}', ''))
                s2 = str(request.form.get(f's2{i}', ''))
                if SpecifyStateProt == 1:
                    prSt1 = str(request.form.get(f'prSt1{i}', ''))
                    prSt2 = str(request.form.get(f'prSt2{i}', ''))
                    s1 += prSt1
                    s2 += prSt2
                if s1 != '':
                    state1_list.append(s1)
                    if s2 == '':
                        s2 = str(request.form.get('s21', ''))
                if s2 != '':
                    state2_list.append(s2)
            state1_list = state1_list[:len(state2_list)]
        separate_plots_pls = float(request.form.get('separate_plots','0'))
        if separate_plots_pls == 1:
            plot_separate = 1
            plot_stacked = 0
            print('separate')
        else:    
            plot_stacked  = 1
            plot_separate = 0
            print('not separate')
        editing_figtitle = float(request.form.get('usealtxaxis', 0))
        editing_xaxis = float(request.form.get('usealtxaxis', 0))
        if editing_xaxis == 1:
            new_x_axis_title = str(request.form.get('altxaxistit','Invalid Text'))
        else:
            if h_or_v == 1 or h_or_v == 2:
                new_x_axis_title = 'Peptides'
            if h_or_v == 3:
                new_x_axis_title = 'Position'
            if h_or_v == 4:
                new_x_axis_title = 'Difference in Deuterium Uptake (Da)'
                if relativeUptakeCalc == 1:
                        new_x_axis_title = 'Relative Difference in Deuterium Uptake'
            if h_or_v == 6:
                new_x_axis_title = 'Position'
        editing_yaxis = float(request.form.get('usealtyaxis', 0))
        if editing_yaxis == 1:
            new_y_axis_title = str(request.form.get('altyaxistit','Invalid Text'))
        else:
            if h_or_v == 1 or h_or_v == 2:
                new_y_axis_title = 'H/D exchange time'
            if h_or_v == 3:
                new_y_axis_title = 'Change in Deuterium Uptake (Da)'
                if relativeUptakeCalc == 1:
                        new_y_axis_title = 'Relative Difference in Deuterium Uptake'
            if h_or_v == 4:
                new_y_axis_title = '-log$_{10}$(p-value)'
            if h_or_v == 6:
                new_y_axis_title = 'Peptide #'
        output_bitmap_dpi = 100
        dif_dpi = float(request.form.get('dif_dpi', 0))
        if dif_dpi != 0:
            output_bitmap_dpi = float(request.form.get('dpi_in'))
        if PDFgeneration == 1:
            output_bitmap_format = 'pdf'
        else:
            output_bitmap_format = 'png'
        download_pymol = int(request.form.get('download_pymol', 0))
        download_chimera = int(request.form.get('download_chimera', 0))
        if editing_figtitle == 1:
            title = str(request.form.get('alttitleword','Invalid Text Entry'))
            newtitle = str(request.form.get('alttitleword','Invalid Text'))
        altwoodsthick = float(request.form.get('altwoodsthick',0))
        if altwoodsthick == 1:
            woodscolthick = float(request.form.get('woodscolthick',3))
            woodsbackthick = float(request.form.get('woodsbackthick',0.5))
        woodsplotYbyglobal = float(request.form.get('woodsplotYbyglobal',1))
        if usealtpad == 1:
            padthick = float(request.form.get('altpad',20))
            spadthick = float(request.form.get('altpads',20))
        domainlabel = float(request.form.get('domainlabel',0))
        if domainlabel == 1:
            numdomain = int(request.form.get('numdomain', 4))
            padthick = padthick + 40 + (font_size_title-32)*2
            for i in range(1, numdomain+4):
                domainName = request.form.get(f'domainName{i}', '')
                dompepst = float(request.form.get(f'dompepst{i}', 0)) - 1
                dompepend = float(request.form.get(f'dompepend{i}', 0))
                domColour = request.form.get(f'domColour{i}', '#000000')
                if domainName:
                    color_sections.append((dompepst, dompepend, domColour, domainName))
        usehmthick = float(request.form.get('usehmthick', 0))
        usestatesep = float(request.form.get('stackedplotdivl',0))
        usebordchange = float(request.form.get('hmplotbord',0))
        usehmtickchange = float(request.form.get('changehmtick',0))
        if usehmthick == 1:
            hmspacerthick = float(request.form.get('hmthickness',4)) 
            hmcolordivide = float(request.form.get('hmcolor',1))
            if hmcolordivide == 1: #White
                hmlcolor = whitec 
            elif hmcolordivide == 2: #Black
                hmlcolor = blackc
            elif hmcolordivide == 3: #Blue
                hmlcolor = dbluec
            elif hmcolordivide == 4: #Red
                hmlcolor = dredc
            elif hmcolordivide == 5: #Grey
                hmlcolor = greyc
        if usestatesep == 1:
            hmsepthick = float(request.form.get('hmsepthickness',4)) 
            hmsepcolor = float(request.form.get('hmsepcolor',2))
            if hmsepcolor == 1: #White
                hmsepcolorc = whitec 
            elif hmsepcolor == 2: #Black
                hmsepcolorc = blackc
            elif hmsepcolor == 3: #Blue
                hmsepcolorc = dbluec
            elif hmsepcolor == 4: #Red
                hmsepcolorc = dredc
            elif hmsepcolor == 5: #Grey
                hmsepcolorc = greyc    
        if usebordchange == 1:
            hmbordthick = float(request.form.get('hmbordthickness',10))
            hmbordcolor = float(request.form.get('hmbordcolor',2))
            if hmbordcolor == 1: #White
                hmbordcolorc = whitec 
            elif hmbordcolor == 2: #Black
                hmbordcolorc = blackc
            elif hmbordcolor == 3: #Blue
                hmbordcolorc = dbluec
            elif hmbordcolor == 4: #Red
                hmbordcolorc = dredc
            elif hmbordcolor == 5: #Grey
                hmbordcolorc = greyc 
        if usehmtickchange == 1:
            hmtw = float(request.form.get('hmtickwidth',4))
            hmtl = float(request.form.get('hmticklength',25))
            hmtickcolor = float(request.form.get('hmtickcolor',2))
            hmtickcolorlabel = float(request.form.get('hmtickcolorlabel',2))
            #Color Ticks
            if hmtickcolor == 1: #White
                hmtcolor = whitec
            elif hmtickcolor == 2: #Black
                hmtcolor = blackc
            elif hmtickcolor == 3: #Blue
                hmtcolor = dbluec
            elif hmtickcolor == 4: #Red
                hmtcolor = dredc
            elif hmtickcolor == 5: #Grey
                hmtcolor = greyc
            #Color Labels
            if hmtickcolorlabel == 1: #White
                hmtcolor_labels = whitec 
            elif hmtickcolorlabel == 2: #Black
                hmtcolor_labels = blackc
            elif hmtickcolorlabel == 3: #Blue
                hmtcolor_labels = dbluec
            elif hmtickcolorlabel == 4: #Red
                hmtcolor_labels = dredc
            elif hmtickcolorlabel == 5: #Grey
                hmtcolor_labels = greyc
        padthick = 20 #padding for titles
        spadthick = 10 #padding for axis labels
        usealtpad = float(request.form.get('usealtpad',20))
        if usealtpad == 1:
            padthick = float(request.form.get('altpad',20))
            spadthick = float(request.form.get('altpads',20))
        what_t_unit = 0
        all_in_min = float(request.form.get('all_in_min',0))
        if all_in_min == 1:
            what_t_unit = float(request.form.get('what_t_unit',0))
        p_threshold = request.form.get('pthresh', 0.05)
        altvolc = request.form.get('altvolc', 0)
        scattercolor = blackc
        if altvolc == 1:
            altvolccol = request.form.get('altvolccol', blackc)
            if is_valid_hex_color(altvolccol) == True:
                scattercolor = altvolccol
            else:
                scattercolor = blackc
        else:
            scattercolor = blackc   
        #Setup colormap color list and bin boundaries:
        if custom_colors != None: #custom colours (input colours and bounds option) 
            if custom_colors:
                colors = custom_colors
                if custom_bounds != None:
                    if custom_bounds:
                        bounds = custom_bounds
                        print('bounds and ccol')
                else:
                    bounds = make_bounds(len(colors))
                    print('no bounds found')
        elif custom_bounds != None: #Just custom bounds (choose col set bound option)
            if custom_bounds:
                print("The user has provided custom_bounds list and no custom_colors list. Using custom_bounds to calculate colors, and ingnoring max_range and num_shades variables.\n")
                b = np.array(custom_bounds)
                # the 3 variables below assume there is a zero value in the custom_bounds list, which divides all data into positive colors and negative colors:
                neg_shades = len(b[b<0])
                pos_shades = len(b[b>0])
                col_mid = []
                # if zero is absent from the custom_bounds list, then set the middle color to c2 (default white), but only if custom_bounds list has at least one value on each side of zero!
                # Also, decrease number of shades on each side of zero by one:
                if len(b[b==0])==0:
                    if neg_shades>0 and pos_shades>0:
                        col_mid = [ mc.to_hex(c2) ]
                        if zerobound == 0:
                            neg_shades = len(b[b<0])-1
                            pos_shades = len(b[b>0])-1
                        if zerobound == 1:
                            neg_shades = len(b[b<0])
                            pos_shades = len(b[b>0])
                if neg_shades > 0:
                    cols1 = [colorFader(c1,c2,x/neg_shades) for x in range(neg_shades+1)]
                else:
                    cols1 = [ mc.to_hex(c2) ]
                if pos_shades>0:
                    cols2 = [colorFader(c2,c3,x/pos_shades) for x in range(pos_shades+1)]
                else:
                    cols2 = [ mc.to_hex(c2) ]
                if zerobound == 0:
                    colors = cols1[:-1] + col_mid + cols2[1:]
                if zerobound == 1:
                    colors = cols1[:-1] + cols2[1:]
                bounds = custom_bounds
                #print(bounds)
                #print('no ccol')
        else: #no custom colours or bounds (auto or set range option)
            #c_missing = '#bdbdbd' # gray;
            num_shades = round(num_shades) # Convert num_shades to integer if it's a float
            print('num_shades:', num_shades)
            print('custb:', custb)
            print('zerobound:', zerobound)
            if zerobound == 0 and custb != 3:
                cols1 = [colorFader(c1,c2,x/num_shades) for x in range(num_shades+1)]
                cols2 = [colorFader(c2,c3,x/num_shades) for x in range(num_shades+1)]
                if funkybound == 1: #only positive
                    colors = cols2[1:] 
                elif funkybound == 2: #only negative
                    colors = cols1
                    print('not custb 3 and fb == 2, colors:', colors)
                else: #positive and negative
                    colors = cols1 + cols2[1:]
            elif zerobound == 0 and custb == 3:
                if funkybound == 1 or funkybound == 2:
                    cols1 = [colorFader(c1,c2,x/num_shades) for x in range(num_shades-1)]
                    cols2 = [colorFader(c2,c3,x/num_shades) for x in range(num_shades)]
                    if funkybound == 1:
                        colors = [c2] + cols2[1:]
                        print('custb 3 and fb 1, colors:', colors)
                    if funkybound == 2:
                        colors = cols1 + [c2]
                        print('custb 3 and fb 2, colors:', colors)
                else:
                    cols1 = [colorFader(c1,c2,x/num_shades) for x in range(num_shades+1)]
                    cols2 = [colorFader(c2,c3,x/num_shades) for x in range(num_shades+1)]
                    colors = cols1 + cols2[1:]
            if zerobound == 1:
                cols1 = [colorFader(c1, c2, (x/num_shades)) for x in range(num_shades+1)]
                cols2 = [colorFader(c2, c3, (x/num_shades)) for x in range(num_shades+1)]
                if funkybound == 1:
                    colors = cols2[1:]
                elif funkybound == 2:
                    colors = cols1[:-1]
                    print('in 0 bound')
                else:
                    colors = cols1[:-1] + cols2[1:]
            bounds = make_bounds(len(colors))
            if zerobound == 1:
                bounds = np.append(bounds, 0)
                bounds = list(set(bounds))
                bounds = np.sort(bounds)
            custom_bounds = bounds.flatten().tolist()
            custom_colors = colors
        ###   Create a Dict of Necessary Variables   ###
        variables = {
            key: value
            for key, value in {
                "fontchoice": locals().get("fontchoice", None),
                "PDFgeneration": locals().get("PDFgeneration", None),
                "alt_c2": locals().get("alt_c2", None),
                "c2": locals().get("c2", None),
                "alt_cmissing": locals().get("alt_cmissing", None),
                "c_missing": locals().get("c_missing", None),
                "color_by_heatmap": locals().get("color_by_heatmap", None),
                "colcutopt": locals().get("colcutopt", None),
                "colcutoff": locals().get("colcutoff", None),
                "woodsdimen": locals().get("woodsdimen", None),
                "woodsx": locals().get("woodsx", None),
                "woodsy": locals().get("woodsy", None),
                "woodscol": locals().get("woodscol", None),
                "woodscolpos": locals().get("woodscolpos", None),
                "woodscolneu": locals().get("woodscolneu", None),
                "woodscolneg": locals().get("woodscolneg", None),
                "nolines": locals().get("nolines", None),
                "font_size": locals().get("font_size", None),
                "font_size_title": locals().get("font_size_title", None),
                "font_size_ticklabel": locals().get("font_size_ticklabel", None),
                "fontcolorc": locals().get("fontcolorc", None),
                "h_or_v": locals().get("h_or_v", None),
                "zerobound": locals().get("zerobound", None),
                "funkybound": locals().get("funkybound", None),
                "globalmax_delta": locals().get("globalmax_delta", None),
                "p_threshold": locals().get("p_threshold", None),
                "custb": locals().get("custb", None),
                "max_range": locals().get("max_range", None),
                "num_shades": locals().get("num_shades", None),
                "alt_col": locals().get("alt_col", None),
                "c1": locals().get("c1", None),
                "c3": locals().get("c3", None),
                "neg_col": locals().get("neg_col", None),
                "pos_col": locals().get("pos_col", None),
                "negshadein": locals().get("negshadein", None),
                "posshadein": locals().get("posshadein", None),
                "dtime": locals().get("dtime", None),
                "dpept": locals().get("dpept", None),
                "dprot": locals().get("dprot", None),
                "renumdict": locals().get("renumdict", None),
                "muthandling": locals().get("muthandling", None),
                "slist": locals().get("slist", None),
                "annotateHM": locals().get("annotateHM", None),
                "chaindictuse": locals().get("chaindictuse", None),
                "numrenum": locals().get("numrenum", None),
                "numchaindict": locals().get("numchaindict", None),
                "numMutProt": locals().get("numMutProt", None),
                "keyMutProt": locals().get("keyMutProt", None),
                "valueMutProt": locals().get("valueMutProt", None),
                "numMutRes": locals().get("numMutRes", None),
                "resPos_raw": locals().get("resPos_raw", None),
                "WTRes": locals().get("WTRes", None),
                "num_timedao": locals().get("num_timedao", None),
                "value": locals().get("value", None),
                "state_of_interest": locals().get("state_of_interest", None),
                "ref_time": locals().get("ref_time", None),
                "numNewTimes": locals().get("numNewTimes", None),
                "valueNT": locals().get("valueNT", None),
                "num_dpept": locals().get("num_dpept", None),
                "dpeppro": locals().get("dpeppro", None),
                "dpepst": locals().get("dpepst", None),
                "dpepend": locals().get("dpepend", None),
                "dpro1": locals().get("dpro1", None),
                "dpro2": locals().get("dpro2", None),
                "dpro3": locals().get("dpro3", None),
                "dpro4": locals().get("dpro4", None),
                "dpro5": locals().get("dpro5", None),
                "dpro6": locals().get("dpro6", None),
                "s1": locals().get("s1", None),
                "s2": locals().get("s2", None),
                "separate_plots_pls": locals().get("separate_plots_pls", None),
                "editing_fig_title": locals().get("editing_fig_title", None),
                "editing_xaxis": locals().get("editing_xaxis", None),
                "new_x_axis_title": locals().get("new_x_axis_title", None),
                "editing_yaxis": locals().get("editing_yaxis", None),
                "new_y_axis_title": locals().get("new_y_axis_title", None),
                "dif_dpi": locals().get("dif_dpi", None),
                "output_bitmap_dpi": locals().get("output_bitmap_dpi", None),
                "altvolc": locals().get("altvolc", None),
                "altvolccol": locals().get("altvolccol", None),
                "title": locals().get("title", None),
                "all_in_min": locals().get("all_in_min", None),
                "what_t_unit": locals().get("what_t_unit", None),
                "usehmthick": locals().get("usehmthick", None),
                "usestatesep": locals().get("usestatesep", None),
                "usebordchange": locals().get("usebordchange", None),
                "usehmtickchange": locals().get("usehmtickchange", None),
                "hmspacerthick": locals().get("hmspacerthick", None),
                "hmcolordivide": locals().get("hmcolordivide", None),
                "hmsepthick": locals().get("hmsepthick", None),
                "hmsepcolor": locals().get("hmsepcolor", None),
                "hmbordthick": locals().get("hmbordthick", None),
                "hmbordcolor": locals().get("hmbordcolor", None),
                "hmtw": locals().get("hmtw", None),
                "hmtl": locals().get("hmtl", None),
                "hmtickcolor": locals().get("hmtickcolor", None),
                "hmtickcolorlabel": locals().get("hmtickcolorlabel", None),
                "usealtpad": locals().get("usealtpad", None),
                "padthick": locals().get("padthick", None),
                "spadthick": locals().get("spadthick", None),
                "domainlabel": locals().get("domainlabel", None),
                "numdomain": locals().get("numdomain", None),
                "domainName": locals().get("domainName", None),
                "dompepst": locals().get("dompepst", None),
                "dompepend": locals().get("dompepend", None),
                "domColour": locals().get("domColour", None),
                "newtitle": locals().get("newtitle", None),
                "altwoodsthick": locals().get("altwoodsthick", None),
                "woodscolthick": locals().get("woodscolthick", None),
                "woodsbackthick": locals().get("woodsbackthick", None),
                "woodsplotYbyglobal": locals().get("woodsplotYbyglobal", None),
                "download_pymol": locals().get("download_pymol", None),
                "download_chimera": locals().get("download_chimera", None),
                "fontcolor": locals().get("fontcolor", None),
                "custom_colors": locals().get("custom_colors", None),
                "custom_bounds": locals().get("custom_bounds", None),
                "negative_bounds": locals().get("negative_bounds", None),
                "positive_bounds": locals().get("positive_bounds", None),
                "renumbering_dict": locals().get("renumbering_dict", None),
                "chain_dict": locals().get("chain_dict", None),
                "mut_id_dict": locals().get("mut_id_dict", None),
                "mutation_dict": locals().get("mutation_dict", None),
                "correlation_data": locals().get("correlation_data", None),
                "drop_times": locals().get("drop_times", None),
                "new_times": locals().get("new_times", None),
                "drop_pept": locals().get("drop_pept", None),
                "drop_prot": locals().get("drop_prot", None),
                "state1_list": locals().get("state1_list", None),
                "state2_list": locals().get("state2_list", None),
                "plot_separate": locals().get("plot_separate", None),
                "plot_stacked": locals().get("plot_stacked", None),
                "output_bitmap_format": locals().get("output_bitmap_format", None),
                "color_sections": locals().get("color_sections", None),
                "hmlcolor": locals().get("hmlcolor", None),
                "hmsepcolorc": locals().get("hmsepcolorc", None),
                "hmbordcolorc": locals().get("hmbordcolorc", None),
                "hmtcolor": locals().get("hmtcolor", None),
                "hmtcolor_labels": locals().get("hmtcolor_labels", None),
                "scattercolor": locals().get("scattercolor", None),
                "blockposwoods": locals().get("blockposwoods", None),
                "wpadthick": locals().get("wpadthick", None)
            }.items()
            if value is not None
        }
        ###   Create and Upload File   ###
        output_file = BLOB_NAMEsv
        with open(output_file, mode="w", newline="") as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(["Variable Name", "Data Type", "Value"])
            for var_name, value in variables.items():
                if value is not None:
                    if isinstance(value, str):
                        writer.writerow([var_name, type(value).__name__, value])
                    else:
                        writer.writerow([var_name, type(value).__name__, repr(value)])
        print(f"CSV file '{output_file}' has been generated.")
        with open(output_file, "rb") as data:
            blob_client = container_client.get_blob_client(BLOB_NAMEsv)
            blob_client.upload_blob(data, overwrite=True)
            print(f"File '{output_file}' uploaded to blob '{BLOB_NAMEsv}' in container '{CONTAINER_NAMEsv}'.")
        ###   Download File to User   ###
        blob_clientsv = blob_service_clientsv.get_blob_client(container=CONTAINER_NAMEsv, blob=BLOB_NAMEsv)
        blob_datasv = blob_clientsv.download_blob()
        blob_contentsv = io.BytesIO(blob_datasv.readall())
        return send_file(
            blob_contentsv,
            as_attachment=True,
            download_name=BLOB_NAMEsv)
    else:
        if scatter_plot == 1:
            scatter_zip_name = generate_unique_filename('scatter_plots', 'zip', timestamp)
            scatter_blob_name = upload_to_blob_storage(container_client, scatter_buffer.getvalue(), scatter_zip_name)
            scatter_blob_stream = download_blob_as_bytes(container_client, scatter_blob_name)
            response = make_response(send_file(scatter_blob_stream, as_attachment=True, download_name=scatter_zip_name))
            response.headers['Content-Type'] = 'application/zip'
            response.headers['Content-Disposition'] = f'attachment; filename={scatter_zip_name}'
            print('scatterplot', scatter_zip_name)
            #flash("Volcano Plots Sent!")
            return response
        else:
            if download_pymol == 1:
                print('pymol')
                response = download_collect(chain_dict, chain_id, colors, bounds, all_deltas, all_states, buffer, heatmap_buffer, output_buffer, mk_pymol)
                #flash("PyMOL and Plot Files Sent!")
                return response
            elif download_chimera == 1:
                print('chimerax')
                #chain_tocol = str(request.form.get('chain_tocol', '#1'))
                response = download_collect_chimerax(chain_dict, chain_id, colors, bounds, all_deltas, all_states, buffer, heatmap_buffer, output_buffer, mk_chimerax)
                #flash("ChimeraX and Plot Files Sent!")
                return response
            elif plot_w != 0:
                zip_filename = generate_unique_filename('WoodsPlots', 'zip', timestamp)
                upload_to_blob_storage(container_client, output_buffer, zip_filename)
                zip_blob_stream = download_blob_as_bytes(container_client, zip_filename)
                response = make_response(send_file(zip_blob_stream, as_attachment=True, download_name=zip_filename))
                response.headers['Content-Type'] = 'application/zip'
                response.headers['Content-Disposition'] = f'attachment; filename={zip_filename}'
                output_buffer.close()
                print('woods', zip_filename)
                #flash("Woods Plots Sent!")
                return response
            elif plot_l != 0:
                zip_filename = generate_unique_filename('LadderPlots', 'zip', timestamp)
                upload_to_blob_storage(container_client, output_buffer, zip_filename)
                zip_blob_stream = download_blob_as_bytes(container_client, zip_filename)
                response = make_response(send_file(zip_blob_stream, as_attachment=True, download_name=zip_filename))
                response.headers['Content-Type'] = 'application/zip'
                response.headers['Content-Disposition'] = f'attachment; filename={zip_filename}'
                output_buffer.close()
                print('ladders', zip_filename)
                #flash("Woods Plots Sent!")
                return response
            elif len(all_prot) > 1:
                heatmaps_zip_name = generate_unique_filename('heatmaps', 'zip', timestamp)
                upload_to_blob_storage(container_client, heatmap_buffer.getvalue(), heatmaps_zip_name)
                zip_blob_stream = download_blob_as_bytes(container_client, heatmaps_zip_name)
                response = make_response(send_file(zip_blob_stream, as_attachment=True, download_name=heatmaps_zip_name))
                response.headers['Content-Type'] = 'application/zip'
                response.headers['Content-Disposition'] = f'attachment; filename={heatmaps_zip_name}'
                print('1 heatmap', heatmaps_zip_name)
                #flash("Heatmap Files Sent!")
                return response
            elif PDFgeneration == 1:
                heatmaps_zip_name = generate_unique_filename('heatmap', 'zip', timestamp)
                upload_to_blob_storage(container_client, heatmap_buffer.getvalue(), heatmaps_zip_name)
                zip_blob_stream = download_blob_as_bytes(container_client, heatmaps_zip_name)
                response = make_response(send_file(zip_blob_stream, as_attachment=True, download_name=heatmaps_zip_name))
                response.headers['Content-Type'] = 'application/zip'
                response.headers['Content-Disposition'] = f'attachment; filename={heatmaps_zip_name}'
                print('pdf', heatmaps_zip_name)
                #flash("Heatmap Files Sent!")
                return response
            else:
                heatmaps_zip_name = generate_unique_filename('heatmap', 'zip', timestamp)
                upload_to_blob_storage(container_client, heatmap_buffer.getvalue(), heatmaps_zip_name)
                zip_blob_stream = download_blob_as_bytes(container_client, heatmaps_zip_name)
                response = make_response(send_file(zip_blob_stream, as_attachment=True, download_name=heatmaps_zip_name))
                response.headers['Content-Type'] = 'application/zip'
                response.headers['Content-Disposition'] = f'attachment; filename={heatmaps_zip_name}'
                print('all else', heatmaps_zip_name)
                #flash("Heatmap Files Sent!")
                return response
    #except:
    #    error = "ERROR HERE"
    #    print('IN THE EXCEPTION')
    #    return render_template('index.html', output=None, error=error)



    #################################################################################################################################################################
    #################################################################################################################################################################
    #############################################################            OTHER FUNCTIONS            #############################################################
    #################################################################################################################################################################
    #################################################################################################################################################################



@app.route('/reset_values', methods=['POST'])
def reset_values():
    return render_template('index.html', output=None)


@app.route('/help_doc', methods=['Post'])
def help_doc():
    AZURE_STORAGE_CONNECTION_STRINGhd = "ThisIsObtainedOnAzure"
    CONTAINER_NAMEhd = "AContainerForTemporaryOutputStorage"
    BLOB_NAMEhd = 'HDgraphiX Help Doc.pdf'
    blob_service_clienthd = BlobServiceClient.from_connection_string(AZURE_STORAGE_CONNECTION_STRINGhd)
    blob_clienthd = blob_service_clienthd.get_blob_client(container=CONTAINER_NAMEhd, blob=BLOB_NAMEhd)
    blob_datahd = blob_clienthd.download_blob()
    blob_contenthd = io.BytesIO(blob_datahd.readall())
    # Send the file to the user
    return send_file(
        blob_contenthd,
        as_attachment=True,
        download_name=BLOB_NAMEhd)


@app.route('/sample_file', methods=['Post']) #DYNAMX
def sample_file():
    AZURE_STORAGE_CONNECTION_STRINGsf = "ThisIsObtainedOnAzure"
    CONTAINER_NAMEsf = "AContainerForTemporaryOutputStorage"
    BLOB_NAMEsf = 'HDgraphiX Sample DynamX.csv'
    blob_service_clientsf = BlobServiceClient.from_connection_string(AZURE_STORAGE_CONNECTION_STRINGsf)
    blob_clientsf = blob_service_clientsf.get_blob_client(container=CONTAINER_NAMEsf, blob=BLOB_NAMEsf)
    blob_datasf = blob_clientsf.download_blob()
    blob_contentsf = io.BytesIO(blob_datasf.readall())
    # Send the file to the user
    return send_file(
        blob_contentsf,
        as_attachment=True,
        download_name=BLOB_NAMEsf)


@app.route('/sample_file_pool', methods=['Post']) #HDEXAMINER POOL
def sample_file_pool():
    AZURE_STORAGE_CONNECTION_STRINGsf = "ThisIsObtainedOnAzure"
    CONTAINER_NAMEsf = "AContainerForTemporaryOutputStorage"
    BLOB_NAMEsf = 'HDgraphiX Sample HDExaminer Pool.csv'
    blob_service_clientsf = BlobServiceClient.from_connection_string(AZURE_STORAGE_CONNECTION_STRINGsf)
    blob_clientsf = blob_service_clientsf.get_blob_client(container=CONTAINER_NAMEsf, blob=BLOB_NAMEsf)
    blob_datasf = blob_clientsf.download_blob()
    blob_contentsf = io.BytesIO(blob_datasf.readall())
    # Send the file to the user
    return send_file(blob_contentsf, as_attachment=True, download_name=BLOB_NAMEsf)


@app.route('/sample_file_summary', methods=['Post']) #HDEXAMINER SUMMARY
def sample_file_summary():
    AZURE_STORAGE_CONNECTION_STRINGsf = "ThisIsObtainedOnAzure"
    CONTAINER_NAMEsf = "AContainerForTemporaryOutputStorage"
    BLOB_NAMEsf = 'HDgraphiX Sample HDExaminer Summary.csv'
    blob_service_clientsf = BlobServiceClient.from_connection_string(AZURE_STORAGE_CONNECTION_STRINGsf)
    blob_clientsf = blob_service_clientsf.get_blob_client(container=CONTAINER_NAMEsf, blob=BLOB_NAMEsf)
    blob_datasf = blob_clientsf.download_blob()
    blob_contentsf = io.BytesIO(blob_datasf.readall())
    # Send the file to the user
    return send_file(blob_contentsf, as_attachment=True, download_name=BLOB_NAMEsf)


@app.route('/sample_file_workbench', methods=['Post']) #HDEXAMINER SUMMARY
def sample_file_workbench():
    AZURE_STORAGE_CONNECTION_STRINGsf = "ThisIsObtainedOnAzure"
    CONTAINER_NAMEsf = "AContainerForTemporaryOutputStorage"
    BLOB_NAMEsf = 'HDgraphiX Sample HDXWorkbench.csv'
    blob_service_clientsf = BlobServiceClient.from_connection_string(AZURE_STORAGE_CONNECTION_STRINGsf)
    blob_clientsf = blob_service_clientsf.get_blob_client(container=CONTAINER_NAMEsf, blob=BLOB_NAMEsf)
    blob_datasf = blob_clientsf.download_blob()
    blob_contentsf = io.BytesIO(blob_datasf.readall())
    # Send the file to the user
    return send_file(blob_contentsf, as_attachment=True, download_name=BLOB_NAMEsf)


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=8000, debug=True)
    #app.run(debug=False)


# In[ ]:





# In[ ]:




