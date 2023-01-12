import matplotlib as mpl 
import numpy as np
import matplotlib.pyplot as plt
import sys

def read_data(data_file):
    with open(data_file, "r") as data:
        for i, line in enumerate(data):
            if i == 0:
                col_names = line.strip().split("\t")
                col_dict = {}
                for name in col_names:
                    col_dict[name] = []
            else:
                line = line.strip().split("\t")
                for i, value in enumerate(line):
                    key = col_names[i]
                    col_dict[key].append(value)
    return col_dict

def get_mito_stress_test(col_dict, group_dict, well_data_dict, output_filename):
    output_rows = ["Basal_Respiration", "Proton_leak", "Max_Respiration", "ATP_Linked", "Reserve_Capacity" ]
    with open(output_filename, 'w') as output_file:
        count = 0
        for group in group_dict.keys():
            #print(group)
            col_labels = []
            cols = set(group_dict[group])
            well_type = list(set(group_dict[group]))
            #print(group_dict[group])
            rows = group_dict[group].count(well_type[0])
            #print(rows,len(cols))
            OCR_matrix = np.zeros((rows, len(cols)))
            #print(OCR_matrix.shape)
            ECAR_matrix = np.zeros((rows, len(cols)))
            time_points = list(set(col_dict["Time"]))
            time_dict = {}
            well_dict = {}
            wells = group_dict[group]
            grouped_wells = list(set(wells))
            grouped_wells.sort()
            time_floats = [float(x) for x in time_points]
            time_floats.sort()
            for well_col, well_name in enumerate(grouped_wells):
                well_dict[well_name] = well_col
            for start, time_point in enumerate(time_floats):
                time_dict[time_point] = start
            #print("time points", time_points)
            row_labels = list(time_points)
            #print("wells", wells)
            #print("col, rows", len(wells), rows)
            for j, well in enumerate(well_type):
                print(well)
                col_labels.append(well)
                well_col = well_dict[well]
                for i in range(rows):
                    label = group + "_" + well + "_" + time_points[i]
                    time, OCR, ECAR = well_data_dict[label]
                    row = time_dict[time]
                    #print(i,j)
                    OCR_matrix[row][well_col] = float(OCR)
                    ECAR_matrix[row][well_col] = float(ECAR)
            print(OCR_matrix)
            non_o2_mito_OCR_matrix = np.subtract(OCR_matrix, OCR_matrix[-1][:])
            #print(non_o2_mito_OCR_matrix)
            non_o2_mito_ECAR_matrix = np.subtract(ECAR_matrix, ECAR_matrix[-1][:])
            ave_change_mito_OCR_matrix = np.zeros((rows, len(cols)))
            ave_change_mito_ECAR_matrix = np.zeros((rows, len(cols)))
            #print(OCR_matrix.shape)
            #print(non_o2_mito_ECAR_matrix.shape)
            for j, well in enumerate(well_type):
                for i in range(rows-1):
                    #print(i,j)
                    time_diff = float(time_points[i+1]) - float(time_points[i])
                    time_diff_average_OCR = time_diff*(non_o2_mito_OCR_matrix[i+1][j] + non_o2_mito_OCR_matrix[i][j])/2
                    time_diff_average_ECAR = time_diff*(non_o2_mito_ECAR_matrix[i+1][j] + non_o2_mito_ECAR_matrix[i][j])/2
                    ave_change_mito_OCR_matrix[i][j] = time_diff_average_OCR
                    ave_change_mito_ECAR_matrix[i][j] = time_diff_average_ECAR
            # next step: Sum over some rows to get basel/proton_leak/max/reserve from OCR
            #print(ave_change_mito_OCR_matrix.shape)
            basal = ave_change_mito_OCR_matrix[:3][:].sum(axis=0)
            proton_leak = ave_change_mito_OCR_matrix[3:6][:].sum(axis=0)
            max_resp = ave_change_mito_OCR_matrix[6:9][:].sum(axis=0)
            atp_link = ave_change_mito_OCR_matrix[9:][:].sum(axis=0)
            res_capacity = np.subtract(max_resp, basal)
            #coupling_efficiency = atp_link/basal*100
            # proton_leak = average
            if count == 0:
                output_file.write(group + "\t")
            elif count > 0:
                output_file.write("\n\n" + group + "\t")
            else:
                print("error: ")
            count += 1
            for well in grouped_wells:
                output_file.write(well + "\t")
            proton_str = ''
            basal_str = ''
            max_resp_str = ''
            atp_link_str = ''
            res_capacity_str = ''
            #print(atp_link.shape)
            #print("shape", basal.shape)
            for j in range(basal.shape[0]):
                basal_val = basal[:][j]
                basal_str += "\t" + str(basal_val)
                proton_val = proton_leak[:][j]
                proton_str += "\t" + str(proton_val)
                max_resp_val = max_resp[:][j]
                max_resp_str += "\t" + str(max_resp_val)
                atp_link_val = atp_link[:][j]
                atp_link_str += "\t" + str(atp_link_val)
                res_capacity_val = res_capacity[:][j]
                res_capacity_str += "\t" + str(res_capacity_val)
            file_output_list = [basal_str, proton_str, max_resp_str, atp_link_str, res_capacity_str]
            for i, row_lab in enumerate(output_rows):
                output_file.write("\n" + output_rows[i] + file_output_list[i])

                    

        
""" # get basal glycolysis rate, glycolytic capacity, non-glycolytic acidification, and glycolytic reserve from ECAR
           basal_glycolysis = ave_change_mito_ECAR_matrix[:3][:].sum(axis=1)
            glycolytic_capacity =  ave_change_mito_ECAR_matrix[3:6][:].sum(axis=1)
            non_glyc_acid = ave_change_mito_ECAR_matrix[6:9][:].sum(axis=1)
            glycolytic_reserve = np.subtract(glycolytic_capacity, basal_glycolysis)"""

def seperate_data(col_dict):
    group_dict = {}
    well_data_dict = {}
    for i, group in enumerate(col_dict["Group"]):
        time = col_dict["Time"][i]
        well = col_dict["Well"][i]
        OCR = col_dict["OCR"][i]
        ECAR = col_dict["ECAR"][i]
        label = group + "_" + well + "_" + time
        well_data_dict[label] = [float(time), float(OCR), float(ECAR)]
        if group not in group_dict.keys():
            group_dict[group] = [well]
        else:
            group_dict[group].append(well)
    return group_dict, well_data_dict    

########
##MAIN##
########

usage = sys.argv[0] + "< File_input > < output_filename >"

if len(sys.argv) != 3:
    print(usage)
    sys.exit()

data_file = sys.argv[1]
output_file = sys.argv[2]

col_dict = read_data(data_file)
group_dict, well_data_dict = seperate_data(col_dict)
get_mito_stress_test(col_dict, group_dict, well_data_dict, output_file)

# example to run
# need data in tab delimited file
# python automate_seahorse.py <title of input.txt> <title of output.txt>

