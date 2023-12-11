########## imports and setup
import numpy as np
import matplotlib.pyplot as plt
import os
import json
import datetime as dt

# import functions from other file
import Voltage_base_analysis as vba
import flow_on_off_cycle_analysis as fca

#plot Options
import matplotlib as mpl
font = {#'family' : 'normal','weight' : 'bold',
        'size'   : 16
        #,'serif':['Helvetica']
        }
mpl.rc('font', **font)

plot_dir = (os.path.dirname(os.path.abspath(__file__)) + os.sep 
            + "output/pressure/")
sc_dir = (os.path.dirname(os.path.abspath(__file__)) + os.sep + "../" 
         + "SC_downloads/")
pressure_data_dir = (os.path.dirname(os.path.abspath(__file__)) + os.sep 
                    + "pressure_data/")
os.makedirs(pressure_data_dir, exist_ok=True)
os.makedirs(plot_dir, exist_ok=True)
#######################

def load_json(file_path):
    with open(file_path) as f:
        raw_dict  = json.load(f)
    return raw_dict

def str_to_dt(date_string):
    date_str_format = "%Y-%m-%d %H:%M:%S"
    date = dt.datetime.strptime(date_string,
                                    date_str_format
                                )
    return date

def pstr_to_mbar(p_str):
    #"74006" always preceeds the relevant pressure data
    # which is stored in the next 6 characters
    temp_str = p_str.split("74006")[1][0:6]
    # characters 0-3 are the mantissa 
    mantissa = temp_str[0:4]
    # char 4 and 5 are the exponent + 20 for mbar
    exp = temp_str[4:6]

    num_mbar = 1e-3*float(mantissa)  * 10**(int(exp) - 20)
    return num_mbar

def pstr_pg9_to_mbar(p_str):
    #"@253ACK1.25E-1"
    #"@253ACK" always preceeds the relevant pressure data
    # which is stored in the next 7 characters
    temp_str = p_str.split("@253ACK")[1][0:7]
    num_mbar = float(temp_str)
    return num_mbar

def pstr_alicat_to_mbar(p_str):
    #"A +00532 +025.14 +25.599 +14.995 +15.000 +025.33     H2 VOV",
    #"A +" always preceeds the relevant pressure data
    # which is stored in the next 5 characters
    temp_str = p_str.split("A +")[1][0:5]
    num_mbar = float(temp_str)
    return num_mbar


def make_p_dict(raw_dict, utc_offset = 2, pg_key = "pg60_pressure_string"):
    dates = np.array([str_to_dt(ls[0])
                    for ls in raw_dict[pg_key]])
    # move to German timezone
    dates = np.array([(
            date.astimezone(dt.timezone(
                dt.timedelta(hours=utc_offset))) 
                + dt.timedelta(hours=utc_offset)
                )
            for date in dates])
    if pg_key in ["pg9_pressure_string", "pg8_pressure_string"]:
        pressures = np.array([pstr_pg9_to_mbar(ls[1]) 
                    for ls in raw_dict[pg_key]])
    elif pg_key in ["alicat1_readout_string"]:
        pressures = np.array([pstr_alicat_to_mbar(ls[1]) 
                    for ls in raw_dict[pg_key]])
    else:
        pressures = np.array([pstr_to_mbar(ls[1]) 
                    for ls in raw_dict[pg_key]])
    p_dict = {"dates": dates, "pressure":pressures}
    return p_dict

def p_dict_date_select(p_dict,start_date, end_date):
    start,end = vba.select_date_indices(p_dict,
                                        start_date,
                                        end_date)
    p_dict_ds = {key:p_dict[key][start:end] for key in p_dict.keys()}
    return p_dict_ds
    
def plot_pressure(p_dict, plotname):

    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    ax1.plot(p_dict["dates"],
            p_dict["pressure"],
            "-",# markersize=1,
            # color = color,
            # label = label
            )
    # ax1.plot(t_avg_series-t_series[0],v_avg_series*1000,
    #         label = f"moving average {mavg_len}s")

    plt.xticks(rotation = 45)

    ax1.set_xlabel(r"Time")
    ax1.set_ylabel(r"pressure [mbar]")

    plt.grid(True)
    #plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname + "_pressure"
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()

def plot_pressure_mult(p_dicts, plotname, labels=[]):


    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    for i,p_dict in enumerate(p_dicts):
        ax1.plot(p_dict["dates"],
                p_dict["pressure"],
                "-",# markersize=1,
                # color = color,
                label = labels[i]
                )
        # ax1.plot(t_avg_series-t_series[0],v_avg_series*1000,
        #         label = f"moving average {mavg_len}s")

    plt.xticks(rotation = 45)

    ax1.set_xlabel(r"Time")
    ax1.set_ylabel(r"pressure [mbar]")

    plt.grid(True)
    plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname + "_pressure"
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)


    ax1.set_yscale("log")   
    plt.legend(shadow=True, ncol = 2, loc="lower left", bbox_to_anchor =(0,1))
    plt.tight_layout()
    plt.savefig(plot_dir + plotname + "_pressure"
                + '_log.{}'.format(format_im),
                format=format_im, dpi=dpi)         
    ax1.cla()



if __name__ =="__main__":
    ########### 2023-01-06 closed ballast valve pump shutdwon
    file_path = pressure_data_dir + "2023-01-06_PGs.json"
    raw_dict = load_json(file_path)
    pg_keys = ["pg8_pressure_string","pg9_pressure_string",
        "pg60_pressure_string","pg80_pressure_string","pg90_pressure_string",
        "pg100_pressure_string", "alicat1_readout_string" ]

    p_dicts = list([])
    print(p_dicts)
    for pg_key in pg_keys:
        print(pg_key)
        p_dict = make_p_dict(raw_dict, utc_offset=1,pg_key= pg_key
        #     p_dict_ds = p_dict_date_select(p_dict,
        #                     start_date = dt.datetime(2022, 11, 1,0, 0,
        #                         tzinfo=dt.timezone(dt.timedelta(hours=1))), 
        #                     end_date = dt.datetime(2022, 11, 7,12, 0,
        #                         tzinfo=dt.timezone(dt.timedelta(hours=1))),
                                                )
        p_dicts.append(p_dict)


    plot_pressure_mult(p_dicts, plotname = "2023-01-06_PGs",
                         labels=[key.split("_string")[0] for key in pg_keys])

#     ########### 2022-11-07 troubleshooting
#     file_path = sc_dir + "2022-11-07_PGs.json"
#     raw_dict = load_json(file_path)
#     p_dict = make_p_dict(raw_dict, utc_offset=1)
#     p_dict_ds = p_dict_date_select(p_dict,
#                         start_date = dt.datetime(2022, 11, 1,0, 0,
#                                 tzinfo=dt.timezone(dt.timedelta(hours=1))), 
#                         end_date = dt.datetime(2022, 11, 7,12, 0,
#                                 tzinfo=dt.timezone(dt.timedelta(hours=1))),
#                                         )
#     plot_pressure(p_dict_ds, plotname = "2022-11-01to07_PG60")

#     p_dict_ds = p_dict_date_select(p_dict,
#                         start_date = dt.datetime(2022, 11, 7,0, 0,
#                                 tzinfo=dt.timezone(dt.timedelta(hours=1))), 
#                         end_date = dt.datetime(2022, 11, 7,12, 0,
#                                 tzinfo=dt.timezone(dt.timedelta(hours=1))),
#                                         )
#     plot_pressure(p_dict_ds, plotname = "2022-11-07_PG60")

#     p_dict_ds = p_dict_date_select(p_dict,
#                         start_date = dt.datetime(2022, 11, 2,12, 0,
#                                 tzinfo=dt.timezone(dt.timedelta(hours=1))), 
#                         end_date = dt.datetime(2022, 11, 2,18, 0,
#                                 tzinfo=dt.timezone(dt.timedelta(hours=1))),
#                                         )
#     plot_pressure(p_dict_ds, plotname = "2022-11-02_PG60")

#     #pg90 (on the dls shroud)
#     p_dict = make_p_dict(raw_dict, utc_offset=1, pg_key="pg90_pressure_string")
#     p_dict_ds = p_dict_date_select(p_dict,
#                         start_date = dt.datetime(2022, 11, 1,0, 0,
#                                 tzinfo=dt.timezone(dt.timedelta(hours=1))), 
#                         end_date = dt.datetime(2022, 11, 7,12, 0,
#                                 tzinfo=dt.timezone(dt.timedelta(hours=1))),
#                                         )
#     plot_pressure(p_dict_ds, plotname = "2022-11-01to07_PG90")

#     p_dict_ds = p_dict_date_select(p_dict,
#                         start_date = dt.datetime(2022, 11, 2,12, 0,
#                                 tzinfo=dt.timezone(dt.timedelta(hours=1))), 
#                         end_date = dt.datetime(2022, 11, 2,18, 0,
#                                 tzinfo=dt.timezone(dt.timedelta(hours=1))),
#                                         )
#     plot_pressure(p_dict_ds, plotname = "2022-11-02_PG90")

#     p_dict_ds = p_dict_date_select(p_dict,
#                         start_date = dt.datetime(2022, 11, 2,12, 0,
#                                 tzinfo=dt.timezone(dt.timedelta(hours=1))), 
#                         end_date = dt.datetime(2022, 11, 2,14, 0,
#                                 tzinfo=dt.timezone(dt.timedelta(hours=1))),
#                                         )
#     plot_pressure(p_dict_ds, plotname = "2022-11-02-12_PG90")

    # ########### 2022-09-29
    # file_path = sc_dir + "2022-09-29_PG60.json"
    # raw_dict = load_json(file_path)
    # p_dict = make_p_dict(raw_dict)
    # p_dict_ds = p_dict_date_select(p_dict,
    #                     start_date = dt.datetime(2022, 9, 29,17, 0,
    #                             tzinfo=dt.timezone(dt.timedelta(hours=2))), 
    #                     end_date = dt.datetime(2022, 9, 30,8, 30,
    #                             tzinfo=dt.timezone(dt.timedelta(hours=2))),
    #                                     )
    # plot_pressure(p_dict_ds, plotname = "2022-09-29_PG60")

    # ########### 2022-09-23
    # file_path = sc_dir + "PG_60_2022-09-23.json"
    # raw_dict = load_json(file_path)
    # p_dict = make_p_dict(raw_dict)
    # p_dict_ds = p_dict_date_select(p_dict,
    #                     start_date = dt.datetime(2022, 9, 23,12, 0,
    #                             tzinfo=dt.timezone(dt.timedelta(hours=2))), 
    #                     end_date = dt.datetime(2022, 9, 23,16, 30,
    #                             tzinfo=dt.timezone(dt.timedelta(hours=2))),
    #                                     )
    # plot_pressure(p_dict_ds, plotname = "2022-09-23_PG60")

    # p_dict_ds = p_dict_date_select(p_dict,
    #                     start_date = dt.datetime(2022, 9, 23,13, 0,
    #                             tzinfo=dt.timezone(dt.timedelta(hours=2))), 
    #                     end_date = dt.datetime(2022, 9, 23,14, 0,
    #                             tzinfo=dt.timezone(dt.timedelta(hours=2))),
    #                                     )
    # plot_pressure(p_dict_ds, plotname = "2022-09-23_13-14_PG60")

    # p_dict_ds = p_dict_date_select(p_dict,
    #                     start_date = dt.datetime(2022, 9, 23,12, 0,
    #                             tzinfo=dt.timezone(dt.timedelta(hours=2))), 
    #                     end_date = dt.datetime(2022, 9, 23,13, 0,
    #                             tzinfo=dt.timezone(dt.timedelta(hours=2))),
    #                                     )
    # plot_pressure(p_dict_ds, plotname = "2022-09-23_12-13_PG60")

    # p_dict_ds = p_dict_date_select(p_dict,
    #                     start_date = dt.datetime(2022, 9, 23,14, 0,
    #                             tzinfo=dt.timezone(dt.timedelta(hours=2))), 
    #                     end_date = dt.datetime(2022, 9, 23,14, 30,
    #                             tzinfo=dt.timezone(dt.timedelta(hours=2))),
    #                                     )
    # plot_pressure(p_dict_ds, plotname = "2022-09-23_1sccm_PG60")

    # p_dict_ds = p_dict_date_select(p_dict,
    #                     start_date = dt.datetime(2022, 9, 23,14, 30,
    #                             tzinfo=dt.timezone(dt.timedelta(hours=2))), 
    #                     end_date = dt.datetime(2022, 9, 23,15, 0,
    #                             tzinfo=dt.timezone(dt.timedelta(hours=2))),
    #                                     )
    # plot_pressure(p_dict_ds, plotname = "2022-09-23_2sccm_PG60")

    # p_dict_ds = p_dict_date_select(p_dict,
    #                     start_date = dt.datetime(2022, 9, 23,15, 0,
    #                             tzinfo=dt.timezone(dt.timedelta(hours=2))), 
    #                     end_date = dt.datetime(2022, 9, 23,15, 30,
    #                             tzinfo=dt.timezone(dt.timedelta(hours=2))),
    #                                     )
    # plot_pressure(p_dict_ds, plotname = "2022-09-23_5sccm_PG60")

    # p_dict_ds = p_dict_date_select(p_dict,
    #                     start_date = dt.datetime(2022, 9, 23,15, 30,
    #                             tzinfo=dt.timezone(dt.timedelta(hours=2))), 
    #                     end_date = dt.datetime(2022, 9, 23,16, 0,
    #                             tzinfo=dt.timezone(dt.timedelta(hours=2))),
    #                                     )
    # plot_pressure(p_dict_ds, plotname = "2022-09-23_10sccm_PG60")

    # ###### test
    # file_path = sc_dir + "2022_06_21_PG60_string.json"
    # raw_dict = load_json(file_path)
    # print(raw_dict.keys())
    # print(raw_dict["pg60_pressure_string"][0:3])
    # date_str_format = "%Y-%m-%d %H:%M:%S"
    # # date_test = dt.datetime.strptime(raw_dict["pg60_pressure_string"][0][0],
    # #                                  date_str_format
    # #                                 )
    # date_test = str_to_dt(raw_dict["pg60_pressure_string"][0][0])
    # # print(date_test)
    # # print(type(date_test))
    # p_test = pstr_to_mbar(raw_dict["pg60_pressure_string"][0][1])
    # print(raw_dict["pg60_pressure_string"][0][1])
    # # print(p_test)
    # # print(type(p_test))

    # p_dict = make_p_dict(raw_dict)
    # # print(p_dict)

    # plot_pressure(p_dict, plotname = "2022_06_21_PG60")
    # p_dict_1sccm_9A = p_dict_date_select(p_dict,
    #             start_date = dt.datetime(2022, 6, 21,19, 10,
    #                                 tzinfo=dt.timezone(dt.timedelta(hours=2))), 
    #             end_date = dt.datetime(2022, 6, 21,20,20,
    #                                 tzinfo=dt.timezone(dt.timedelta(hours=2)))
    #                                     )
    # plot_pressure(p_dict_1sccm_9A, plotname = "2022_06_21_PG60_1sccm_9A")
    # p_dict_02sccm_9A = p_dict_date_select(p_dict,
    #             start_date = dt.datetime(2022, 6, 21,17, 50,
    #                                 tzinfo=dt.timezone(dt.timedelta(hours=2))), 
    #             end_date = dt.datetime(2022, 6, 21,19,20,
    #                                 tzinfo=dt.timezone(dt.timedelta(hours=2)))
    #                                     )
    # plot_pressure(p_dict_02sccm_9A, plotname = "2022_06_21_PG60_02sccm_9A")
    

